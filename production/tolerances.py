import os
import argparse
import numpy as np
import astropy.units as u
from astropy.table import Table
from marxs.design.tolerancing import (wiggle, moveglobal, moveindividual,
                                      varyattribute, varyperiod,
                                      generate_6d_wigglelist)
from marxs.analysis.gratings import CaptureResAeff_CCDgaps
from marxs.missions.arcus import arcus
from marxs.missions.arcus.defaults import DefaultPointing, DefaultSource
from marxs.missions.athena.spo import ScatterPerChannel, SPOChannelMirror
from marxs.missions.arcus import ralfgrating
from marxs.missions.mitsnl.catgrating import NonParallelCATGrating
import marxs.missions.mitsnl.tolerances as tol
from marxs.optics import aperture
from marxs.design.tolerancing import run_tolerances_for_energies2
from marxs.missions.arcus.utils import config

from utils import get_path

conf = arcus.defaultconf


parser = argparse.ArgumentParser(description='Run tolerancing simulations.')
parser.add_argument('--n_photons', default=200000, type=int,
                    help='Number of photons per simulation')
parser.add_argument('--n_baseline_budget', default=100, type=int,
                    help='Number of baseline mislignment simulations to run ' +
                    '(This parameter has no effect is baseline_budget is not in the list of scenarios.)')
parser.add_argument('-s',  '--scenario', action="extend", nargs="+", type=str,
                    help='Specify which scenarios to to run. If argument is not set, all test will be run.')
parser.add_argument('-e', '--exclude', action="extend", nargs="+", type=str,
                    help='Exclude specific scenarios from running')
args = parser.parse_args()


src = DefaultSource(energy=0.5)
wave = np.array([15., 25., 37.]) * u.Angstrom
energies = wave.to(u.keV, equivalencies=u.spectral())

def run_n_errorbudgets(align, conf, n=50, n_photons=200_000):
    '''Run Arcus R/Aeff simulations for a particular set of alignment tolerances

    Parameters
    ----------
    align : list
        Error budget in the Randall-Smith form
    conf : dict
        Configuration dictionary. In particular, the alignment tolerance table in
        that dict determines the aligments for the runs.
    n : int
        Number of simulations. The first simulation is always run with a perfect
        instrument for comparison, so ``n=50`` will get 49 simulations with random
        misalignments.
    n_photons : int
        Number of photons for each simulation.

    Returns
    -------
    tab : `astropy.table.Table`
        Table with results. The first row is a run with a perfect instrument for
        comparison, the remaining rows are runs with random realizations of the
        alignment tolerances in `conf`.
    '''
    out = []

    for i in range(n):
        print('Run tolerance budget: {}/{}'.format(i, n))

        arcus.reformat_randall_errorbudget(align, globalfac=None)
        conf['alignmentbudget'] = align

        if i == 0:
            arc = arcus.PerfectArcus(channels='1')
        else:
            arc = arcus.Arcus(channels='1', conf=conf)

        for e in energies:
            src.energy = e
            photons_in = src.generate_photons(n_photons * u.s)
            photons_in = DefaultPointing()(photons_in)
            photons = arc(photons_in)

            out.append(analyzer(photons))
            out[-1]['energy'] = e.to(u.keV, equivalencies=u.spectral()).value
            out[-1]['run'] = i

    tab = Table([{d: out[i][d].value
                  if isinstance(out[i][d], u.Quantity) else out[i][d]
                  for d in out[i]} for i in range(len(out))])

    tab['energy'].unit = u.keV
    tab['wave'] = tab['energy'].to(u.Angstrom, equivalencies=u.spectral())
    return tab

translation_list = [0., .1, .2, .4, .7, 1., 2., 5., 10.] * u.mm
rotation_list = [0., 2., 5., 10., 15., 20., 25., 30., 40., 50., 60., 120., 180.] * u.arcmin
changeglobal, changeindividual = generate_6d_wigglelist(translation_list,
                                                        rotation_list)


scatter = np.array([0, .5, 1., 2., 4., 6., 8.])
scatter = np.hstack([np.vstack([scatter, np.zeros_like(scatter)]),
                     np.vstack([np.zeros_like(scatter[1:]), scatter[1:]])])
scatter = np.deg2rad(scatter / 3600.).T

instrumfull = arcus.PerfectArcus(channels='1')

# In general, we measure the resolving power from an ideal, circular detector. That's very, very
# close to detectors on the rowland circle, and it saves us from dealing with chip gaps.
# However, if we actually want to wiggle the detectors, then we actually need to measure
# the results on the detectors, thus the specific analyzer for that case.
analyzer = CaptureResAeff_CCDgaps(A_geom=instrumfull.elements[0].area.to(u.cm**2),
                          dispersion_coord='circ_phi',
                          orders=np.arange(-12, 5),
                          aeff_filter_col='CCD')
analyzer_det = CaptureResAeff_CCDgaps(A_geom=instrumfull.elements[0].area.to(u.cm**2),
                          dispersion_coord='proj_x',
                          orders=np.arange(-12, 5),
                          aeff_filter_col='CCD')

def filter_noCCD(photons):
    photons['probability'][~np.isfinite(photons['det_x'])] = 0
    return photons

class Arcus(arcus.PerfectArcus):
    def post_process(self):
        return []

    def __init__(self):
        super().__init__(channels='1')
        self.elements.insert(0, DefaultPointing())
        self.elements.append(filter_noCCD)

arcus_eff_tab = Table.read(os.path.join(config['data']['caldb_inputdata'],
                                        'gratings', 'efficiency.csv'), format='ascii.ecsv')

def increase_aperture_size(instrum, pars):
    '''
    Increase size of aperture to make sure light reaches SPOs.
    In practice thermal precolimators etc.will impose further restrictions.
    '''
    max_size_increase = np.max(translation_list).to(u.mm).value
    for elem in instrum.elements_of_class(aperture.RectangleAperture):
        elem.pos4d[0, 1] += max_size_increase
        elem.pos4d[1, 2] += max_size_increase

runs = {'jitter': (DefaultPointing, analyzer, varyattribute,
                  [{'jitter': j} for j in np.array([0.1, 0.2, 0.25, 0.5, 1., 1.5, 2., 5., 10., 20.]) * u.arcsec]),
        'scatter': (ScatterPerChannel, analyzer, varyattribute,
                   [{'inplanescatter': a, 'perpplanescatter':b} for a, b in scatter]),
        'detector_global': (arcus.DetCamera, analyzer_det, moveglobal, changeglobal),
        'detector_individual': (arcus.DetCamera, analyzer_det, moveindividual, changeglobal),
        'CAT_global': (ralfgrating.CATfromMechanical, analyzer, moveglobal, changeglobal),
        'CAT_window': (ralfgrating.CATfromMechanical, analyzer, wiggle, changeindividual),
        'CAT_individual': (ralfgrating.CATWindow, analyzer, wiggle, changeindividual),
        'CAT_period': (NonParallelCATGrating, analyzer, varyperiod,
                       [{'period_mean': 0.0002, 'period_sigma': s} for s in np.logspace(-6, -2, 13) * 0.0002]),
        'CAT_flatness': (NonParallelCATGrating, analyzer, varyattribute,
                         [{'order_selector': tol.OrderSelectorWavy(wavysigma=s, tab=arcus_eff_tab)} for s in np.deg2rad([0., .1, .2, .4, .6, .8, 1.])]),
        'CAT_buckling': (NonParallelCATGrating, analyzer, varyattribute,
                         [{'order_selector': tol.OrderSelectorTopHat(tophatwidth=s, tab=arcus_eff_tab)} for s in np.deg2rad([0., .25, .5, .75, 1., 1.5, 2., 3., 5])]),
        'blazegradient': (NonParallelCATGrating, analyzer, varyattribute,
                          [{'d_blaze_mm': s} for s in np.deg2rad(np.array([-2, -1.5, -1, -.5, -.17, 0., .17, .5, 1., 1.5, 2.]) / 30)]),
        'SPOs_global': (SPOChannelMirror, analyzer, moveglobal, changeglobal, increase_aperture_size),
        'SPOs_individual': (SPOChannelMirror, analyzer, wiggle, changeindividual, increase_aperture_size),
        'baseline_budget': (None),
            }
'''
Format for the runs: dict with entries:  name: (element, wigglefunc, parameters, preparefunc)
name : string -
    Name of run, also use as filename
element : marx simulation elements
analyser : Marxs analyser funcstion
wigglefunc : function
parameters : list
preparefunc (optional) : function
    Will be executed before the test run, use this to modify the instrument in preparation
'''

if args.scenario is None:
    scenarios = runs.keys()
else:
    scenarios = args.scenario

if args.exclude is not None:
    scenarios = [e for e in scenarios if e not in args.exclude]

print('Running the following scenarios:', scenarios)

for outfile in scenarios:
    pars = runs[outfile]
    instrum = Arcus()
    if outfile == 'baseline_budget':
        align = arcus.align_requirement_smith
        tab = run_n_errorbudgets(align, conf, n=args.n_baseline_budget,
                                 n_photons=args.n_photons)
    else:
        if len(pars) > 4:
            pars[4](instrum, pars)
        tab = run_tolerances_for_energies2(src, energies, instrum,
                                           pars[0], pars[2], pars[3],
                                           pars[1],
                                           t_source=args.n_photons * u.s)

    # For column with dtype object
    # This happens only when the input is the orderselector, so we can special
    # special case that here
    if 'order_selector' in tab.colnames:
        o0 = tab['order_selector'][0]
        if hasattr(o0, 'sigma'):
            tab['sigma'] = [o.sigma for o in tab['order_selector']]
        elif hasattr(o0, 'tophatwidth'):
            tab['tophatwidth'] = [o.tophatwidth for o in tab['order_selector']]
        tab.remove_column('order_selector')

    outfull = os.path.join(get_path('tolerances'), outfile + '.fits')
    tab.write(outfull, overwrite=True)
    print('Writing {}'.format(outfull))
