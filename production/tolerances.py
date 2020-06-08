import os
import argparse
from copy import deepcopy
import numpy as np
import astropy.units as u
from astropy import table
from arcus import tolerances
from marxs.simulator import Sequence
from marxs.design.tolerancing import (wiggle, moveglobal, moveindividual,
                                      varyattribute, varyorderselector,
                                      varyperiod,
                                      CaptureResAeff,
                                      run_tolerances,
                                      generate_6d_wigglelist,
                                      run_tolerances_for_energies)
from marxs.optics import CATGrating
from arcus.defaults import DefaultSource, DefaultPointing
from arcus.arcus import PerfectArcus
from arcus.ralfgrating import CATWindow
from arcus.spo import ScatterPerChannel
import arcus.tolerances as tol
import arcus
from utils import get_path

parser = argparse.ArgumentParser(description='Run tolerancing simulations.')
parser.add_argument('--n_photons', default=200000, type=int)
args = parser.parse_args()


src = DefaultSource(energy=0.5)
wave = np.array([15., 25., 37.]) * u.Angstrom
energies = wave.to(u.keV, equivalencies=u.spectral())

changeglobal, changeindividual = generate_6d_wigglelist([0., .1, .2, .4, .7, 1., 2., 5., 10.] * u.mm,
                                                        [0., 2., 5., 10., 15., 20., 25., 30., 40., 50., 60., 120., 180.] * u.arcmin)


scatter = np.array([0, .5, 1., 2., 4., 6., 8.])
scatter = np.hstack([np.vstack([scatter, np.zeros_like(scatter)]),
                     np.vstack([np.zeros_like(scatter[1:]), scatter[1:]])])
scatter = np.deg2rad(scatter / 3600.).T

instrumfull = PerfectArcus(channels='1')
analyzer = CaptureResAeff(A_geom=instrumfull.elements[0].area.to(u.cm**2),
                          dispersion_coord='proj_x',
                          orders=np.arange(-12, 5))


def run_for_energies(instrum_before, wigglefunc, wiggleparts, parameters,
                     instrum, outfile, reset=None):
    dettab = run_tolerances_for_energies(src, energies,
                                         instrum_before, instrum,
                                         wigglefunc, wiggleparts, parameters,
                                         analyzer, reset=reset,
                                         t_source=args.n_photons)
    # For column with dtype object
    # This happens only when the input is the orderselector, so we can special
    # special case that here
    if 'order_selector' in dettab.colnames:
        o0 = dettab['order_selector'][0]
        if hasattr(o0, 'sigma'):
            dettab['sigma'] = [o.sigma for o in dettab['order_selector']]
        elif hasattr(o0, 'tophatwidth'):
            dettab['tophatwidth'] = [o.tophatwidth for o in dettab['order_selector']]
        dettab.remove_column('order_selector')
    outfull = os.path.join(get_path('tolerances'), outfile)
    dettab.write(outfull, overwrite=True)
    print('Writing {}'.format(outfull))


def filter_noCCD(photons):
    photons['probability'][~np.isfinite(photons['det_x'])] = 0
    return photons


class Arcus(PerfectArcus):
    def post_process(self):
        return []

    def __init__(self):
        super().__init__(channels='1')
        self.elements.insert(0, DefaultPointing())
        self.elements.append(filter_noCCD)

reset_6d = {'dx': 0., 'dy': 0., 'dz': 0., 'rx': 0., 'ry': 0., 'rz': 0.}


# # jitter
# def dummy(p):
#     '''Function needs something here, but nothing happens'''
#     return p

# instrum = Arcus()
# run_for_energies(dummy, varyattribute, instrum.elements[0],
#                  [{'jitter': j} for j in np.array([0.5, 1., 1.5, 2., 5., 10., 20., 30., 60.]) * u.arcsec],
#                  instrum,
#                  'jitter.fits')

# # SPO scatter
# instrum = Arcus()
# run_for_energies(Sequence(elements=instrum.elements[:2]), varyattribute,
#                  instrum.elements[2].elements[1],
#                  [{'inplanescatter': a, 'perpplanescatter':b} for a, b in scatter],
#                  Sequence(elements=instrum.elements[2:]),
#                  'scatter.fits')

# # detectors
# instrum = Arcus()
# run_for_energies(Sequence(elements=instrum.elements[:9]), moveglobal,
#                  instrum.elements[9],
#                  changeglobal,
#                  Sequence(elements=instrum.elements[9:]),
#                  'detector_global.fits')


# instrum = Arcus()
# run_for_energies(Sequence(elements=instrum.elements[:9]), moveindividual,
#                  instrum.elements[9],
#                  changeglobal,
#                  Sequence(elements=instrum.elements[9:]),
#                  'detector_individual.fits')

# # CATs
# instrum = Arcus()
# run_for_energies(Sequence(elements=instrum.elements[:3]), moveglobal,
#                  instrum.elements[3].elements[0],
#                  changeglobal,
#                  Sequence(elements=instrum.elements[3:]),
#                  'CAT_global.fits')

# # individual CATs are the elements of the CATWindows
# instrum = Arcus()
# run_for_energies(Sequence(elements=instrum.elements[:3]), wiggle,
#                  instrum.elements_of_class(CATWindow),
#                  changeindividual,
#                  Sequence(elements=instrum.elements[3:]),
#                  'CAT_individual.fits')
# # Windows are the elements of CATfromMechanical (which is instrum.elements[3])
# instrum = Arcus()
# run_for_energies(Sequence(elements=instrum.elements[:3]), wiggle,
#                  instrum.elements[3].elements[0],
#                  changeindividual,
#                  Sequence(elements=instrum.elements[3:]),
#                  'CAT_window.fits')

# # Period Variation
# instrum = Arcus()
# run_for_energies(Sequence(elements=instrum.elements[:3]), varyperiod,
#                  instrum.elements_of_class(CATGrating),
#                  [{'period_mean': 0.0002, 'period_sigma': s} for s in np.logspace(-6, -2, 13) * 0.0002],
#                  Sequence(elements=instrum.elements[3:]),
#                  'CAT_period.fits')


# # CAT surfaceflatness
# instrum = Arcus()
# run_for_energies(Sequence(elements=instrum.elements[:3]), varyattribute,
#                  instrum.elements_of_class(CATGrating),
#                  [{'order_selector': tol.OrderSelectorWavy(wavysigma=s)} for s in np.deg2rad([0., .1, .2, .4, .6, .8, 1.])],
#                  Sequence(elements=instrum.elements[3:]),
#                  'CAT_flatness.fits')

# # CAT buckeling
# run_for_energies(Sequence(elements=instrum.elements[:3]), varyattribute,
#                  instrum.elements_of_class(CATGrating),
#                  [{'order_selector': tol.OrderSelectorTopHat(tophatwidth=s)} for s in np.deg2rad([0., .25, .5, .75, 1., 1.5, 2., 3., 5])],
#                  Sequence(elements=instrum.elements[3:]),
#                  'CAT_buckeling.fits')

# # SPOs
# # increase size of aperture to make sure light reaches SPOs.
# # In practice thermal precolimators etc.will impose further restrictions.
# instrum = Arcus()
# instrum.elements[1].elements[0].pos4d[0, 1] += np.max(trans_steps)
# instrum.elements[1].elements[0].pos4d[1, 2] += np.max(trans_steps)

# run_for_energies(Sequence(elements=instrum.elements[:2]), moveglobal,
#                  instrum.elements[2].elements[0],
#                  changeglobal,
#                  Sequence(elements=instrum.elements[2:]),
#                  'SPOs_global.fits')


# run_for_energies(Sequence(elements=instrum.elements[:2]), wiggle,
#                  instrum.elements[2].elements[0],
#                  changeindividual,
#                  Sequence(elements=instrum.elements[2:]),
#                  'SPOs_individual.fits')

# Run default tolerance budget a few times
n_budget = 50
out = []

conf = deepcopy(arcus.arcus.defaultconf)

for i in range(n_budget):
    print('Run default tolerance budget: {}/{}'.format(i, n_budget))
    align = deepcopy(arcus.arcus.align_requirement_smith)
    arcus.arcus.reformat_randall_errorbudget(align, globalfac=None)
    conf['alignmentbudget'] = align
    if i == 0:
        arc = PerfectArcus(channels='1')
    else:
        arc = arcus.arcus.Arcus(channels='1', conf=conf)

    for e in energies:
        src.energy = e.to(u.keV).value
        photons_in = src.generate_photons(args.n_photons)
        photons_in = DefaultPointing()(photons_in)
        photons = arc(photons_in)
        # good = (photons['probability'] > 0) & (photons['CCD'] > 0)
        # out([i, src.energy], photons[good], n_photons)
        out.append(analyzer(photons))
        out[-1]['energy'] = e.value
        out[-1]['run'] = i

tab = table.Table([{d: out[i][d].value
                    if isinstance(out[i][d], u.Quantity) else out[i][d]
                    for d in out[i]} for i in range(len(out))])

tab['energy'].unit = u.keV
tab['wave'] = tab['energy'].to(u.Angstrom, equivalencies=u.spectral())

outfull = os.path.join(get_path('tolerances'), 'baseline_budget.fits')
tab.write(outfull, overwrite=True)
print('Writing {}'.format(outfull))


# This one should be part of step 1, but needs to be run at the end of the
# script because it monkey-patches the definition of Arcus and we
# do not want to mess up any of the other runs.
instrum = PerfectArcus(channels='1')  # Just to get Aeff below
import arcus.ralfgrating as rg

class CATL1L2Stack(rg.FlatStack):
    elements = [rg.GeneralLinearNonParallelCAT,
                rg.CATGratingL1,
                rg.L2,
                rg.RandomGaussianScatter]
    keywords = [{'order_selector': rg.globalorderselector,
                 'd': 0.0002,
                 'd_blaze_mm': 1,
                 'blaze_center': 0.},
                {'d': 0.005,
                 'order_selector': rg.l1orderselector,
                 'groove_angle': np.pi / 2.},
                {},
                {'scatter': rg.l2diffraction}]
    def __init__(self, **kwargs):
        kwargs['elements'] = self.elements
        kwargs['keywords'] = self.keywords
        super().__init__(**kwargs)


class CATWindow(rg.Parallel):

    id_col = 'facet'

    def __init__(self, **kwargs):
        kwargs['id_col'] = self.id_col
        kwargs['elem_class'] = CATL1L2Stack
        super().__init__(**kwargs)


rg.CATWindow = CATWindow

instrum = Arcus()
run_for_energies(Sequence(elements=instrum.elements[:3]), varyattribute,
                 instrum.elements_of_class(rg.GeneralLinearNonParallelCAT),
                 [{'d_blaze_mm': s} for s in np.deg2rad(np.array([-2, -1.5, -1, -.5, -.17, 0., .17, .5, 1., 1.5, 2.]) / 30)],
                 Sequence(elements=instrum.elements[3:]),
                 'blazegradient.fits')
