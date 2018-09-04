from __future__ import print_function
import os
from copy import deepcopy
import numpy as np
import astropy.units as u
from astropy import table
from arcus import tolerances as tol
from marxs.simulator import Sequence
from marxs.design.tolerancing import (wiggle, moveglobal, moveindividual,
                                      varyattribute, varyorderselector, CaptureResAeff,
                                      run_tolerances)
from arcus.defaults import DefaultSource, DefaultPointing
from arcus.arcus import PerfectArcus
from arcus.ralfgrating import CATWindow, CATGratingwithL1
from arcus.spo import ScatterPerChannel
import arcus
from utils import get_path


n_photons = 200000
#n_photons = 20000
src = DefaultSource(energy=0.5)
wave = np.array([15., 25., 37.]) * u.Angstrom
energies = wave.to(u.keV, equivalencies=u.spectral())

trans_steps = np.array([0., .1, .2, .4, .7, 1., 2., 5., 10.])
rot_steps = np.deg2rad([0., 2., 5., 10., 15., 20., 25., 30., 40., 50., 60., 120., 180.]) / 60.
n_trans = len(trans_steps)
n_rot = len(rot_steps)
n = 3 * n_trans
changeglobal = np.zeros((n_trans * 2 * 3 + n_rot  * 2 * 3, 6))
for i in range(3):
    changeglobal[i * n_trans: (i + 1) * n_trans, i] = trans_steps
    changeglobal[n + i * n_rot: n + (i + 1) * n_rot, i + 3] = rot_steps


half = changeglobal.shape[0] / 2
changeglobal[int(half):, :] = - changeglobal[:int(half), :]
changeindividual = changeglobal[: int(half), :]
# Remove multiple entries with [0,0,0,0, ...]
changeglobal = np.unique(changeglobal, axis=0)
changeindividual = np.unique(changeindividual, axis=0)



scatter = np.array([0, .5, 1., 2., 4., 6., 8.])
scatter = np.hstack([np.vstack([scatter, np.zeros_like(scatter)]),
                     np.vstack([np.zeros_like(scatter[1:]), scatter[1:]])])
scatter = np.deg2rad(scatter / 3600.).T

instrumfull = PerfectArcus(channels='1')
analyser = tol.CaptureResAeff(Ageom=instrum.elements[0].area.to(u.cm**2),
                              dispersion_coord='proj_x', orders=np.arange(-12, 5))


def run_for_energies(instrum_before, wigglefunc, wiggleparts, parameters,
                     outfull, instrum):
    outtabs = []
    for i, e in enumerate(energies):
        src.energy = e.to(u.keV).value
        photons_in = src.generate_photons(n_photons)
        photons_in = instrum_before(photons_in)
        tab = table.Table(run_tolerances(photons_in, instrum, wigglefunc, wiggleparts,
                                         parameters, analyzer))
        tab['energy'] = e
        tab['wave'] = wave[i]
        outtabs.append(tab)
    dettab = table.vstack(outtabs)
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

# jitter
def dummy(p):
    '''Function needs something here, but nothing happens'''
    return p

instrum = Arcus()
run_for_energies(dummy, varyattribute, instrum.elements[0],
                 [{'jitter': j} for j in np.array([0.5, 1., 1.5, 2., 5., 10., 20., 30., 60.]) * u.arcsec],
                 Sequence(elements=instrum.elements[1:]),
                 'jitter.fits')

# SPO scatter
instrum = Arcus()
run_for_energies(Sequence(elements=instrum.elements[:2]), varyattribute,
                 instrum.elements[2].elements,
                 [{'inplanescatter': a, 'perpplanescatter':b} for a, b in scatter],
                 Sequence(elements=instrum.elements[2:]),
                 'scatter.fits')

# detectors
def tabtolist(tab6):
    l = [{'dx': a, 'dy': b, 'dz': c, 'rx': d, 'ry': e, 'rz': f}
         for (a,b, c, d, e, f) in tab6]
    return l

instrum = Arcus()
run_for_energies(instrum.elements[:6], moveglobal,
                 instrum.elements[6],
                 tab2list(changeglobal),
                 Sequence(elements, instrum.elements[6:]),
                 'detector_global.fits')


instrum = Arcus()
run_for_energies(instrum.elements[:6], moveindividual,
                 instrum.elements[6],
                 tab2list(changeglobal),
                 Sequence(elements, instrum.elements[6:]),
                 'detector_individual.fits')

# CATs
instrum = Arcus()
run_for_energies(instrum.elements[:3], moveglobal,
                 instrum.elements[3],
                 tab2list(changeglobal),
                 Sequence(elements, instrum.elements[3:]),
                 'CAT_global.fits')

# individual CATs are the elements of the CATWindows
instrum = Arcus()
run_for_energies(instrum.elements[:3], wiggle,
                 instrum.elements_of_class(CATWindow),
                 tab2list(changeindividual),
                 Sequence(elements, isstrum.elements[3:]),
                 'CAT_individual.fits')
# Windows are the elements of CATfromMechanical (which is instrum.elements[3])
instrum = Arcus()
run_for_energies(instrum.elements[:3], wiggle,
                 instrum.elements[3],
                 tab2list(changeindividual),
                 Sequence(elements, instrum.elements[3:]),
                 'CAT_window.fits')

# Period Variation
instrum = Arcus()
run_for_energies(instrum.elements[:3], varyperiod,
                 instrum.elements_of_class(CATGratingwithL1),
                 [{'period_mean': 0.0002, 'period_sigma': s} for s in np.logspace(-6, -2, 13) * 0.0002],
                 Sequence(elements, instrum.elements[3:]),
                 'CAT_period.fits')


# CAT surfaceflatness
instrum = Arcus()
run_for_energies(instrum.elements[:3], varyperiod,
                 instrum.elements_of_class(CATGratingwithL1),
                 [{'orderselector': tol.OrderSelectorWavy, 'wavysigma': s} for s in np.deg2rad([0., .1, .2, .4, .6, .8, 1.])],
                 Sequence(elements, instrum.elements[3:]),
                 'CAT_flatness.fits')

# CAT buckeling
run_for_energies(instrum.elements[:3], varyperiod,
                 instrum.elements_of_class(CATGratingwithL1),
                 [{'orderselector': tol.OrderSelectorTopHat, 'tophatwidth': s} for s in np.deg2rad([0., .25, .5, .75, 1., 1.5, 2., 3., 5])],
                 Sequence(elements, instrum.elements[3:]),
                 'CAT_buckeling.fits')

# SPOs
# increase size of aperture to make sure light reaches SPOs.
# In practice thermal precolimators etc.will impose further restrictions.
instrum = Arcus()
instrum.elements[1].elements[0].pos4d[0, 1] += np.max(trans_steps)
instrum.elements[1].elements[0].pos4d[1, 2] += np.max(trans_steps)

run_for_energies(instrum.elements[:2], moveglobal,
                 instrum.elements[2],
                 tab2list(changeglobal),
                 Sequence(elements, instrum.elements[2:]),
                 'SPOs_global.fits')

moveglobal(instrum.elements[2], 0, 0, 0, 0, 0, 0)

run_for_energies(instrum.elements[:2], wiggle,
                 instrum.elements[2],
                 tab2list(changeindividual),
                 Sequence(elements, instrum.elements[2:]),
                 'SPOs_individual.fits')

# Run default tolerance budget a few times
n_budget = 50
out = tol.CaptureResAeff(2, Ageom=instrum.elements[0].area.to(u.cm**2))

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
        photons_in = src.generate_photons(n_photons)
        photons_in = pnt(photons_in)
        photons = arc(photons_in)
        good = (photons['probability'] > 0) & (photons['CCD'] > 0)
        out([i, src.energy], photons[good], n_photons)

out.tab['energy'] = out.tab['Parameters'].data[:, 1] * u.keV
out.tab['wave'] = out.tab['energy'].to(u.Angstrom, equivalencies=u.spectral())
out.tab['run'] = out.tab['Parameters'].data[:, 0]

outfull = os.path.join(get_path('tolerances'), 'baseline_budget.fits')
out.tab.write(outfull, overwrite=True)
print('Writing {}'.format(outfull))


# This one should be part of step 1, but needs to be run at the end of the
# script because it monkey-patches the definition of Arcus and we
# do not want to mess up any of the other runs.
instrum = PerfectArcus(channels='1')  # Just to get Aeff below
import arcus.ralfgrating as rg
rg.CATWindow.elem_class = rg.NonParallelCATGrating
rg.CATWindow.extra_elem_args['d_blaze_mm'] = 1. # any dummy value

instrum = Arcus()
run_for_energies(instrum.elements[:3], varyattribute,
                 instrum.elements_of_class(CATGratingwithL1),
                 [{'d_blaze_mm': s} for s in np.deg2rad(np.array([-2, -1.5, -1, -.5, -.17, 0., .17, .5, 1., 1.5, 2.]) / 30)],
                 Sequence(elements, instrum.elements[3:]),
                 'blazegradient.fits')

# More detailed analysis of Ralf's worst case gratings
rg.CATWindow.elem_class = rg.GeneralLinearNonParallelCAT
# add a rotation to correct for the blaze problem

for i, blazeargs in enumerate([(0.036, -0.8), (0.036, 0.8),
                               (-0.036, -0.8), (-0.036, 0.8)]):
    rg.CATWindow.extra_elem_args['d_blaze_mm'] = np.deg2rad(blazeargs[0])
    rg.CATWindow.extra_elem_args['blaze_center'] = np.deg2rad(blazeargs[1])
    instrum = Arcus()
    run_for_energies(instrum.elements[:3], moveindividual,
                 instrum.elements_of_class(CATGratingwithL1),
                 [{'ry': s} for s innp.deg2rad(np.arange(-1.2, 1.3, .2)) ],
                 Sequence(elements, instrum.elements[3:]),
                 f'CAT_blaze_detail{i}.fits')
