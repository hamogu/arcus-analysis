from __future__ import print_function
import os
from copy import deepcopy
import numpy as np
import astropy.units as u
from astropy import table
from arcus import tolerances as tol
from marxs.simulator import Sequence
from arcus.defaults import DefaultSource, DefaultPointing
from arcus.arcus import PerfectArcus
from arcus.ralfgrating import CATWindow, CATGratingwithL1
from arcus.spo import ScatterPerChannel
import arcus
from utils import get_path


n_photons = 200000
#n_photons = 20000
src = DefaultSource(energy=0.5)
pnt = DefaultPointing()
wave = np.array([15., 25., 37.]) * u.Angstrom
energies = wave.to(u.keV, equivalencies=u.spectral())

jitter_steps = np.array([0.5, 1., 1.5, 2., 5., 10., 20., 30., 60.]) * u.arcsec
trans_steps = np.array([0., .1, .2, .4, .7, 1., 2., 5., 10.])
rot_steps = np.deg2rad([0., 2., 5., 10., 15., 20., 25., 30., 40., 50., 60., 120., 180.]) / 60.
# Ignore the fact that I'm running many iterations with pure zeros.
# Can be fixed easily in higher numpy version where np.unique( ... axis=) works
n_trans = len(trans_steps)
n_rot = len(rot_steps)
n = 3 * n_trans
changeglobal = np.zeros((n_trans * 2 * 3 + n_rot  * 2 * 3, 6))
for i in range(3):
    changeglobal[i * n_trans: (i + 1) * n_trans, i] = trans_steps
    changeglobal[n + i * n_rot: n + (i + 1) * n_rot, i + 3] = rot_steps


half = changeglobal.shape[0] / 2
changeglobal[half:, :] = - changeglobal[:half, :]
changeindividual = changeglobal[: half, :]

sigma_period = np.logspace(-6, -2, 13) * 0.0002
changeperiod = np.zeros((len(sigma_period), 2))
changeperiod[:, 0] = 0.0002
changeperiod[:, 1] = sigma_period

catflatness = np.deg2rad([[0., .1, .2, .4, .6, .8, 1.]]).T
catbuckeling = np.deg2rad([[0., .25, .5, .75, 1., 1.5, 2., 3., 5]]).T

scatter = np.array([0, .5, 1., 2., 4., 6., 8.])
scatter = np.hstack([np.vstack([scatter, np.zeros_like(scatter)]),
                     np.vstack([np.zeros_like(scatter[1:]), scatter[1:]])])
scatter = np.deg2rad(scatter / 3600.).T


def rename_axes(tab):
    '''Rename axes from x->z, y->x, and z->y in the translation and rotation

    Some elements of Arcus are generated with the xyz2zxy matrix which changes
    the order of the axes. This is the case for example for the SPOs where
    the xyz2zxy matrix is part of the global definition. In other cases,
    e.g. the CATs, the pos4d matrices are generated from the
    Rowland torus, this is instead applied to individual elements.
    Due to the order in which the misalignment matrices are applied, the SPOs
    need to have the columns renamed, while the other parts do not.
    '''
    rename = [('tx', 'temp'), ('ty', 'tx'), ('tz', 'ty'), ('temp', 'tz'),
              ('rx', 'temp'), ('ry', 'rx'), ('rz', 'ry'), ('temp', 'rz')]
    for a, b in rename:
        tab.rename_column(a, b)


def run_for_energies(energies,
                     instrum_before, wigglefunc, wigglepars, instrum_after,
                     outfile,
                     parameters=['tx', 'ty', 'tz', 'rx', 'ry', 'rz'],
                     xyz2zxy=False):
    wave = energies.to(u.Angstrom, equivalencies=u.spectral())
    outtabs = []
    for i, e in enumerate(energies):
        src.energy = e.to(u.keV).value
        photons_in = src.generate_photons(n_photons)
        photons_in = pnt(photons_in)

        out = tol.CaptureResAeff(len(parameters), Ageom=instrum_before.elements[0].area.to(u.cm**2))
        tol.singletolerance(photons_in,
                            instrum_before,
                            wigglefunc,
                            wigglepars,
                            instrum_after,
                            out)
        out.tab['energy'] = e
        out.tab['wave'] = wave[i]
        outtabs.append(out.tab)
    dettab = table.vstack(outtabs)
    if len(parameters) == 1:
        dettab[parameters[0]] = dettab['Parameters'].data
    else:
        for i, c in enumerate(parameters):
            dettab[c] = dettab['Parameters'].data[:, i]
    if xyz2zxy:
        rename_axes(dettab)

    outfull = os.path.join(get_path('tolerances'), outfile)
    dettab.write(outfull, overwrite=True)
    print('Writing {}'.format(outfull))

# jitter
instrum = PerfectArcus(channels='1')
wave = energies.to(u.Angstrom, equivalencies=u.spectral())
outtabs = []
for i, e in enumerate(energies):
    src.energy = e.to(u.keV).value
    photons_in = src.generate_photons(n_photons)

    out = tol.CaptureResAeff(1, Ageom=instrum.elements[0].area.to(u.cm**2))
    for j, jit in enumerate(jitter_steps):
        print('Working on jitter {}/{}'.format(j, len(jitter_steps)))
        jitterpnt = DefaultPointing(jitter=jit)
        p_out = jitterpnt(photons_in.copy())
        p_out = instrum(p_out)
        ind = np.isfinite(p_out['det_x']) & (p_out['probability'] > 0)
        out(jit, p_out[ind], n_photons)
        out.tab['energy'] = e
        out.tab['wave'] = wave[i]
    outtabs.append(out.tab)
dettab = table.vstack(outtabs)
dettab.write(os.path.join(get_path('tolerances'), 'jitter.fits'),
             overwrite=True)

# SPO scatter
instrum = PerfectArcus(channels='1')
run_for_energies(energies=energies,
                 instrum_before=instrum.elements[0],
                 wigglefunc=tol.ScatterVariation(instrum.elements[1],
                                             ScatterPerChannel),
                 wigglepars=scatter,
                 instrum_after=Sequence(elements=instrum.elements[2:]),
                 outfile='scatter.fits',
                 parameters=['inplanescatter', 'perpplanescatter'])


# detectors
instrum = PerfectArcus(channels='1')
run_for_energies(energies=energies,
                 instrum_before=Sequence(elements=instrum.elements[:5]),
                 wigglefunc=tol.WiggleGlobalParallel(instrum.elements[5]),
                 wigglepars=changeglobal,
                 instrum_after=Sequence(elements=instrum.elements[6:]),
                 outfile='detector_global.fits')

instrum = PerfectArcus(channels='1')
run_for_energies(energies=energies,
                 instrum_before=Sequence(elements=instrum.elements[:5]),
                 wigglefunc=tol.WiggleIndividualElements(instrum.elements[5]),
                 wigglepars=changeindividual,
                 instrum_after=Sequence(elements=instrum.elements[6:]),
                 outfile='detector_individual.fits')

# CATs
instrum = PerfectArcus(channels='1')
run_for_energies(energies=energies,
                 instrum_before=Sequence(elements=instrum.elements[:2]),
                 wigglefunc=tol.WiggleGlobalParallel(instrum.elements[2]),
                 wigglepars=changeglobal,
                 instrum_after=Sequence(elements=instrum.elements[3:]),
                 outfile='CAT_global.fits')

instrum = PerfectArcus(channels='1')
# individual CATs are the elements of the CATWindows
run_for_energies(energies=energies,
                 instrum_before=Sequence(elements=instrum.elements[:2]),
                 wigglefunc=tol.WiggleIndividualElements(instrum.elements[2],
                                                         CATWindow),
                 wigglepars=changeindividual,
                 instrum_after=Sequence(elements=instrum.elements[3:]),
                 outfile='CAT_individual.fits')
# Windows are the elements of CATfromMechanical (which is instrum.elements[2])
instrum = PerfectArcus(channels='1')
run_for_energies(energies=energies,
                 instrum_before=Sequence(elements=instrum.elements[:2]),
                 wigglefunc=tol.WiggleIndividualElements(instrum.elements[2]),
                 wigglepars=changeindividual,
                 instrum_after=Sequence(elements=instrum.elements[3:]),
                 outfile='CAT_window.fits')

# Period Variation
instrum = PerfectArcus(channels='1')
run_for_energies(energies=energies,
                 instrum_before=Sequence(elements=instrum.elements[:2]),
                 wigglefunc=tol.PeriodVariation(instrum.elements[2], CATGratingwithL1),
                 wigglepars=changeperiod,
                 instrum_after=Sequence(elements=instrum.elements[3:]),
                 outfile='CAT_period.fits',
                 parameters=['nominal', 'sigma'])

# CAT surfaceflatness
instrum = PerfectArcus(channels='1')
run_for_energies(energies=energies,
                 instrum_before=Sequence(elements=instrum.elements[:2]),
                 wigglefunc=tol.CATFlatnessVariation(instrum.elements[2], CATGratingwithL1),
                 wigglepars=catflatness,
                 instrum_after=Sequence(elements=instrum.elements[3:]),
                 outfile='CAT_flatness.fits',
                 parameters=['sigma'])
# CAT buckeling
instrum = PerfectArcus(channels='1')
run_for_energies(energies=energies,
                 instrum_before=Sequence(elements=instrum.elements[:2]),
                 wigglefunc=tol.CATFlatnessVariation(instrum.elements[2],
                                                     parallel_class=CATGratingwithL1,
                                                     orderselector=tol.OrderSelectorTopHat),
                 wigglepars=catbuckeling,
                 instrum_after=Sequence(elements=instrum.elements[3:]),
                 outfile='CAT_buckeling.fits',
                 parameters=['fullwidth'])


# SPOs
# increase size of aperture to make sure light reaches SPOs.
# In practice thermal precolimators etc.will impose further restrictions.
instrum = PerfectArcus(channels='1')
instrum.elements[0].elements[0].pos4d[0, 1] += np.max(trans_steps)
instrum.elements[0].elements[0].pos4d[1, 2] += np.max(trans_steps)

run_for_energies(energies=energies,
                 instrum_before=Sequence(elements=instrum.elements[:1]),
                 wigglefunc=tol.WiggleGlobalParallel(instrum.elements[1]),
                 wigglepars=changeglobal,
                 instrum_after=Sequence(elements=instrum.elements[2:]),
                 outfile='SPOs_global.fits')

instrum = PerfectArcus(channels='1')
instrum.elements[0].elements[0].pos4d[0, 1] += np.max(trans_steps)
instrum.elements[0].elements[0].pos4d[1, 2] += np.max(trans_steps)
run_for_energies(energies=energies,
                 instrum_before=Sequence(elements=instrum.elements[:1]),
                 wigglefunc=tol.WiggleIndividualElements(instrum.elements[1]),
                 wigglepars=changeindividual,
                 instrum_after=Sequence(elements=instrum.elements[2:]),
                 outfile='SPOs_individual.fits')


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
catblazegradients = np.deg2rad(np.array([-2, -1.5, -1, -.5, -.17, 0., .17, .5, 1., 1.5, 2.]) / 30)

outtabs = []
for i, e in enumerate(energies):
    src.energy = e.to(u.keV).value
    photons_in = src.generate_photons(n_photons)
    photons_in = pnt(photons_in)
    out = tol.CaptureResAeff(1, Ageom=instrum.elements[0].area.to(u.cm**2))
    for j, grad in enumerate(catblazegradients):
        print('Working on CAT gradient {}/{}'.format(j, len(catblazegradients)))
        rg.CATWindow.extra_elem_args['d_blaze_mm'] = grad
        instrum = PerfectArcus(channels='1')
        p_out = instrum(photons_in.copy())
        ind = np.isfinite(p_out['det_x']) & (p_out['probability'] > 0)
        out(grad, p_out[ind], n_photons)
        out.tab['energy'] = e
        out.tab['wave'] = wave[i]
    outtabs.append(out.tab)
dettab = table.vstack(outtabs)
dettab.write(os.path.join(get_path('tolerances'), 'blazegradient.fits'),
             overwrite=True)


# More detailed analysis of Ralf's worst case gratings
rg.CATWindow.elem_class = rg.GeneralLinearNonParallelCAT
# add a rotation to correct for the blaze problem
blazecorrect = np.zeros((13, 6))
blazecorrect[:, 4] = np.deg2rad(np.arange(-1.2, 1.3, .2))

for i, blazeargs in enumerate([(0.036, -0.8), (0.036, 0.8),
                               (-0.036, -0.8), (-0.036, 0.8)]):
    rg.CATWindow.extra_elem_args['d_blaze_mm'] = np.deg2rad(blazeargs[0])
    rg.CATWindow.extra_elem_args['blaze_center'] = np.deg2rad(blazeargs[1])
    instrum = PerfectArcus(channels='1')
    run_for_energies(energies=energies,
                     instrum_before=Sequence(elements=instrum.elements[:2]),
                     wigglefunc=tol.MoveIndividualParallel(instrum.elements[2],
                                                           CATWindow),
                     wigglepars=blazecorrect,
                     instrum_after=Sequence(elements=instrum.elements[3:]),
                     outfile='CAT_blaze_detail{}.fits'.format(i))
