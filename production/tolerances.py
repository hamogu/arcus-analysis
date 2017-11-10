from __future__ import print_function
import os
import numpy as np
import astropy.units as u
from astropy import table
from arcus import tolerances as tol
from marxs.simulator import Sequence
from arcus.defaults import DefaultSource, DefaultPointing
from arcus.arcus import PerfectArcus
from utils import get_path


n_photons = 20000
src = DefaultSource(energy=0.5)
pnt = DefaultPointing()
wave = np.array([15., 25., 35.]) * u.Angstrom
energies = wave.to(u.keV, equivalencies=u.spectral())

trans_steps = np.array([0., .1, .5, 1., 2., 5., 10., 20., 50.])
rot_steps = np.deg2rad([0., 0.1, 0.25, 0.5, 1., 2., 5., 10., 30., 60., 120., 180.]) / 60.
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


instrum = PerfectArcus(channels='1')


def rename_axes(tab):
    '''Rename axes from x->z, y->x, and z->y in the translation and rotation

    Some elements of Arcus are generated with the xyz2zxy matrix which changes
    The order of the axes. This is the case for example for the SPOs where
    the xyz2zxy matrix is part of the global definition. Inother cases,
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
                     outfile, xyz2zxy=False):
    wave = energies.to(u.Angstrom, equivalencies=u.spectral())
    outtabs = []
    for i, e in enumerate(energies):
        src.energy = e.to(u.keV).value
        photons_in = src.generate_photons(n_photons)
        photons_in = pnt(photons_in)

        out = tol.CaptureResAeff(Ageom=instrum_before.elements[0].area.to(u.cm**2))
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
    for i, c in enumerate(['tx', 'ty', 'tz', 'rx', 'ry', 'rz']):
        dettab[c] = dettab['Parameters'].data[:, i]
    if xyz2zxy:
        rename_axes(dettab)
    dettab.write(os.path.join(get_path('tolerances'), outfile),
                 overwrite=True)

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
pos4d = np.eye(4)
pos4d[:, 3] = np.array(instrum.elements[2].elements[0].elem_pos).mean(axis=0)[:, 3]
instrum.elements[2].elements[0].move_center(pos4d)
run_for_energies(energies=energies,
                 instrum_before=Sequence(elements=instrum.elements[:2]),
                 wigglefunc=tol.WiggleGlobalParallel(instrum.elements[2]),
                 wigglepars=changeglobal,
                 instrum_after=Sequence(elements=instrum.elements[3:]),
                 outfile='CAT_global.fits')

instrum = PerfectArcus(channels='1')
run_for_energies(energies=energies,
                 instrum_before=Sequence(elements=instrum.elements[:2]),
                 wigglefunc=tol.WiggleIndividualElements(instrum.elements[2]),
                 wigglepars=changeindividual,
                 instrum_after=Sequence(elements=instrum.elements[3:]),
                 outfile='CAT_individual.fits')

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
                 outfile='SPOs_global.fits',
                 xyz2zxy=True)

instrum = PerfectArcus(channels='1')
instrum.elements[0].elements[0].pos4d[0, 1] += np.max(trans_steps)
instrum.elements[0].elements[0].pos4d[1, 2] += np.max(trans_steps)
run_for_energies(energies=energies,
                 instrum_before=Sequence(elements=instrum.elements[:1]),
                 wigglefunc=tol.WiggleIndividualElements(instrum.elements[1]),
                 wigglepars=changeindividual,
                 instrum_after=Sequence(elements=instrum.elements[2:]),
                 outfile='SPOs_individual.fits',
                 xyz2zxy=True)
