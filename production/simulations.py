import os
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

from marxs.source import PointSource, JitterPointing
from marxs.simulator import Sequence

import arcus.arcus as arcus

n_photons = 1e4
wave = np.arange(8., 50., 0.5) * u.Angstrom
energies = wave.to(u.keV, equivalencies=u.spectral()).value
outpath = '../../../../Dropbox/ARCUS/rays/semi-compact'

mypointing = JitterPointing(coords=SkyCoord(30, 30., unit='deg'),
                            jitter=arcus.jitter_sigma)

instrum = Sequence(elements=[arcus.aper4,
                             arcus.mirror4,
                             arcus.gas4,
                             arcus.filtersandqe,
                             arcus.det, arcus.projectfp,
                             arcus.tagversion])

for i, e in enumerate(energies):
    print '{0}/{1}'.format(i + 1, len(energies))
    mysource = PointSource(coords=SkyCoord(30., 30., unit='deg'),
                           energy=e, flux=1.)

    photons = mysource.generate_photons(n_photons)
    photons = mypointing(photons)
    photons = instrum(photons)

    photons.write(os.path.join(outpath,
                               'wave{0:05.2f}.fits'.format(wave.value[i])),
                  overwrite=True)
