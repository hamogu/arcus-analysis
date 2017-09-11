import os
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

from marxs.source import PointSource, JitterPointing

from arcus import Arcus, jitter_sigma
from arcus.defaults import DefaultSource, DefaultPointing

from utils import get_path

n_photons = 1e4
wave = np.arange(8., 50., 0.5) * u.Angstrom
energies = wave.to(u.keV, equivalencies=u.spectral()).value
outpath = get_path('raygrid')

mypointing = DefaultPointing()

instrum = Arcus()

for i, e in enumerate(energies):
    print '{0}/{1}'.format(i + 1, len(energies))
    mysource = DefaultSource(energy=e)

    photons = mysource.generate_photons(n_photons)
    photons = mypointing(photons)
    photons = instrum(photons)

    photons.write(os.path.join(outpath,
                               'wave{0:05.2f}.fits'.format(wave.value[i])),
                  overwrite=True)
