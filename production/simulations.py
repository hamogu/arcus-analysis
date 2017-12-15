import os
import time
import numpy as np
import astropy.units as u
from arcus import Arcus
from arcus.defaults import DefaultSource, DefaultPointing

from utils import get_path

n_photons = 4e5
wave = np.arange(8., 50., 0.15) * u.Angstrom
energies = wave.to(u.keV, equivalencies=u.spectral()).value
outpath = get_path('raygrid')

mypointing = DefaultPointing()

instrum = Arcus()

for i, e in enumerate(energies):
    print '{0}/{1} = {2}'.format(i + 1, len(energies), time.ctime())
    mysource = DefaultSource(energy=e)

    photons = mysource.generate_photons(n_photons)
    photons = mypointing(photons)
    photons = instrum(photons)

    photons.write(os.path.join(outpath,
                               'wave{0:05.2f}.fits'.format(wave.value[i])),
                  overwrite=True)
