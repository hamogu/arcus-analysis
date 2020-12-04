from os.path import join
import time
import numpy as np
import astropy.units as u
from arcus.arcus import (Arcus, PerfectArcus, defaultconf)
from arcus.defaults import DefaultSource, DefaultPointing
from arcus.analysis.grid import aeffRfromraygrid, csv_per_order
from utils import get_path

n_photons = 400000
wave = np.arange(1.5, 50., 0.15) * u.Angstrom

energies = wave.to(u.keV, equivalencies=u.spectral()).value

mypointing = DefaultPointing()

for instrum, path in zip([#Arcus(),
                          PerfectArcus()
                      ],
                         [#'raygrid',
                          'raygrid-perfect'
                      ]):
    outpath = get_path(path)

    for i in range(len(wave)):
        print('{0}/{1} = {2}'.format(i + 1, len(energies), time.ctime()))
        mysource = DefaultSource(energy=energies[i])

        photons = mysource.generate_photons(n_photons)
        photons = mypointing(photons)
        photons = instrum(photons)

        photons.write(join(outpath,
                            'wave{0:05.2f}.fits'.format(wave.value[i])),
                      overwrite=True)
    aeffRfromraygrid(outpath, instrum.elements[0], defaultconf,
                     join(get_path('arcus'), path + 'RAeff.fits'))

csv_per_order(join(get_path('arcus'), 'raygridRAeff.fits'), 'Aeff4',
              join(get_path('figures'), 'Aeff.csv'))
csv_per_order(join(get_path('arcus'), 'raygridRAeff.fits'), 'R4',
              join(get_path('figures'), 'R.csv'))
csv_per_order(join(get_path('arcus'), 'raygrid-perfectRAeff.fits'), 'Aeff4',
              join(get_path('figures'), 'Aeff-perfect.csv'))
csv_per_order(join(get_path('arcus'), 'raygrid-perfectRAeff.fits'), 'R4',
              join(get_path('figures'), 'R-perfect.csv'))
