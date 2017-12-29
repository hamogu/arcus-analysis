from os.path import join as pjoin
import time
import numpy as np
import astropy.units as u
from arcus.arcus import Arcus, PerfectArcus, defaultconf
from arcus.defaults import DefaultSource, DefaultPointing
from arcus.analysis.grid import aeffRfromraygrid, csv_per_order
from utils import get_path

n_photons = 4e5
wave = np.arange(8., 50., 0.15) * u.Angstrom
energies = wave.to(u.keV, equivalencies=u.spectral()).value
outpath = get_path('raygrid')

mypointing = DefaultPointing()

for instrum, path in zip([Arcus(), PerfectArcus()],
                         ['raygrid', 'raygrid-perfect']):
    outpath = get_path(path)

    for i, e in enumerate(energies):
        print '{0}/{1} = {2}'.format(i + 1, len(energies), time.ctime())
        mysource = DefaultSource(energy=e)

        photons = mysource.generate_photons(n_photons)
        photons = mypointing(photons)
        photons = instrum(photons)

        photons.write(pjoin(outpath,
                            'wave{0:05.2f}.fits'.format(wave.value[i])),
                      overwrite=True)
    aeffRfromraygrid(outpath, instrum.elements[0], defaultconf,
                     pjoin(get_path('arcus'), path + 'RAeff.fits'))

csv_per_order(pjoin(get_path('arcus'), 'raygridRAeff.fits'), 'Aeff4',
              pjoin(get_path('figures'), 'Aeff.csv'))
csv_per_order(pjoin(get_path('arcus'), 'raygridRAeff.fits'), 'R4',
              pjoin(get_path('figures'), 'R.csv'))
csv_per_order(pjoin(get_path('arcus'), 'raygrid-perfectRAeff.fits'), 'Aeff4',
              pjoin(get_path('figures'), 'Aeff-perfect.csv'))
csv_per_order(pjoin(get_path('arcus'), 'raygrid-perfectRAeff.fits'), 'R4',
              pjoin(get_path('figures'), 'R-perfect.csv'))
