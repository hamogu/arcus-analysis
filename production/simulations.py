from os.path import join as pjoin
import time
import numpy as np
import astropy.units as u
from marxs.analysis import ProjectOntoPlane
from arcus.arcus import (Arcus, PerfectArcus, defaultconf,
                         DetCamera,
                         CircularDetector)
from arcus.defaults import DefaultSource, DefaultPointing
from arcus.analysis.grid import aeffRfromraygrid, csv_per_order
from arcus import xyz2zxy
from utils import get_path

n_photons = 400000
wave = np.arange(1.5, 50., 0.15) * u.Angstrom

energies = wave.to(u.keV, equivalencies=u.spectral()).value

mypointing = DefaultPointing()


# Add a circular detector to standard Arcus definitions
class ArcusDet(Arcus):
    def add_detectors(self, conf):
        twostrips = DetCamera(conf)
        proj = ProjectOntoPlane(orientation=xyz2zxy[:3, :3])
        detcirc = CircularDetector(conf['rowland_detector'], 'circ')
        return [twostrips, proj, detcirc]


class PerfectArcusDet(PerfectArcus):
    def add_detectors(self, conf):
        twostrips = DetCamera(conf)
        proj = ProjectOntoPlane(orientation=xyz2zxy[:3, :3])
        detcirc = CircularDetector(conf['rowland_detector'], 'circ')
        return [twostrips, proj, detcirc]


for instrum, path in zip([ArcusDet(), PerfectArcusDet()],
                         ['raygrid-small', 'raygrid-perfect-small']):
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
