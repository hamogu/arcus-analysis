from os.path import join
import time

import numpy as np
import astropy.units as u
import transforms3d
from marxs.optics import CircularDetector, RectangleAperture
from marxs.optics.scatter import RadialMirrorScatter
from marxs.optics.mirror import PerfectLens
from marxs.optics.grating import FlatGrating
from marxs.missions.mitsnl.catgrating import L1

from arcus.instrument.arcus import Arcus, DetCamera
from arcus.instrument.spo import ScatterPerChannel
from arcus.defaults import DefaultSource, DefaultPointing
from arcus.instrument.constants import xyz2zxy

from utils import get_path

n_photons = 5e5
wave = np.arange(1.5, 60., 1.) * u.Angstrom
wave = np.array([1.25, 70, 85, 100, 120]) * u.Angstrom
#wave = np.arange(50.5, 61.9, 1.) * u.Angstrom
#wave = np.array([42]) * u.Angstrom

energies = wave.to(u.keV, equivalencies=u.spectral())

mypointing = DefaultPointing()
channel = '1'

'''400000 photons with all columns saved is about 162 MB. So, for
larger runs, defintely need to limit the columns saved to those that I
acually need.
With the settings below and keep_column in the end, I can bring that
down to ~6 MB.
'''
ScatterPerChannel.loc_coos_name = None
RadialMirrorScatter.inplanescattercol = None
RadialMirrorScatter.perpplanescattercol = None
PerfectLens.loc_coos_name = None
FlatGrating.loc_coos_name = None
FlatGrating.blaze_name = None
L1.blaze_name = None


class MyArcus(Arcus):
    filter_and_qe_class = None

    def add_detectors(self, conf):
        # rotatate such that phi=0 is at the bottom
        rot = transforms3d.axangles.axangle2mat(np.array([0, 1, 0]), np.pi)
        rot2 = transforms3d.axangles.axangle2mat(np.array([1, 0, 0]), np.pi)
        circdet = CircularDetector(orientation=xyz2zxy[:3, :3] @ rot @ rot2,
                                   zoom=[conf['rowland_detector'].r,
                                         conf['rowland_detector'].r,
                                         500],
                                   position=[0, 0, np.sqrt(conf['rowland_detector'].r**2 - conf['d']**2)],
                                   )
        circdet.display['opacity'] = 0.1
        circdet.detpix_name = [None, None]
        circdet.loc_coos_name = ['circ_phi', 'circ_y']
        return [circdet]
        #reset = marxs.simulator.simulator.Propagator(distance=-100.)
        #twostrips = DetCamera(conf)
        #return [circdet, reset, twostrips]


instrum = MyArcus(channels=[channel])
outpath = get_path('arfrmf')

for i, (en, wav) in enumerate(zip(energies, wave)):
    print('{0}/{1} = {2}'.format(i + 1, len(energies), time.ctime()))
    mysource = DefaultSource(energy=en)

    photons = mysource.generate_photons(n_photons * u.s)
    photons = mypointing(photons)
    photons.remove_columns(['time', 'polangle', 'ra', 'dec'])
    photons = instrum(photons)
    photons.meta['ENERGY'] = (en.value, en.unit)
    photons.meta['WAVELEN'] = (wav.value, wav.unit)
    photons.meta['CHANNEL'] = channel
    for aper in instrum.elements_of_class(RectangleAperture):
        area = aper.area.to(u.cm**2)
        photons.meta[f'AGEOM_{aper.id_num}'] = (area.value, area.unit)
    aper = instrum.elements[0]
    area = aper.area.to(u.cm**2)
    photons.meta['A_GEOM'] = (area.value, area.unit)
    # Number of photons is saved in meta['EXPOSURE'] with this setup
    # so no need to waste disk space for photons with probability=0
    photons.keep_columns(['probability', 'order', 'order_L1',
                          'circ_phi', 'circ_y'])
    photons = photons[photons['probability'] > 0]
    photons.write(join(outpath,
                           'wave{0:05.2f}.fits'.format(wav.value)),
                      overwrite=True)
