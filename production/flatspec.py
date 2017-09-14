from os.path import join as pjoin
import datetime
import numpy as np
import astropy.units as u

from arcus import Arcus
from arcus.defaults import DefaultSource, DefaultPointing

from utils import get_path

n_photons = 1e5
outdir = get_path('arcus')

# define position and spectrum of source
mysource = DefaultSource(energy={'energy': np.array([0.25, 1.7]),
                                 'flux': np.ones(2)})
jitterpointing = DefaultPointing()
fixedpointing = DefaultPointing(jitter=0. * u.rad)
arc = Arcus()

for i in range(10):
    print 'jitter: {:03d}'.format(i), ' : ', datetime.datetime.now()
    # Ignore geometric area and set number of photons by hand.
    photons = mysource.generate_photons(n_photons)
    photons = jitterpointing(photons)
    photons = arc(photons)

    photons.write(pjoin(outdir, 'flatspecjitter{:03d}.fits'.format(i)),
                  overwrite=True)

for i in range(10):
    print 'fixed: {:03d}'.format(i), ' : ', datetime.datetime.now()
    # Ignore geometric area and set number of photons by hand.
    photons = mysource.generate_photons(n_photons)
    photons = fixedpointing(photons)
    photons = arc(photons)

    photons.write(pjoin(outdir, 'flatspecfixed{:03d}.fits'.format(i)),
                  overwrite=True)
