from __future__ import print_function

import os
import time
import numpy as np
from astropy.table import Table
import astropy.table
import astropy.units as u
from marxs.source import DiskSource
from marxs.source.source import poisson_process
from arcus.defaults import DefaultPointing, coords
from arcus import Arcus

from utils import get_path

bkgspec = Table.read('../inputdata/Diffuse_mdl_2.csv', comment='#')
bkgspec.rename_column('whi', 'energy')

# restrict table to ARCUS energy range
bkgspec = bkgspec[(bkgspec['energy'] > 0.25) &
                  (bkgspec['energy'] < 1.5)]

# units of flux is  phot/cm^2/s/bin/arcmin^2

requested_time = 5e4
simulated_time = 0
i = 0

instrum = Arcus()
a1 = 1. * u.deg
# circle with r = a1 has about pi * r **2 arcmin**2
# Could do that for real, but here it's good enough.
fluxperarea = bkgspec['flux'].sum() * np.pi * 3600
diffusesrc = DiskSource(a_outer=a1,
                        flux=poisson_process(fluxperarea),
                        energy=bkgspec,
                        geomarea=instrum.elements[0].area,
                        coords=coords)
pointing = DefaultPointing()

# Simulations are big, so it's more efficient to run them in several
# smaller runs
while simulated_time < requested_time:
    photons = diffusesrc.generate_photons(1e3)
    photons = pointing(photons)
    photons = instrum(photons)

    # cut away photons that miss the detector
    pdet = photons[np.isfinite(photons['det_x'])]
    # Fudge the probabilities before drawing to mimick a longer observations
    maxprob = np.max(pdet['probability'])
    pdet['probability'] /= maxprob
    pdet['time'] /= maxprob
    pdet.meta['EXPOSURE'] = (pdet.meta['EXPOSURE'][0] / maxprob,
                             pdet.meta['EXPOSURE'][1])
    pobs = pdet[pdet['probability'] > np.random.uniform(size=len(pdet))]
    # add offset to time
    pobs['time'] += simulated_time
    # Save in case of crash
    pobs.write(os.path.join(get_path('rays'), 'bkg_{}.fits'.format(i)),
               overwrite=True)
    simulated_time += pdet.meta['EXPOSURE'][0]
    print('{}: Finished simulation {} reaching {} of {} s'.format(time.ctime(), i, simulated_time, requested_time))
    i += 1

tablist = [Table.read(os.path.join(get_path('rays'),
                                   'bkg_{}.fits'.format(j))) for j in range(i)]
fulltab = astropy.table.vstack(tablist)
fulltab.meta['EXPOSURE'] = simulated_time
fulltab.write(os.path.join(get_path('rays'), 'longbkg.fits'.format(i)),
              overwrite=True)
