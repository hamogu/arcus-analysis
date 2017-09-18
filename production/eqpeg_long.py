import os
import numpy as np
from astropy.table import Table
import astropy.table

from arcus import Arcus
from arcus.defaults import DefaultSource, DefaultPointing

from utils import get_path

EQPegAspec = Table.read('../../inputdata/EQPegA_flux.tbl', format='ascii',
                        names=['energy', 'flux'])
# restrict table to ARCUS energy range
EQPegAspec = EQPegAspec[(EQPegAspec['energy'] > 0.25) &
                        (EQPegAspec['energy'] < 1.5)]

coord = astropy.coordinates.SkyCoord.from_name("EQ Peg A")

# define position and spectrum of source
arc = Arcus()
mysource = DefaultSource(coords=coord, energy=EQPegAspec,
                         geomarea=arc.elements[0].area,
                         flux=(EQPegAspec['flux'][1:] * np.diff(EQPegAspec['energy'])).sum())
mypointing = DefaultPointing(coords=coord)

photons = mysource.generate_photons(1e5)
photons = mypointing(photons)
photons = arc(photons)

photons.write(os.path.join(get_path('rays'), 'EQPegA100ks.fits'),
              overwrite=True)
