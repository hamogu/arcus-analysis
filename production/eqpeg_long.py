import numpy as np
from astropy.table import Table
import astropy.table

from marxs.source import PointSource, JitterPointing
import arcus.arcus as arcus

EQPegAspec = Table.read('../../inputdata/EQPegA_flux.tbl', format='ascii',
                        names=['energy', 'flux'])
# restrict table to ARCUS energy range
EQPegAspec = EQPegAspec[(EQPegAspec['energy'] > 0.25) & (EQPegAspec['energy'] < 1.5)]

coord = astropy.coordinates.SkyCoord.from_name("EQ Peg A")

# define position and spectrum of source
mysource = PointSource(coords=coord, energy=EQPegAspec,
                       geomarea=arcus.aper4.area,
                       flux = (EQPegAspec['flux'][1:] * np.diff(EQPegAspec['energy'])).sum())
mypointing = JitterPointing(coords=coord, jitter=arcus.jitter_sigma)

photons = mysource.generate_photons(1e5)
photons = mypointing(photons)
photons = arcus.arcus_extra_det4(photons)

photons.write('/melkor/d1/guenther/Dropbox/ARCUS/rays/EQPegA100ks.fits', overwrite=True)
