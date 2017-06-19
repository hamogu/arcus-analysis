from astropy.table import Table
import astropy.table

from marxs.source import PointSource, JitterPointing
import arcus.arcus as arcus

EQPegAspec = Table.read('../inputdata/EQPegA_flux.tbl', format='ascii',
                        names=['energy', 'flux'])
# restrict table to ARCUS energy range
EQPegAspec = EQPegAspec[(EQPegAspec['energy'] > 0.25) & (EQPegAspec['energy'] < 1.5)]

coord = astropy.coordinates.SkyCoord.from_name("EQ Peg A")

# define position and spectrum of source
mysource = PointSource(coords=coord, energy=EQPegAspec, geomarea=arcus.aper4.area)
mypointing = JitterPointing(coords=coord, jitter=arcus.jitter_sigma)

photons = mysource.generate_photons(1e4)
photons = mypointing(photons)
photons = arcus.arcus4(photons)

photons.write('/melkor/d1/guenther/Dropbox/ARCUS/rays/EQPegA10ks.fits', overwrite=True)
