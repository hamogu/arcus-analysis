from __future__ import print_function
from os.path import join as pjoin
import time
import numpy as np

from marxs.math.utils import h2e, norm_vector
from astropy.table import Table
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord

from arcus.arcus import ArcusForSIXTE
from arcus.defaults import DefaultPointing, DefaultSource

from utils import get_path

n_photons = 5e5

wave = np.arange(8., 50., 0.5) * u.Angstrom
energies = wave.to(u.keV, equivalencies=u.spectral()).value

outdir = get_path('forsixte')

pointing_offsets = [0 * u.degree, 0.1 * u.degree]

arc = ArcusForSIXTE()

spectab = Table([energies], names=['ENERGY'])
spectab['ENERGY'].unit = u.keV
spectab['N_PHOT'] = n_photons
spectab.meta['HDUNAME'] = 'SPECTRUM'


def write_joerntables(photons, outdir, ie, indx, offx, indy, offy):
    if not np.all(photons['energy'] == photons['energy'][0]):
        raise ValueError("need monoenergetic simulations")
    orders = set(photons['order'][np.isfinite(photons['order'])])
    for o in orders:
        filename = pjoin(outdir, '{0}_{1}_{2}_{3}.fits'.format(ie, indx, indy, int(np.abs(o))))
        ind = (photons['order'] == o)
        tab = photons[ind]
        tab.write(filename, overwrite=True)
        # easier to open again to add keywords then use
        # fits interface directly above
        hdulist = fits.open(filename, mode='update')
        hdulist[0].header['ENERGY'] = (photons[0]['energy'], 'energy in keV')
        hdulist[0].header['ORDER'] = (o, 'diffraction order')
        hdulist[0].header['OFFX'] = (offx.to(u.rad).value,
                                     'offset from optical axis in radian')
        hdulist[0].header['OFFY'] = (offy.to(u.rad).value,
                                     'offset from optical axis in radian')
        hdulist[0].header['NPHOTONS'] = (n_photons, 'Number of photons per simulation')
        hdulist.close()

for ix, offx in enumerate(pointing_offsets):
    for iy, offy in enumerate(pointing_offsets):
        for ie, e in enumerate(energies):
            print('{} {} {} - {}'.format(ix, iy, ie, time.ctime()))
            mysource = DefaultSource(coords=SkyCoord(0. * u.rad, 0. * u.rad), energy=e)
            photons = mysource.generate_photons(n_photons)
            offsetcoord = SkyCoord(offx, offy)
            mypointing = DefaultPointing(coords=offsetcoord, jitter=0.)
            photons = mypointing(photons)
            photons = arc(photons)
            # Reformat and delete columns not required for SIXTE to save space
            photons['POS'] = h2e(photons['pos'])
            # Make sure direction is normalized
            photons['DIR'] = norm_vector(h2e(photons['dir']))
            photons.rename_column('probability', 'weight')
            photons.rename_column('aperture', 'channel')
            photons.keep_columns(['POS', 'DIR', 'time', 'weight', 'ra', 'dec',
                                  'channel', 'order', 'energy'])
            photons.meta['A_GEOM'] = (arc.elements[0].area.to(u.cm**2).value, 'Geometric opening area in cm')
            filename = '{0}_{1}_{2}.fits'.format(ie, ix, iy)
            # Drop photons that are absorbed or miss grating
            photons = photons[np.isfinite(photons['order']) & (photons['weight'] > 0.)]
            photons.write(pjoin(outdir, filename), overwrite=True)
            # Add spectrum hdu
            spechdu = fits.table_to_hdu(spectab[[ie]])  # double [[]] to get table, not row
            with fits.open(pjoin(outdir, filename), 'append') as hdulist:
                hdulist.append(spechdu)
                hdulist.flush()
