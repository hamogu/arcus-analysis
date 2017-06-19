import numpy as np

from marxs.source import PointSource, FixedPointing
import astropy.units as u
from astropy.coordinates import SkyCoord
import arcus
import arcus.arcus

from mayavi import mlab
from marxs import visualization
from marxs.visualization.mayavi import plot_object, plot_rays
import json
%matplotlib

n_photons = 1e4
wave = np.arange(8., 50., 0.5) * u.Angstrom
energies = wave.to(u.keV, equivalencies=u.spectral()).value
outfile = '../results/aeff.fits'

e = 0.5

mysource = PointSource(coords=SkyCoord(30. * u.deg, 30. * u.deg),
                       energy=e, flux=1.)
photons = mysource.generate_photons(n_photons)

mypointing = FixedPointing(coords=SkyCoord(30 * u.deg, 30. * u.deg))
photons = mypointing(photons)

photons = arcus.arcus.arcus4(photons)


ind = (photons['probability'] > 0)
posdat = visualization.utils.format_saved_positions(arcus.arcus.keeppos4)[ind, :, :]
fig = mlab.figure()
obj = plot_object(arcus.arcus.arcus4, viewer=fig)
rays = plot_rays(posdat, scalar=photons['probability'][ind])
