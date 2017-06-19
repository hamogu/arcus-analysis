import numpy as np

from marxs.source import PointSource, FixedPointing
import marxs
import astropy.units as u
from astropy.coordinates import SkyCoord
import arcus
import arcus.arcus

from mayavi import mlab
from marxs.math.utils import h2e
from marxs import visualization
from marxs.visualization.mayavi import plot_object, plot_rays


n_photons = 1e4
wave = np.arange(8., 50., 0.5) * u.Angstrom
#energies = np.arange(.2, 1.9, .01)
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
rays = plot_rays(posdat, scalar=photons['energy'][ind])


d = np.dstack(keeppos.data)
d = np.swapaxes(d, 1, 2)
d = h2e(d)

marxs.visualization.mayavi.plot_rays(d, scalar=photons['order'], viewer=fig)
arcus.arcus.plot(format="mayavi", viewer=fig)

theta, phi = np.mgrid[-0.2 + np.pi:0.2 + np.pi:60j, -1:1:60j]
arcus.rowland.plot(theta=theta, phi=phi, viewer=fig, format='mayavi')


photonsm = mysource.generate_photons(n_photons)
photonsm = mypointing(photonsm)

keepposm = marxs.simulator.KeepCol('pos')
arcus.arcus_extra_det_m.postprocess_steps = [keepposm]

photonsm = arcus.arcus_extra_det_m(photonsm)
d = np.dstack(keepposm.data)
d = np.swapaxes(d, 1, 2)
d = h2e(d)

marxs.visualization.mayavi.plot_rays(d, scalar=photonsm['order'], viewer=fig)
arcus.arcusm.plot(format="mayavi", viewer=fig)


from astropy.stats import sigma_clipped_stats
ind = np.isfinite(photons['order']) & (photons['CCD_ID'] >=0)
