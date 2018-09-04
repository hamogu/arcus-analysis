import numpy as np

import marxs
import arcus
import arcus.arcus

from mayavi import mlab
from marxs.math.utils import h2e
from marxs import visualization
from marxs.visualization.mayavi import plot_object, plot_rays
%matplotlib

from arcus.arcus import Arcus
from arcus.defaults import DefaultSource, DefaultPointing


n_photons = 1e4
e = 0.5

mypointing = DefaultPointing()

mysource = DefaultSource(energy=e)

photons = mysource.generate_photons(n_photons)
photons = mypointing(photons)

instrum = Arcus()
photons = instrum(photons)


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
