# coding: utf-8
import os
from copy import deepcopy
import transforms3d
import numpy as np
from astropy.coordinates import SkyCoord
import marxs
from marxs.visualization.mayavi import plot_rays, plot_object
from marxs import optics, source, simulator
from marxs.design import rowland

from mayavi import mlab

x3dpath = '../../temp/x3d/'

get_ipython().magic(u'matplotlib')

# Make a model that is not to scale for display purposes.
# Pick a small Rowland radius and make diffraction angle much larger.


demo_rowland = rowland.RowlandTorus(500, 500)

demo_aperture = optics.CircleAperture(position=[1200, 0, 0],
                                      zoom=[1, 100, 100])
demo_mirror = optics.FlatStack(position=[1100, 0, 0], zoom=[20, 100, 100],
                               elements=[optics.PerfectLens,
                                         optics.RadialMirrorScatter],
                               keywords=[{'focallength': 1100},
                                         {'inplanescatter': 1e-3,
                                          'perpplanescatter': 1e-4}])

demo_gas_kwargs = {"rowland": demo_rowland, 'd_element': 25,
                   "x_range": [800, 1000],
                   "radius": [20, 100],
                   "elem_class": marxs.optics.FlatGrating,
                   "elem_args": {'d': 1e-5, 'zoom': [1, 10, 10],
                                 'order_selector': marxs.optics.OrderSelector([-1, 0, 1])}
                  }

demo_gas = rowland.GratingArrayStructure(**demo_gas_kwargs)
demo_gas_sub1 = rowland.GratingArrayStructure(phi=[-0.5 + 0.5 * np.pi, .5 + 0.5 * np.pi], **demo_gas_kwargs)
demo_gas_sub2 = rowland.GratingArrayStructure(phi=[-0.5 + 1.5 * np.pi, .5 + 1.5 * np.pi], **demo_gas_kwargs)

demo_det_kwargs = {'d_element': 10.5, 'elem_class': optics.FlatDetector,
                  'elem_args': {'zoom': [1, 5, 5], 'pixsize': 0.01}}
demo_det = rowland.RowlandCircleArray(demo_rowland, theta=[np.pi - 0.8, np.pi + 0.8], **demo_det_kwargs)

detfp = optics.FlatDetector(pixsize=1., zoom=[1, 500, 500])
detfp.display = deepcopy(detfp.display)
detfp.display['opacity'] = 0.2

demo_det.elements[0].display['color'] = 'orange'

demo_onaxis_full = simulator.Sequence(elements=[demo_aperture, demo_mirror,
                                                demo_gas, demo_det, detfp])

demo_onaxis_sub = simulator.Sequence(elements=[demo_aperture, demo_mirror,
                                               demo_gas_sub1, demo_gas_sub2,
                                               demo_det, detfp])


def make_x3dplot(instrument):
    keeppos = simulator.KeepCol('pos')
    instrument.postprocess_steps = [keeppos]
    target = SkyCoord(30., 30., unit='deg')
    star = source.PointSource(coords=target, energy=.5, flux=1.)
    pointing = source.FixedPointing(coords=target)
    photons = star.generate_photons(5000)
    photons = pointing(photons)
    photons = instrument(photons)
    ind = (photons['probability'] >= 0) & (photons['facet'] >= 0)
    posdat = keeppos.format_positions()[ind, :, :]
    pp = photons[ind]
    fig = mlab.figure()
    plot_object(instrument, viewer=fig)
    plot_rays(posdat, scalar=pp['order'])
    return fig

fig = make_x3dplot(demo_onaxis_full)
demo_rowland.display['coo2'] = np.linspace(-.2, .2, 60)
obj = plot_object(demo_rowland, viewer=fig)
mlab.savefig(os.path.join(x3dpath, 'toy_chandralike.x3d'))

fig = make_x3dplot(demo_onaxis_sub)
obj = plot_object(demo_rowland, viewer=fig)
mlab.savefig(os.path.join(x3dpath, 'toy_subaper.x3d'))

# ## CAT gratings
#
# - Are illuminated at an angle, thus they work better with a tilted Rowland torus.
# - Diffract almost exculsively to one side.
# - Diffract efficiently into higher orders.
# - Need CCDs to image those high orders, but also want to see zeroth order for wavelength calibration.


alpha = 0.4
beta = 0.8
R, r, pos4d = rowland.design_tilted_torus(1e3, alpha, beta)


demo_cat_rowland = rowland.RowlandTorus(R, r, pos4d=pos4d)
demo_cat_aper = optics.RectangleAperture(position=[1200, 0, 0], zoom=[1, 30, 100])

blazeang = 10
blazemat = transforms3d.axangles.axangle2mat(np.array([0, 0, 1]), np.deg2rad(-blazeang))

demo_cat_gas_kwargs = deepcopy(demo_gas_kwargs)
demo_cat_gas_kwargs['rowland'] = demo_cat_rowland
demo_cat_gas_kwargs["elem_class"] = marxs.optics.CATGrating
demo_cat_gas_kwargs['elem_args']['order_selector'] = optics.OrderSelector([0, -1, -2])
demo_cat_gas_kwargs['elem_args']['orientation'] = blazemat

z_offset_spectra = 0 # later if I want to offset to one side of the detector
demo_cat_gas_kwargs['normal_spec'] = np.array([0, 0., -z_offset_spectra, 1.])

demo_gas_1 = rowland.GratingArrayStructure(phi=[-0.5 + 0.5 * np.pi, .5 + 0.5 * np.pi], **demo_cat_gas_kwargs)
demo_gas_2 = rowland.GratingArrayStructure(phi=[-0.5 + 1.5 * np.pi, .5 + 1.5 * np.pi], **demo_cat_gas_kwargs)

demo_cat_gas = simulator.Sequence(elements=[demo_gas_1, demo_gas_2])
demo_cat_det = rowland.RowlandCircleArray(demo_cat_rowland, theta=[np.pi - 0.8, np.pi + 1.2], **demo_det_kwargs)

demo_cat = simulator.Sequence(elements=[demo_cat_aper, demo_mirror,
                                        demo_cat_gas, demo_cat_det, detfp])

fig = make_x3dplot(demo_cat)
obj = plot_object(demo_cat_rowland, viewer=fig)
mlab.savefig(os.path.join(x3dpath, 'toy_cat.x3d'))

# ## Double tilted Rowland design

Rm, rm, pos4dm = rowland.design_tilted_torus(1e3, - alpha, -beta)
rowlandm = rowland.RowlandTorus(Rm, rm, pos4d=pos4dm)
d = r * np.sin(alpha)
# Relative to z=0 in the center of the CCD strip
shift_optical_axis_2 = np.eye(4)
shift_optical_axis_2[1, 3] = 2. * d

rowlandm.pos4d = np.dot(shift_optical_axis_2, rowlandm.pos4d)

demo_cat_aper2 = optics.RectangleAperture(pos4d=np.dot(shift_optical_axis_2, demo_cat_aper.pos4d))
demo_cat_aperm = optics.MultiAperture(elements=[demo_cat_aper, demo_cat_aper2])
demo_mirrorm = optics.FlatStack(pos4d=np.dot(shift_optical_axis_2, demo_mirror.pos4d),
                                elements=[optics.PerfectLens, optics.RadialMirrorScatter],
                                keywords = [{'focallength': 1100},
                                         {'inplanescatter': 1e-3, 'perpplanescatter': 1e-4}])

blazematm = transforms3d.axangles.axangle2mat(np.array([0, 0, 1]), np.deg2rad(blazeang))

demo_cat_gas_kwargsm = deepcopy(demo_cat_gas_kwargs)
demo_cat_gas_kwargsm['rowland'] = rowlandm
demo_cat_gas_kwargsm['normal_spec'] = np.array([0, 2 * d, 0, 1.])
demo_cat_gas_kwargsm['x_range'] = [-800, -1000]
demo_cat_gas_kwargsm['elem_args']['orientation'] = blazematm

demo_gas_1m = rowland.GratingArrayStructure(phi=[-0.5 + 0.5 * np.pi, .5 + 0.5 * np.pi], **demo_cat_gas_kwargsm)
demo_gas_2m = rowland.GratingArrayStructure(phi=[-0.5 + 1.5 * np.pi, .5 + 1.5 * np.pi], **demo_cat_gas_kwargsm)

demo_cat_gasm = simulator.Sequence(elements=[demo_gas_1, demo_gas_2,
                                             demo_gas_1m, demo_gas_2m])

demo_catm = simulator.Sequence(elements=[demo_cat_aperm, demo_mirror,
                                         demo_mirrorm,
                                         demo_cat_gasm, demo_cat_det, detfp])

fig = make_x3dplot(demo_catm)
obj = plot_object(demo_cat_rowland, viewer=fig)
obj = plot_object(rowlandm, viewer=fig)
mlab.savefig(os.path.join(x3dpath, 'toy_cat2.x3d'))

fig = make_x3dplot(demo_catm)
demo_rowland.display['coo2'] = np.linspace(-np.pi, np.pi, 60)
obj = plot_object(demo_cat_rowland, viewer=fig)
obj = plot_object(rowlandm, viewer=fig)
mlab.savefig(os.path.join(x3dpath, 'toy_cat2_tori.x3d'))
