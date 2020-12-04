''' '''
from __future__ import print_function
import time
import os
import numpy as np
import astropy.units as u
import arcus.arcus
from arcus.defaults import DefaultSource, DefaultPointing
from arcus.generate_rowland import make_rowland_from_d_BF_R_f
from arcus.ralfgrating import (CATL1L2Stack, RectangularGrid,
                               InterpolateRalfTable,
                               RalfQualityFactor,
                               catsupport,
                               catsupportbars)
from transforms3d.axangles import axangle2aff, axangle2mat
from marxs.simulator import Sequence

from utils import get_path
outpath = get_path('grid2designtorus')

n_photons = 10000
wave = np.arange(8., 50.1, 1.) * u.Angstrom
energies = wave.to(u.keV, equivalencies=u.spectral()).value

mypointing = DefaultPointing()
mysource = DefaultSource()

# Make input photon list with grid of discrete energies
photons_in = mysource.generate_photons(n_photons * len(wave))
photons_in = mypointing(photons_in)

for i, e in enumerate(energies):
    photons_in['energy'][i * n_photons: (i + 1) * n_photons] = e


arr_R = np.arange(5900., 6001., 100.)
arr_d_BF = np.arange(500., 701., 50.)
arr_blaze = np.arange(1.2, 2.21, 0.2)


class CATGratings(Sequence):
    order_selector_class = InterpolateRalfTable
    gratquality_class = RalfQualityFactor
    grid_width_x = 200
    grid_width_y = 300

    def __init__(self, conf, channels=['1', '2', '1m', '2m'], **kwargs):

        elements = []

        self.order_selector = self.order_selector_class()
        self.gratquality = self.gratquality_class()
        blazemat = axangle2mat(np.array([0, 0, 1]),
                               np.deg2rad(-conf['blazeang']))
        blazematm = axangle2mat(np.array([0, 0, 1]),
                                np.deg2rad(conf['blazeang']))

        gratinggrid = {'d_element': 32., 'z_range': [1e4, 1.4e4],
                       'elem_class': CATL1L2Stack,
                       'elem_args': {'zoom': [1., 13.5, 13.]},
                       'parallel_spec': np.array([1., 0., 0., 0.])
                   }
        for chan in channels:
            gratinggrid['rowland'] = conf['rowland_' + chan]
            b = blazematm if 'm' in chan else blazemat
            gratinggrid['elem_args']['orientation'] = b
            gratinggrid['normal_spec'] = conf['pos_opt_ax'][chan].copy()
            xm, ym = conf['pos_opt_ax'][chan][:2].copy()
            sig = 1 if '1' in chan else -1
            x_range = [-self.grid_width_x + xm,
                       +self.grid_width_x + xm]
            y_range = [sig * (600 - ym - self.grid_width_y),
                       sig * (600 - ym + self.grid_width_y)]
            y_range.sort()
            elements.append(RectangularGrid(x_range=x_range, y_range=y_range,
                                            id_num_offset=arcus.arcus.id_num_offset[chan],
                                            **gratinggrid))
        elements.extend([catsupport, catsupportbars, self.gratquality])
        super(CATGratings, self).__init__(elements=elements, **kwargs)


class Arcus(arcus.arcus.PerfectArcus):
    gratings_class = CATGratings


def tiltCATGratings(instrum, defaultblaze, blaze):
    '''Tilt Arcus gratings.

    In the early days of this script, my code would place the gratings
    on the Rowlandtorus depending on the input values in the conf
    dictionary.
    However, with the more mature design, the position and location of the
    gratings is read from CATfromMechanical. This function will twist those
    gratings to keep he position consistent, but experiment with different
    blaze angles.
    '''
    rot = np.deg2rad(blaze - defaultblaze)
    grats = instrum.elements_of_class(CATL1L2Stack)
    for g in grats:
        rotmat = axangle2aff(g.geometry['e_z'][:3], rot)
        newpos = rotmat @ g.geometry.pos4d
        g.geometry.pos4d = newpos
        for ge in g.elements:
            g.geometry.pos4d = newpos


for R in arr_R:
    for d_BF in arr_d_BF:
        conf = make_rowland_from_d_BF_R_f(d_BF, R, 11880.)
        for k in ['phi_det_start', 'n_CCDs']:
            conf[k] = arcus.arcus.defaultconf[k]
        for blaze in arr_blaze:
            conf['blazeang'] = blaze
            filename = '{:04.0f}_{:04.1f}_{:05.0f}.fits'.format(d_BF, blaze, R)
            filepath = os.path.join(outpath, filename)
            if os.path.isfile(filepath):
                continue
            else:
                photons_in.write(filepath)
                print('d_BF: {:4.0f} mm - blaze: {:4.1f} deg - R: {:5.0f} mm - {}'.format(d_BF, blaze, R, time.ctime()))
                instrum = Arcus(channels=['1'], conf=conf)
                #tiltCATGratings(instrum, arcus.arcus.defaultconf['blazeang'], blaze)
                photons_out = instrum(photons_in.copy())
                photons_out.meta['TORUS_R'] = R
                photons_out.meta['CIRCLE_R'] = conf['rowland_1'].r
                photons_out.meta['MAX_F'] = conf['f']
                photons_out.meta['BLAZE'] = blaze
                photons_out.meta['D_CHAN'] = d_BF
                photons_out.meta['N_PHOT'] = n_photons
                photons_out.meta['A_GEOM'] = instrum.elements[0].area.to(u.cm**2).value
                # Keep only those columns absolutely needed for analysis
                # to reduce file size
                photons_out.keep_columns(['energy', 'probability', 'order',
                                          'circ_phi', 'circ_y'])
                photons_out.write(filepath, overwrite=True)
