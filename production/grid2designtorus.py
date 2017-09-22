from __future__ import print_function
import time
import os
import numpy as np
import astropy.units as u
from marxs.design import RowlandTorus
from arcus import Arcus
import arcus.arcus
from arcus.defaults import DefaultSource, DefaultPointing

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


class ArcusCirc(Arcus):
    def add_detectors(self, conf):
        return [arcus.arcus.CircularDetector(conf['rowland_central'], 'circ'),
                arcus.arcus.FocalPlaneDet()]


def make_rowland(d_BF, R, f=11880.):
    r = 0.5 * np.sqrt(f**2 + d_BF**2)
    alpha = np.arctan2(d_BF, f)
    pos = [(r + R) * np.sin(alpha), 0, f - (r + R) * np.cos(alpha)]
    orient = [[-np.sin(alpha), np.cos(alpha), 0],
              [0., 0., 1],
              [np.cos(alpha), np.sin(alpha), 0]]
    geometry = {'d_BF': d_BF,
                'd': 0.5 * d_BF,
                'f': f,
                'offset_spectra': 5.,
                'rowland_central': RowlandTorus(R=R, r=r, position=pos, orientation=orient)}

    # Now offset that Rowland torus in a z axis by a few mm.
    # Shift is measured from a focal point that hits the center of the CCD strip.
    geometry['shift_optical_axis_1'] = np.eye(4)
    geometry['shift_optical_axis_1'][1, 3] = - geometry['offset_spectra']

    geometry['rowland'] = RowlandTorus(R, r, pos4d=geometry['rowland_central'].pos4d)
    geometry['rowland'].pos4d = np.dot(geometry['shift_optical_axis_1'],
                                       geometry['rowland'].pos4d)

    # Not really needed for a 1 channel run, but some functions in arcus.py
    # expect those keywords to be present.
    geometry['shift_optical_axis_2'] = np.eye(4)
    geometry['shift_optical_axis_2'][0, 3] = d_BF
    geometry['shift_optical_axis_2'][1, 3] = + geometry['offset_spectra']

    # Optical axis 2 relative to optical axis 1
    geometry['shift_optical_axis_12'] = np.dot(np.linalg.inv(geometry['shift_optical_axis_1']),
                                               geometry['shift_optical_axis_2'])


    return geometry


arr_R = np.arange(5800., 6001., 100.)
arr_d_BF = np.arange(600., 951., 50.)
arr_blaze = np.arange(1.2, 2.21, 0.1)

for R in arr_R:
    for d_BF in arr_d_BF:
        conf = make_rowland(d_BF, R)
        for blaze in arr_blaze:
            conf['blazeang'] = blaze
            filename = '{:04.0f}_{:04.1f}_{:05.0f}.fits'.format(d_BF, blaze, R)
            filepath = os.path.join(outpath, filename)
            if os.path.isfile(filepath):
                continue
            else:
                photons_in.write(filepath)
                print('d_BF: {:4.0f} mm - blaze: {:4.1f} deg - R: {:5.0f} mm - {}'.format(d_BF, blaze, R, time.ctime()))
                instrum = ArcusCirc(channels=['1'], conf=conf)
                photons_out = instrum(photons_in.copy())
                photons_out.meta['TORUS_R'] = R
                photons_out.meta['CIRCLE_R'] = conf['rowland'].r
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
