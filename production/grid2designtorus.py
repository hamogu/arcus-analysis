from __future__ import print_function
import time
import os
import numpy as np
import astropy.units as u
from arcus import Arcus
import arcus.arcus
from arcus.defaults import DefaultSource, DefaultPointing
from arcus.generate_rowland import make_rowland_from_d_BF_R_f

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



arr_R = np.arange(5900., 6001., 100.)
arr_d_BF = np.arange(500., 701., 50.)
arr_blaze = np.arange(1.2, 2.21, 0.2)

for R in arr_R:
    for d_BF in arr_d_BF:
        conf = make_rowland_from_d_BF_R_f(d_BF, R)
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
