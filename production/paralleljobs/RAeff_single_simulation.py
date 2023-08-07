# Licensed under GPL version 3 - see LICENSE.rst
'''Run a single monoenenergetic simulation of Arcus.

This is indented for running grids of with different energies
to build up effective area curves or ARF/RMF. Such a calculation
is "embarrassingly parallel" and can easily be run as a number of
different Python processes. To facilitate that use, this is written
as a Python script, not an importable module.

Further Optimization
====================
Further optimization of this script is possible to reduce file sizes
and memory needs during a run. However, at this point it seems that
that is not worth the effort, so ideas are just kept below as notes
should this become important in the future.

Do no generate local coordinate names
-------------------------------------
The easiest way to achieve that is probably to change the defaults in
the classes before creating the instrument objects:

   >>> from marxs.optics.scatter import RadialMirrorScatter
   >>> RadialMirrorScatter.inplanescattercol = None
   >>> RadialMirrorScatter.perpplanescattercol = None

And similarly for `` ScatterPerChannel.loc_coos_name = None``,
``PerfectLens.loc_coos_name = None``,
``FlatGrating.loc_coos_name = None``
``FlatGrating.blaze_name = None``,
``L1.blaze_name = None``.

Have fewer detectors
--------------------
`marxs.missions.arcus.arcus.Arcus` has several detectors for analysis,
but we probably don't need all of them here. Could save some time by
not having all of them.
'''

import sys
from os.path import join, exists
import time
import argparse
import numpy as np
import astropy.units as u
from marxs.missions.arcus.arcus import (Arcus, PerfectArcus,
                                        defaultconf as conf)
from marxs.missions.arcus.defaults import DefaultSource, DefaultPointing
from marxs.optics import RectangleAperture
from marxs.missions.arcus.utils import id_num_offset

parser = argparse.ArgumentParser(
                    description = 'Run a single simulation in a wavelength grid that can later we used to generate R/Aeff curves or even ARFs and RMFs.\n' + 
                       'The interface is specifically designed to make is easy to iterate over this by changing just the "i" parameter' +
                       'which helps with scheduling systems list slurm job arrays.')
parser.add_argument('outpath', help='Path to write simulation output') 
parser.add_argument('i', type=int, help='index in wavelength array to select wavelength for this simulation') 
parser.add_argument('-n', '--n_photons', default=100000, type=int,
                    help='number of photons')
parser.add_argument('--wave_lo', default=1.5, type=float,
                    help='lower end of wavelength range in Ang')
parser.add_argument('--wave_hi', default=60., type=float,
                    help='upper end of wavelength range in Ang')
parser.add_argument('--wave_step', default=0.15, type=float,
                    help='step of wavelength grid in Ang')
parser.add_argument('--channels', nargs='*',
                    default=list(id_num_offset.keys()),
                    help='Select channels for the simulation. Default: all 4 channels')
parser.add_argument('--perfect',
                    action='store_true',
                    help='Use perfect Arcus without default misalignments') 
default_keep_cols = ['probability', 'xou', 'order', 'order_L1',
                     'circ_phi', 'circ_y', 'CCD', 'proj_x', 'proj_y']
parser.add_argument('--keepcols', nargs='*',
                    default= default_keep_cols,
                    help=f'Select which columns to save. Default: {" ".join(default_keep_cols)}')
parser.add_argument('--skip_existing',
                    action='store_true',
                    help='Skip existing files? Default is to re-run and overwrite')

args = parser.parse_args()


wave = np.arange(args.wave_lo, args.wave_hi, args.wave_step) * u.Angstrom
energies = wave.to(u.keV, equivalencies=u.spectral())
en = energies[args.i]
wav = wave[args.i]
filename = join(args.outpath, 'wave{0:05.2f}.fits'.format(wave.value[args.i]))
if args.skip_existing and exists(filename):
   print('Skipping existing file {filename}.')
   sys.exit(0)

pointing = DefaultPointing()
instrument = PerfectArcus(channels=args.channels) if args.perfect else Arcus(channels=args.channels)
print(f'{time.ctime()} - begin simulation for {wave[args.i]:05.2f} Ang')
source = DefaultSource(energy=en)
photons = source.generate_photons(args.n_photons * u.s)
photons = pointing(photons)
photons = instrument(photons)

photons.meta['ENERGY'] = (en.value, en.unit)
photons.meta['WAVELEN'] = (wav.value, wav.unit)
photons.meta['CHANNEL'] = args.channels
for aper in instrument.elements_of_class(RectangleAperture):
   area = aper.area.to(u.cm**2)
   photons.meta[f'AGEOM_{aper.id_num}'] = (area.value, area.unit)
aper = instrument.elements[0]
area = aper.area.to(u.cm**2)
photons.meta['A_GEOM'] = (area.value, area.unit)
photons.meta['alignmen'] = 'perfect' if args.perfect else 'applied'
# There is no way we need all those columns. Waste of space.
photons.keep_columns(
 [
 'probability',
 'xou',
 'order',
 'order_L1',
 'circ_phi',
 'circ_y',
 'CCD',
 'proj_x',
 'proj_y',
 ])
# Number of photons is saved in meta['EXPOSURE'] with this setup
# so no need to waste disk space for photons with probability=0
photons = photons[photons['probability'] > 0]

photons.write(filename, overwrite=True)
print(f'{time.ctime()} - end simulation for {wave[args.i]:05.2f} Ang')