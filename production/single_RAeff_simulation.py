from os.path import join
import time
import argparse
import numpy as np
import astropy.units as u
from marxs.missions.arcus.arcus import (Arcus, PerfectArcus)
from marxs.missions.arcus.defaults import DefaultSource, DefaultPointing
#from arcus.analysis.grid import aeffRfromraygrid, csv_per_order

parser = argparse.ArgumentParser(
                    description = 'Run a single simulation in a wavelength grid that can later we used to generate R/Aeff curves or even ARFs and RMFs.\n' + 
                       'The interface is specifically designed to make is easy to iterate over this by changing just the "i" parameter' +
                       'which helps with scheduling systems list slurm job arrays.',
                    epilog = 'Text at the bottom of help')
parser.add_argument('outpath', help='Path to write simulation output') 
parser.add_argument('i', type=int, help='index in wavelength array to select wavelength for this simulation') 
parser.add_argument('-n', '--n_photons', default=100000, type=int,
                    help='number of photons')
parser.add_argument('--wave_lo', default=1.5, type=float,
                    help='lower end of wavelength range in Ang')
parser.add_argument('--wave_hi', default=50., type=float,
                    help='upper end of wavelength range in Ang')
parser.add_argument('--wave_step', default=0.15, type=float,
                    help='step size')
parser.add_argument('--perfect',
                    action='store_true',
                    help='Use perfect Arcus without default misalignments') 
args = parser.parse_args()

wave = np.arange(args.wave_lo, args.wave_hi, args.wave_step) * u.Angstrom
energies = wave.to(u.keV, equivalencies=u.spectral())
pointing = DefaultPointing()

instrument = PerfectArcus() if args.perfect else Arcus()
print(f'{time.ctime()} - begin simulation for {wave[args.i]:05.2f} Ang')
source = DefaultSource(energy=energies[args.i])
photons = source.generate_photons(args.n_photons * u.s)
photons = pointing(photons)
photons = instrument(photons)

# There is no way we need all those columns. Waste of space.
photons.keep_columns(
 [
 'energy',   # Could turn into a keyword since it's the same for each row
 'probability',
 'xou',
 'facet',  # Probably don't need, but keep for now
 'order',
 'order_L1',
 'circ_phi',
 'circ_y',
 'CCD',
 'proj_x',
 'proj_y',
 ])
photons.write(join(args.outpath,
                   'wave{0:05.2f}.fits'.format(wave.value[args.i])))
print(f'{time.ctime()} - end simulation for {wave[args.i]:05.2f} Ang')