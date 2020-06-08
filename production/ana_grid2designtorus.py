from __future__ import print_function
import time
import os
from glob import glob
import numpy as np
from astropy import table
from astropy.io import fits
from arcus.analyze_design import summarize_file, get_wave, load_prepare_table

from utils import get_path
outpath = get_path('grid2designtorus')
outfile = os.path.join(outpath, 'summary.fits')

orders = np.arange(-15, 4)
filelist = glob(os.path.join(outpath, '*.fits'))

filelist.sort()
outtab = None

wave0 = get_wave(load_prepare_table(filelist[0]))

for f in filelist:
    # Output file from previous run - skip
    if 'summary.fits' in f:
        continue
    print('{} - {}'.format(f, time.ctime()))
    out, wave = summarize_file(f, orders)
    out = table.Table(out)
    if not np.allclose(wave0.value, wave.value):
        raise Exception('Wavegrid not consistent between files in grid {}'.format(f))
    if outtab is None:
        outtab = out
        wavetab = table.Table([wave], names=['wave'])
        wavetab['wave'].unit = wave.unit
        ordertab = table.Table([orders])
    else:
        outtab = table.vstack([outtab, out], join_type='exact')
    # write after every step
    # - can start analysis while this script is still running
    # - save parts in case of crash
    # - Maybe hook up some clever multi-tasking here?
    outtab.write(outfile, overwrite=True)
    # Really it would be enough to do that once, but it does not take long
    # and it means I have those extensions available any time.
    # takes about 1 s per run...
    hdulist = fits.open(outfile, 'update')
    hdulist.append(fits.table_to_hdu(wavetab))
    hdulist.append(fits.table_to_hdu(ordertab))
    hdulist.flush()
    hdulist.close()
