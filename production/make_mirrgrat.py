import datetime
import os
import glob

import numpy as np
from astropy.table import QTable
from astropy.io import fits
import astropy.units as u

from marxs.missions.arcus.utils import config
from marxs.missions.arcus.arfrmf import filtersqe

wavefile = {}
for f in glob.glob('wave*'):
    wavefile[fits.getval(f, 'wavelen', ext=('EVENTS', 1))] = f

wave = list(wavefile.keys())
wave.sort()

# For now, just fake it anyway
widths = [5] * u.arcmin


orders0 = np.arange(-11, 3)
orders = orders0[orders0 != 0]


aefforder = QTable({'wave': wave * u.Angstrom})
for o in orders0:
    aefforder[str(o)] = np.zeros(len(wave)) * u.cm**2

for iwav, wav in enumerate(wave):
    evt = QTable.read(wavefile[wav])
    # The simulations are run with filters and the CCD QE in place, but
    # we don't want the to be in the mirr_grat.tab in the end
    # as mkarf applies that separately.
    # So, need to divide out here.
    filtqe = filtersqe(([wav] * u.Angstrom).to(u.keV, equivalencies=u.spectral()))
    aeffexp = evt.meta['A_GEOM'] * u.cm**2 / evt.meta['EXPOSURE']
    for o in orders0:
        ind = evt['order'] == o
        aefforder[str(o)][iwav] = aeffexp * evt['probability'][ind].sum() / filtqe


aefforder.meta = evt.meta.copy()
# remove keys that make no sense for aggregate table
# This uses an exclude-list rather than an include-list to include
# any generic keywords that are added in the future
for k in ['EXTNAME', 'ENERGY', 'WAVELEN']:
    del aefforder.meta[k]

aefforder.meta['DATE'] = datetime.datetime.now().isoformat()[:10]
#aefforder.write('mirr_grat.tab', format='ascii.ecsv', overwrite=True)

aefforder.write(os.path.join(config['data']['caldb_inputdata'], 'aeff',
                             'mirr_grat.tab'),
                format='ascii.ecsv', overwrite=True)
