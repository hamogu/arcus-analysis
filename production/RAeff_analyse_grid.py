import numpy as np
from marxs.missions.arcus.analyze_grid import aeffRfromraygrid
from marxs.missions.arcus.arcus import defaultconf as conf

out = aeffRfromraygrid('raysRAeffperfect/',
                       conf,
                       orders=np.arange(-20, 5),
                      allow_inconsistent_meta=False)
out.write('raygrid-perfectRAeff.fits')

out = aeffRfromraygrid('raysRAeff/',
                       conf,
                       orders=np.arange(-20, 5),
                      allow_inconsistent_meta=False)
out.write('raygridRAeff.fits')