'''Copy and paste in IPython for a quick 3D plot'''
import os
import numpy as np
from astropy.table import Table
from scipy.interpolate import interp1d
from arcus.defaults import DefaultSource, DefaultPointing
import astropy

from mayavi import mlab
from marxs.visualization.mayavi import plot_object, plot_rays

%matplotlib

from arcus.arcus import ArcusForPlot, Arcus
arc = ArcusForPlot()
fig = mlab.figure()
out = plot_object(arc, viewer=fig)


EQPegAspec = Table.read('../inputdata/EQPegA_flux.tbl', format='ascii',
                        names=['energy', 'flux'])
# restrict table to ARCUS energy range
EQPegAspec = EQPegAspec[(EQPegAspec['energy'] > 0.25) &
                        (EQPegAspec['energy'] < 1.5)]

coord = astropy.coordinates.SkyCoord.from_name("EQ Peg A")

mysource = DefaultSource(coords=coord, energy=EQPegAspec,
                         geomarea=arc.elements[0].area,
                         flux=(EQPegAspec['flux'][1:] * np.diff(EQPegAspec['energy'])).sum())
mypointing = DefaultPointing(coords=coord)

photons = mysource.generate_photons(2e3)
photons = mypointing(photons)
photons = arc(photons)

ind = (photons['probability'] > 0) & (photons['facet'] >=0) & np.isfinite(photons['detpix_x'])
posdat = arc.postprocess_steps[0].format_positions(atol=1.)[ind, :, :]

rays = plot_rays(posdat, scalar=photons['energy'][ind],
                 kwargssurface={'opacity': .5,
                                'line_width': 1,
                                 'colormap': 'blue-red'})

animpoints = '''
comment time azimuth elevation distance foc_x foc_y foc_z roll
"overview"                  0. 118. 52. 17000. 280. -750. 7000. 183
"Zoom on FA"                3. 134. 65.  3000.   0 0 12000 166
"FA Rotation 0"             5. 134. 85.  3000.   0 0 12000 166
"FA Rotation 1"             7. 224. 90.  3000.   0 0 12000 166
"FA Rotation 3"            11. 120. 90.  3000.   0 0 12000 166
"FA Rotation 4"            13. 134. 90.  3000.   0 0 12000 166
"Begin fly down"           15. 134.  5.  3000.   0 0 12000 166
"View of both detectors"   20. 134.  5.   160. 420. -10 430 166
"Close up on one detector" 23. 180. 11.    80. 430 -10 100 166
'''

animtab = Table.read(animpoints, format='ascii', guess=False)

def populate_anim_viewpoint_table(tab, fps=24):
    time = np.arange(tab['time'][0], tab['time'][-1], 1./fps)
    # would be enough to save the interpolation functions, but this is
    # probably OK, too, for short animations.
    out = Table()
    for col in ['azimuth', 'elevation', 'distance', 'roll', 'foc_x', 'foc_y', 'foc_z']:
        inter = interp1d(tab['time'], tab[col])
        out[col] = inter(time)
    xyz = ['foc_x', 'foc_y', 'foc_z']
    out['focalpoint'] = np.vstack(list([out[n] for n in xyz])).T
    out.remove_columns(xyz)
    return out

tab = populate_anim_viewpoint_table(animtab)

@mlab.animate(delay=int(1000./24))
def anim(animtab):
    for row in animtab:
        mlab.view(**{c:r for c, r in zip(row.colnames, row)})
        fig.scene.render()
        print(row['azimuth'])
        yield

a = anim(tab) # Starts the animation.
