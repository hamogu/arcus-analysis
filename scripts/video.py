'''Copy and paste in IPython for a quick 3D plot'''
import os
import collections
import functools
import numpy as np
from astropy.table import Table
from scipy.interpolate import interp1d
from arcus.defaults import DefaultSource, DefaultPointing
import astropy
import astropy.units as u

from mayavi import mlab
from marxs.visualization.mayavi import plot_object, plot_rays

%matplotlib

from arcus.arcus import ArcusForPlot, Arcus
arc = ArcusForPlot()
fig = mlab.figure(size=(800, 640))
arc.elements[0].display['shape'] = 'None'
outinstrum = plot_object(arc, viewer=fig)


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

photons = mysource.generate_photons(5e3)
photons = mypointing(photons)
photons = arc(photons)
photons['wave'] = photons['energy'].to(u.Angstrom, equivalencies=u.spectral())
ind = (photons['probability'] > 0) & (photons['facet'] >=0) & np.isfinite(photons['detpix_x'])
posdat = arc.postprocess_steps[0].format_positions(atol=1.)[ind, :, :]

#rays = plot_rays(posdat, scalar=photons['energy'][ind],
#                 kwargssurface={'opacity': .5,
#                                'line_width': 1,
#                                 'colormap': 'blue-red'})

animpoints = '''
comment time azimuth elevation distance foc_x foc_y foc_z roll
"overview"                  0. 118. 52. 17000. 280. -750. 7000. 183
"Zoom on FA"                3. 134. 65.  3000.   0 0 12000 166
"FA Rotation 0"             5. 134. 85.  3000.   0 0 12000 166
"FA Rotation 1"             7. 224. 90.  3000.   0 0 12000 166
"FA Rotation 3"            11. 120. 90.  3000.   0 0 12000 166
"FA Rotation 4"            12. 134. 90.  3000.   0 0 12000 166
"Begin fly down"           14. 134.  5.  3000.   0 0 12000 166
"View of both detectors"   19. 134.  5.   160. 420. -10 430 166
"Close up on one detector" 23. 180. 11.    80. 430 -10 100 166
"Close up on one detector" 28. 180. 11.    80. 430 -10 100 166
'''

animtab = Table.read(animpoints, format='ascii', guess=False)

def populate_anim_viewpoint_table(tab, fps=24):
    time = np.arange(tab['time'][0], tab['time'][-1], 1./fps)
    # would be enough to save the interpolation functions, but this is
    # probably OK, too, for short animations.
    out = Table()
    for col in ['azimuth', 'elevation', 'distance', 'roll', 'foc_x', 'foc_y', 'foc_z', 'time']:
        inter = interp1d(tab['time'], tab[col])
        out[col] = inter(time)
    xyz = ['foc_x', 'foc_y', 'foc_z']
    out['focalpoint'] = np.vstack(list([out[n] for n in xyz])).T

    return out

tab = populate_anim_viewpoint_table(animtab, fps=35)

def set_all_in_list(listin, attribute, value):
    for l in listin:
        if l is None:
            pass
        elif isinstance(l, collections.Iterable):
            set_all_in_list(l, attribute, value)
        else:
            setattr(l, attribute, value)

ccd = photons['CCD'][ind]

@mlab.animate(delay=int(1000./24))
def anim(animtab):
    # legendshown = False
    hide_FA = False
    start = 0
    for i, row in enumerate(animtab):
        print(row['time'])
        # Add one ray every other frame
        if (row['time'] < 23.) & (np.mod(i, 2) == 0):
            end = start + 1
            x = posdat[start, :, 0]
            y = posdat[start, :, 1]
            z = posdat[start, :, 2]
            #s = np.ones_like(x) * photons['wave'][start]
            out = mlab.plot3d(x, y, z, color=(.99,.99,.99),
                              #vmin=8,
                              #vmax=50,
                              tube_radius=None,
                              line_width=1)

        elif (row['time'] >= 23.):
            end = start + 120
            # If multiple things are plotted, don't render in between
            fig.scene.disable_render = True
            for i in range(start, end):
                # At t > 23 we are zoomed in. Don't plot photon that will
                # be outside the FOV anyway
                if (ccd[i] < 11.5) | (ccd[i] > 14.5):
                    continue
                    # Don't plot top of lines which are not in the FOv anyway
                x = posdat[i, 2:, 0]
                y = posdat[i, 2:, 1]
                z = posdat[i, 2:, 2]
                # s = np.ones_like(x) * photons['wave'][i]
                out = mlab.plot3d(x, y, z, color=(.99,.99,.99),
                                  #vmin=8,
                                  #vmax=50,
                                  tube_radius=None,
                                  line_width=1.)
            fig.scene.disable_render = False
            fig.scene.render()

        start = end
        # if legendshown is False:
        #     mlab.scalarbar(object=out, title='wavelength [Angstrom]')
        #     legendshown = True
        if (row['time'] > 18) & (hide_FA is False):
            set_all_in_list(outinstrum[:3], 'visible', False)
        mlab.view(azimuth=row['azimuth'], elevation=row['elevation'],
                  distance=row['distance'], roll=row['roll'],
                  focalpoint=row['focalpoint'])
        yield

a = anim(tab) # Starts the animation.
