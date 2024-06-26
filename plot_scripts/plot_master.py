from __future__ import print_function

import time
import os
import json
import argparse

import numpy as np
import astropy.units as u

from astropy.table import Table

from mayavi import mlab
from marxs import simulator
from marxs.visualization.mayavi import plot_object, plot_rays

from arcus import arcus
from arcus import boom
from arcus.defaults import DefaultSource, DefaultPointing


class ArcusPlot(object):
    n_photons = 10000
    wave = np.arange(8., 50., 0.5) * u.Angstrom
    energies = wave.to(u.keV, equivalencies=u.spectral()).value
    instrument = arcus.ArcusForPlot()

    @property
    def filename(self):
        return self.__class__.__name__.lower()

    @property
    def source(self):
        return DefaultSource(energy={'energy': self.energies[::-1],
                                      'flux': np.ones_like(self.energies) / self.energies ** 2}
        )

    pointing = DefaultPointing()

    def __init__(self, outpath):
        self.outpath = outpath

    def simulation(self):
        self.keeppos = simulator.KeepCol('pos')
        self.keepprob = simulator.KeepCol('probability')
        self.instrument.postprocess_steps = [self.keeppos, self.keepprob]
        source = self.source
        self.photons = source.generate_photons(self.n_photons)
        self.photons = self.pointing(self.photons)
        self.photons = self.instrument(self.photons)

    def get_filename(self, index):
        '''Can be overwritten for classes that generate multiple outputs'''
        return self.filename

    def get_jsondata(self, index):
        '''Can be overwritten for classes that generate multiple outputs'''
        return self.jsondata()

    def jsondata(self):
        self.data['path'] = os.path.join(self.outpath, self.filename)
        return self.data

    def save_html_data(self, index=None):
        '''Save text like page title etc. that is defined in class.

        This writes a json file at the same location as the other output.
        The script that makes the website can then parse that json file and
        fill in the values in the website template.
        This way, the code that generates the figure and the text for the
        caption are in the same place and the website generator does not have
        to read in this class again.

        Parameters
        ----------
        index : int
            Index for the plot number. Only used for classes that generate
            more than one plot (in most cases each class is responsible for
            exactly one plot)
        '''
        obj = self.get_jsondata(index)
        filename = os.path.join(self.outpath, self.get_filename(index) + '.json')
        with open(filename, 'w') as f:
            json.dump(obj, f, indent=2)


class X3d(ArcusPlot):
    pass


class XBasicFlat(X3d):
    plot_col_color = 'energy'

    data = {'name': 'xbasicflat',
            'caption': 'Ray-trace components',
            'title': '',
            'figcaption': '''
            <p>
            Basic design for ARCUS, showing a ray-trace for a simple point
            source at infinity with a flat spectrum. ARCUS has four channels
            positioned in two pairs. Each pair of
            channels shares an optical axis.
            The two optical axes are shown here as thin blue
            cylinders.
            </p>

            <p>
            Rays start in the aperture, which consists of four rectangles
            located above the SPO channels. Each channel has a number of SPO
            modules, shown in green.
            The modules have different dimensions depending on their distance
            from their respective optical axis.
            Photons bounce of their mirrors twice in a Wolter type I like
            geometry. However, in this simulation the SPOs are somewhat
            simplified such that the
            reflection actually happens in a single plane, shown in white.
            Rays are imaged onto detectors (yellow).
            </p>
            '''
            }

    def plot(self):
        ind = (self.photons['probability'] > 0)
        posdat = self.keeppos.format_positions()[ind, :, :]
        fig = mlab.figure()
        obj = plot_object(self.instrument, viewer=fig)
        rays = plot_rays(posdat, scalar=self.photons[self.plot_col_color][ind])
        mlab.savefig(os.path.join(self.outpath, self.filename + '.x3d'))


class XSingleEnergy(XBasicFlat):
    plot_col_color = 'order'

    data = {'name': 'singlenergy',
            'caption': 'Single energy',
            'title': '',
            'figcaption': '''
            <p>
            Basic design for ARCUS, showing a ray-trace for an on-axis point
            source with a monochromatic emission at 0.5 keV = 24.8 Angstrom.
            ARCUS has four channels
            positioned in two pairs. Each pair of
            channels shares an optical axis.
            The two optical axes are shown here as thin blue
            cylinders.
            </p>

            <p>
            Rays start in the aperture, which consists of four rectangles
            located above the SPO channels. Each channel has a number of SPO
            modules, shown in green.
            The modules have different dimensions depending on their distance
            from their respective optical axis.
            Photons bounce of their mirrors twice in a Wolter type I like
            geometry. However, in this simulation the SPOs are somewhat
            simplified such that the
            reflection actually happens in a single plane, shown in white.
            Rays are imaged onto detectors (yellow).
            </p>

            <p>
            A closer look at the detectors shows that several
            grating orders are detected at the same time.
            </p>

            '''}

    @property
    def source(self):
        return DefaultSource(energy=0.5, flux=1.)


class XEQPegA(XBasicFlat):
    plot_col_color = 'order'
    n_photons = 2000

    data = {'name': 'emissionlines',
            'caption': 'Emission Line',
            'title': '',
            'figcaption': '''
            <p>
            This simulation shows rays from a stellar source. The spectrum
            consists of a continuum with strong emission lines
            superimposed. The particular spectrum for this simulation
            is taken from the active star EQ Peg A, but the spectrum is
            really representative of any collisionally excited plasma.
            </p>

            <p>
            Rays start in the aperture, which consists of four rectangles
            located above the SPO channels. Each channel has a number of SPO
            modules, shown in green.
            The modules have different dimensions depending on their distance
            from their respective optical axis.
            Photons bounce of their mirrors twice in a Wolter type I like
            geometry. However, in this simulation the SPOs are somewhat
            simplified such that the
            reflection actually happens in a single plane, shown in white.
            Rays are imaged onto detectors (yellow).
            </p>

            <p>
            The coloring of the rays reflects the grating order that each ray
            is dispersed into. Because the continuum covers a wide energy range,
            almost any point on the detector receives photons from several
            different orders. However, these photons have different energies
            and can later be sorted apart based on the CCD energy resolution.
            Every individual emission line is typically seen in 2-3 order
            at different positions on the detector. For most lines, only one
            of those orders is strong enough that it really stands out in
            this simulation, because only a limited number of rays is shown
            here for simplicity.
            </p>

            '''}

    @property
    def source(self):
        EQPegAspec = Table.read(os.path.join(os.path.dirname(__file__), '../inputdata/EQPegA_flux.tbl'), format='ascii',
                        names=['energy', 'flux'])
        # restrict table to ARCUS energy range
        EQPegAspec = EQPegAspec[(EQPegAspec['energy'] > 0.25) &
                                (EQPegAspec['energy'] < 1.5)]

        # define position and spectrum of source
        mysource = DefaultSource(energy=EQPegAspec,
                                 geomarea=self.instrument.elements[0].area,
                                 flux=(EQPegAspec['flux'][1:] * np.diff(EQPegAspec['energy'])).sum())
        return mysource


class Boom(XBasicFlat):
    plot_col_color = 'hitrod'

    data = {'name': 'boom',
            'caption': 'Absorption by the boom',
            'title': '',
            'figcaption': '''
            <p>
            This simulation uses a input spectrum with a flat spectrum
            (flat in wavelength space).
            The elements of the 3-sided boom intersect some of the photons.
            The ray path are colored to distinguish photons that make
            it to the detector (gray) and photons that hit an element of
            the boom at some point (red). On the top and bottom end the
            boom will be attached to some type of socket, which is not
            shown here.
            </p>

            <p>
            This figure is meant to display the shape of the boom. The
            related notebook goes in more detail to compare booms
            of different sizes and rotation angles.
            </p>
            '''}

    def simulation(self):
        super(Boom, self).simulation()
        self.myboom = boom.ThreeSidedBoom(orientation=arcus.xyz2zxy[:3, :3],
                                          position=[arcus.defaultconf['d'], 0, 1000])
        self.photons = self.myboom(self.photons)

    def plot(self):
        self.myboom.elements[0].display['color'] = 'blue'
        ind = (self.photons['probability'] > 0)
        posdat = self.keeppos.format_positions()[ind, :, :]
        fig = mlab.figure()
        obj = plot_object(self.instrument, viewer=fig)
        obj = plot_object(self.myboom, viewer=fig)
        rays = plot_rays(posdat, scalar=np.asarray(self.photons[self.plot_col_color], dtype=float)[ind])
        mlab.savefig(os.path.join(self.outpath, self.filename + '.x3d'))


all_plots = ['XSingleEnergy', 'XEQPegA', 'XSingleEnergy']

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
    Run plotting scripts to make output for website.
    ''')
    parser.add_argument('outpath',
                        help='base directory for output')
    parser.add_argument('-p', '--plot', nargs='+',
                        help='name specific plots that should be generated')
    args = parser.parse_args()

    plot = args.plot if args.plot is not None else all_plots
    for plotname in plot:
        a = globals()[plotname](args.outpath)
        print('Running {}'.format(a.filename))
        a.simulation()
        a.plot()
        a.save_html_data()
    # make sure all files are written before the window closes at the
    # end of the script
    time.sleep(1)
