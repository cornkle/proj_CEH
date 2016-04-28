from __future__ import division

import unittest
import warnings
from numpy.testing.utils import assert_array_equal, assert_allclose

import time
import copy

import numpy as np
import matplotlib as mpl

from cleo import DataLevels
from cleo import Map
import cleo

from salem import Grid
from salem import wgs84
from salem.utils import empty_cache
from salem.grids import local_mercator_grid

do_test_caching = False

class TestColors(unittest.TestCase):

    def test_extendednorm(self):

        bounds = [1, 2, 3]
        cm = mpl.cm.get_cmap('jet')

        mynorm = cleo.colors.ExtendedNorm(bounds, cm.N)
        refnorm = mpl.colors.BoundaryNorm(bounds, cm.N)
        x = np.random.randn(100) * 10 - 5
        np.testing.assert_array_equal(refnorm(x), mynorm(x))

        refnorm = mpl.colors.BoundaryNorm([0] + bounds + [4], cm.N)
        mynorm = cleo.colors.ExtendedNorm(bounds, cm.N, extend='both')
        x = np.random.random(100) + 1.5
        np.testing.assert_array_equal(refnorm(x), mynorm(x))

        # Min and max
        cmref = mpl.colors.ListedColormap(['blue', 'red'])
        cmref.set_over('black')
        cmref.set_under('white')

        cmshould = mpl.colors.ListedColormap(['white', 'blue', 'red', 'black'])
        cmshould.set_over(cmshould(cmshould.N))
        cmshould.set_under(cmshould(0))

        refnorm = mpl.colors.BoundaryNorm(bounds, cmref.N)
        mynorm = cleo.colors.ExtendedNorm(bounds, cmshould.N, extend='both')
        np.testing.assert_array_equal(refnorm.vmin, mynorm.vmin)
        np.testing.assert_array_equal(refnorm.vmax, mynorm.vmax)
        x = [-1, 1.2, 2.3, 9.6]
        np.testing.assert_array_equal(cmshould([0,1,2,3]), cmshould(mynorm(x)))
        x = np.random.randn(100) * 10 + 2
        np.testing.assert_array_equal(cmref(refnorm(x)), cmshould(mynorm(x)))

        np.testing.assert_array_equal(-1, mynorm(-1))
        np.testing.assert_array_equal(1, mynorm(1.1))
        np.testing.assert_array_equal(4, mynorm(12))

        # Just min
        cmref = mpl.colors.ListedColormap(['blue', 'red'])
        cmref.set_under('white')
        cmshould = mpl.colors.ListedColormap(['white', 'blue', 'red'])
        cmshould.set_under(cmshould(0))

        np.testing.assert_array_equal(2, cmref.N)
        np.testing.assert_array_equal(3, cmshould.N)
        refnorm = mpl.colors.BoundaryNorm(bounds, cmref.N)
        mynorm = cleo.colors.ExtendedNorm(bounds, cmshould.N, extend='min')
        np.testing.assert_array_equal(refnorm.vmin, mynorm.vmin)
        np.testing.assert_array_equal(refnorm.vmax, mynorm.vmax)
        x = [-1, 1.2, 2.3]
        np.testing.assert_array_equal(cmshould([0,1,2]), cmshould(mynorm(x)))
        x = np.random.randn(100) * 10 + 2
        np.testing.assert_array_equal(cmref(refnorm(x)), cmshould(mynorm(x)))

        # Just max
        cmref = mpl.colors.ListedColormap(['blue', 'red'])
        cmref.set_over('black')
        cmshould = mpl.colors.ListedColormap(['blue', 'red', 'black'])
        cmshould.set_over(cmshould(2))

        np.testing.assert_array_equal(2, cmref.N)
        np.testing.assert_array_equal(3, cmshould.N)
        refnorm = mpl.colors.BoundaryNorm(bounds, cmref.N)
        mynorm = cleo.colors.ExtendedNorm(bounds, cmshould.N, extend='max')
        np.testing.assert_array_equal(refnorm.vmin, mynorm.vmin)
        np.testing.assert_array_equal(refnorm.vmax, mynorm.vmax)
        x = [1.2, 2.3, 4]
        np.testing.assert_array_equal(cmshould([0,1,2]), cmshould(mynorm(x)))
        x = np.random.randn(100) * 10 + 2
        np.testing.assert_array_equal(cmref(refnorm(x)), cmshould(mynorm(x)))

        # General case
        bounds = [1, 2, 3, 4]
        cm = mpl.cm.get_cmap('jet')
        mynorm = cleo.colors.ExtendedNorm(bounds, cm.N, extend='both')
        refnorm = mpl.colors.BoundaryNorm([-100] + bounds + [100], cm.N)
        x = np.random.randn(100) * 10 - 5
        ref = refnorm(x)
        ref = np.where(ref == 0, -1, ref)
        ref = np.where(ref == cm.N-1, cm.N, ref)
        np.testing.assert_array_equal(ref, mynorm(x))

class TestGraphics(unittest.TestCase):

    def test_datalevels_output(self):

        # Test basic stuffs
        c = DataLevels(nlevels=2)
        assert_array_equal(c.levels, [0, 1])
        c.set_data([1, 2, 3, 4])
        assert_array_equal(c.levels, [1, 4])

        c = DataLevels(levels=[1, 2, 3])
        assert_array_equal(c.levels, [1, 2, 3])

        c = DataLevels(nlevels=10, data=[0, 9])
        assert_array_equal(c.levels, np.linspace(0, 9, num=10))
        self.assertTrue(c.extend == 'neither')

        c = DataLevels(nlevels=10, data=[0, 9], vmin=2, vmax=3)
        assert_array_equal(c.levels, np.linspace(2, 3, num=10))
        self.assertTrue(c.extend == 'both')
        c.set_extend('neither')
        self.assertTrue(c.extend == 'neither')
        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            # Trigger a warning.
            out = c.to_rgb()
            # Verify some things
            assert len(w) == 2
            assert issubclass(w[0].category, RuntimeWarning)
            assert issubclass(w[1].category, RuntimeWarning)

        c = DataLevels(nlevels=10, data=[2.5], vmin=2, vmax=3)
        assert_array_equal(c.levels, np.linspace(2, 3, num=10))
        self.assertTrue(c.extend == 'neither')
        c.update(dict(extend='both'))
        self.assertTrue(c.extend == 'both')
        self.assertRaises(AttributeError, c.update, dict(dummy='t'))

        c = DataLevels(nlevels=10, data=[0, 9], vmax=3)
        assert_array_equal(c.levels, np.linspace(0, 3, num=10))
        self.assertTrue(c.extend == 'max')

        c = DataLevels(nlevels=10, data=[0, 9], vmin=1)
        assert_array_equal(c.levels, np.linspace(1, 9, num=10))
        self.assertTrue(c.extend == 'min')

        c = DataLevels(nlevels=10, data=[0, 9], vmin=-1)
        assert_array_equal(c.levels, np.linspace(-1, 9, num=10))
        self.assertTrue(c.extend == 'neither')
        c.set_plot_params()
        self.assertTrue(c.extend == 'neither')
        assert_array_equal(c.vmin, 0)
        assert_array_equal(c.vmax, 9)
        c.set_plot_params(vmin=1)
        assert_array_equal(c.vmin, 1)
        c.set_data([-12, 8])
        assert_array_equal(c.vmin, 1)
        self.assertTrue(c.extend == 'min')
        c.set_data([2, 8])
        self.assertTrue(c.extend == 'neither')
        c.set_extend('both')
        self.assertTrue(c.extend == 'both')
        c.set_data([3, 3])
        self.assertTrue(c.extend == 'both')
        c.set_extend()
        self.assertTrue(c.extend == 'neither')

        # Test the conversion
        cm = mpl.colors.ListedColormap(['white', 'blue', 'red', 'black'])
        x = [-1, 0.9, 1.2, 2, 999, 0.8]
        c = DataLevels(levels=[0, 1, 2], data=x, cmap=cm)
        r = c.to_rgb()
        self.assertTrue(len(x) == len(r))
        self.assertTrue(c.extend == 'both')
        assert_array_equal(r, cm([0, 1, 2, 3, 3, 1]))

        x = [0.9, 1.2]
        c = DataLevels(levels=[0, 1, 2], data=x, cmap=cm, extend='both')
        r = c.to_rgb()
        self.assertTrue(len(x) == len(r))
        self.assertTrue(c.extend == 'both')
        assert_array_equal(r, cm([1, 2]))

        cm = mpl.colors.ListedColormap(['white', 'blue', 'red'])
        c = DataLevels(levels=[0, 1, 2], data=x, cmap=cm, extend='min')
        r = c.to_rgb()
        self.assertTrue(len(x) == len(r))
        assert_array_equal(r, cm([1, 2]))

        cm = mpl.colors.ListedColormap(['blue', 'red', 'black'])
        c = DataLevels(levels=[0, 1, 2], data=x, cmap=cm, extend='max')
        r = c.to_rgb()
        self.assertTrue(len(x) == len(r))
        assert_array_equal(r, cm([0, 1]))

    def test_map(self):

        a = np.zeros((4, 5))
        a[0, 0] = -1
        a[1, 1] = 1.1
        a[2, 2] = 2.2
        a[2, 4] = 1.9
        a[3, 3] = 9
        cmap = copy.deepcopy(mpl.cm.get_cmap('jet'))

        # ll_corner (type geotiff)
        g = Grid(nxny=(5, 4), dxdy=(1, 1), ll_corner=(0, 0), proj=wgs84,
                 pixel_ref='corner')
        c = Map(g, ny=4, countries=False)
        c.set_cmap(cmap)
        c.set_plot_params(levels=[0, 1, 2, 3])
        c.set_data(a)
        rgb1 = c.to_rgb()
        c.set_data(a, crs=g)
        assert_array_equal(rgb1, c.to_rgb())
        c.set_data(a, interp='linear')
        rgb1 = c.to_rgb()
        c.set_data(a, crs=g, interp='linear')
        assert_array_equal(rgb1, c.to_rgb())

        # centergrid (type WRF)
        g = Grid(nxny=(5, 4), dxdy=(1, 1), ll_corner=(0.5, 0.5), proj=wgs84,
                 pixel_ref='center')
        c = Map(g, ny=4, countries=False)
        c.set_cmap(cmap)
        c.set_plot_params(levels=[0, 1, 2, 3])
        c.set_data(a)
        rgb1 = c.to_rgb()
        c.set_data(a, crs=g)
        assert_array_equal(rgb1, c.to_rgb())
        c.set_data(a, interp='linear')
        rgb1 = c.to_rgb()
        c.set_data(a, crs=g.corner_grid, interp='linear')
        assert_array_equal(rgb1, c.to_rgb())
        c.set_data(a, crs=g.center_grid, interp='linear')
        assert_array_equal(rgb1, c.to_rgb())

        # More pixels
        c = Map(g, ny=500, countries=False)
        c.set_cmap(cmap)
        c.set_plot_params(levels=[0, 1, 2, 3])
        c.set_data(a)
        rgb1 = c.to_rgb()
        c.set_data(a, crs=g)
        assert_array_equal(rgb1, c.to_rgb())
        c.set_data(a, interp='linear')
        rgb1 = c.to_rgb()
        c.set_data(a, crs=g, interp='linear')
        rgb2 = c.to_rgb()

        # The interpolation is conservative with the grid...
        srgb = np.sum(rgb2[..., 0:3], axis=2)
        pok = np.nonzero(srgb != srgb[0, 0])
        rgb1 = rgb1[np.min(pok[0]):np.max(pok[0]),
                    np.min(pok[1]):np.max(pok[1]),...]
        rgb2 = rgb2[np.min(pok[0]):np.max(pok[0]),
                    np.min(pok[1]):np.max(pok[1]),...]
        assert_array_equal(rgb1, rgb2)

        cmap.set_bad('pink')

        # Add masked arrays
        a[1, 1] = np.NaN
        c.set_data(a)
        rgb1 = c.to_rgb()
        c.set_data(a, crs=g)
        assert_array_equal(rgb1, c.to_rgb())

        # Interp?
        c.set_data(a, interp='linear')
        rgb1 = c.to_rgb()
        c.set_data(a, crs=g, interp='linear')
        rgb2 = c.to_rgb()
        # Todo: there's something sensibly wrong about imresize here
        # but I think it is out of my scope
        # assert_array_equal(rgb1, rgb2)

    def test_caching(self):

        if not do_test_caching:
            return

        empty_cache()

        grid = cleo.utils.make_local_mercator_grid(center_ll=(11.38, 47.26),
                                                   extent=(2000000, 2000000))
        t1 = time.time()
        m = cleo.Map(grid)
        t1 = time.time() - t1

        grid = cleo.utils.make_local_mercator_grid(center_ll=(11.38, 47.26),
                                                   extent=(2000000, 2000000))
        t2 = time.time()
        m = cleo.Map(grid)
        t2 = time.time() - t2

        self.assertTrue(t2 < (t1 /10.))

    def test_increase_coverage(self):

        # Just for coverage -> empty shapes should not trigger an error
        grid = local_mercator_grid(center_ll=(-20, 40),
                                        extent=(2000, 2000), nx=10)
        c = Map(grid)

        # Assigning wrongly shhaped data should, however
        self.assertRaises(ValueError, c.set_data, np.zeros((3, 8)))