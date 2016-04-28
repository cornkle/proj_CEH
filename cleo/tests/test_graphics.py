from __future__ import division

import warnings
from nose.tools import assert_true, assert_raises

import os
import copy
import shutil

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.testing.decorators import image_comparison
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import shapely.geometry as shpg
import geopandas as gpd
import nose

from cleo import DataLevels
from cleo import Map
from cleo import files
import cleo

import salem
from salem import Grid
from salem import wgs84
from salem.utils import get_demo_file
from salem.datasets import GeoNetcdf, GoogleCenterMap, GoogleVisibleMap

# Globals
current_dir = os.path.dirname(os.path.abspath(__file__))
testdir = os.path.join(current_dir, 'tmp')

def _create_dummy_shp(fname):

    if not os.path.exists(testdir):
        os.makedirs(testdir)

    e_line = shpg.LinearRing([(1.5, 1), (2., 1.5), (1.5, 2.), (1, 1.5)])
    i_line = shpg.LinearRing([(1.4, 1.4), (1.6, 1.4), (1.6, 1.6), (1.4, 1.6)])
    p1 = shpg.Polygon(e_line, [i_line])
    p2 = shpg.Polygon([(2.5, 1.3), (3., 1.8), (2.5, 2.3), (2, 1.8)])
    p3 = shpg.Point(0.5, 0.5)
    p4 = shpg.Point(1, 1)
    df = gpd.GeoDataFrame()
    df['name'] = ['Polygon', 'Line']
    df['geometry'] = gpd.GeoSeries([p1, p2])
    of = os.path.join(testdir, fname)
    df.crs = {'init': 'epsg:4326'}
    df.to_file(of)
    return of


@image_comparison(baseline_images=['test_extendednorm'],
                  extensions=['png'])
def test_extendednorm():

    a = np.zeros((4, 5))
    a[0,0] = -9999
    a[1,1] = 1.1
    a[2,2] = 2.2
    a[2,4] = 1.9
    a[3,3] = 9999999

    cm = mpl.cm.get_cmap('jet')
    bounds = [0,1,2,3]
    norm = cleo.colors.ExtendedNorm(bounds, cm.N, extend='both')

    #fig, (ax1, ax2) = plt.subplots(1, 2)
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    imax = ax1.imshow(a, interpolation='None', norm=norm, cmap=cm,
                     origin='lower');
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.2)
    plt.colorbar(imax, cax=cax, extend='both')

    ti = cm(norm(a))
    ax2.imshow(ti, interpolation='None', origin='lower')
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="5%", pad=0.2)
    cbar = mpl.colorbar.ColorbarBase(cax, extend='both', cmap=cm,
                                     norm=norm)
    fig.tight_layout()


@image_comparison(baseline_images=['test_datalevels'],
                  extensions=['png'])
def test_datalevels():

    plt.close()

    a = np.zeros((4, 5))
    a[0, 0] = -1
    a[1, 1] = 1.1
    a[2, 2] = 2.2
    a[2, 4] = 1.9
    a[3, 3] = 9

    cm = copy.copy(mpl.cm.get_cmap('jet'))
    cm.set_bad('pink')

    # fig, axes = plt.subplots(nrows=3, ncols=2)
    fig = plt.figure()
    ax = iter([fig.add_subplot(3, 2, i) for i in [1,2,3,4,5,6]])

    # The extended version should be automated
    c = DataLevels(levels=[0,1,2,3], data=a, cmap=cm)
    c.visualize(next(ax), title='levels=[0,1,2,3]')

    # Without min
    a[0, 0] = 0
    c = DataLevels(levels=[0,1,2,3], data=a, cmap=cm)
    c.visualize(next(ax), title='modified a for no min oob')

    # Without max
    a[3, 3] = 0
    c = DataLevels(levels=[0,1,2,3], data=a, cmap=cm)
    c.visualize(next(ax), title='modified a for no max oob')

    # Forced bounds
    c = DataLevels(levels=[0,1,2,3], data=a, cmap=cm, extend='both')
    c.visualize(next(ax), title="extend='both'")

    # Autom nlevels
    a[0, 0] = -1
    a[3, 3] = 9
    c = DataLevels(nlevels=127, vmin=0, vmax=3, data=a, cmap=cm)
    c.visualize(next(ax), title="Auto levels with oob data")

    # Missing data
    a[3, 0] = np.NaN
    c = DataLevels(nlevels=127, vmin=0, vmax=3, data=a, cmap=cm)
    c.visualize(next(ax), title="missing data")

    plt.tight_layout()


@image_comparison(baseline_images=['test_datalevels_visu_h',
                                   'test_datalevels_visu_v'],
                  extensions=['png'])
def test_datalevels_visu():

    a = np.array([-1., 0., 1.1, 1.9, 9.])
    cm = mpl.cm.get_cmap('RdYlBu_r')

    dl = DataLevels(a, cmap=cm, levels=[0, 1, 2, 3])
    dl.visualize(orientation='horizontal', add_values=True)

    dl = DataLevels(a.reshape((5,1)), cmap=cm, levels=[0, 1, 2, 3])
    dl.visualize(orientation='vertical', add_values=True)


@image_comparison(baseline_images=['test_simple_map'],
                  extensions=['png'])
def test_simple_map():

    a = np.zeros((4, 5))
    a[0, 0] = -1
    a[1, 1] = 1.1
    a[2, 2] = 2.2
    a[2, 4] = 1.9
    a[3, 3] = 9
    a_inv = a[::-1, :]
    fs = _create_dummy_shp('fs.shp')

    # UL Corner
    g1 = Grid(nxny=(5, 4), dxdy=(1, -1), ul_corner=(-1, 3), proj=wgs84,
             pixel_ref='corner')
    c1 = Map(g1, ny=4, countries=False)

    # LL Corner
    g2 = Grid(nxny=(5, 4), dxdy=(1, 1), ll_corner=(-1, -1), proj=wgs84,
             pixel_ref='corner')
    c2 = Map(g2, ny=4, countries=False)

    # Settings
    for c, data in zip([c1, c2], [a_inv, a]):
        c.set_cmap(mpl.cm.get_cmap('jet'))
        c.set_plot_params(levels=[0, 1, 2, 3])
        c.set_data(data)
        c.set_shapefile(fs)
        c.set_lonlat_countours(interval=0.5)

    fig = plt.figure(figsize=(9, 8))
    ax1 = fig.add_subplot(321)
    ax2 = fig.add_subplot(322)
    c1.visualize(ax1)
    c2.visualize(ax2)

    # UL Corner
    c1 = Map(g1, ny=400, countries=False)
    c2 = Map(g2, ny=400, countries=False)
    # Settings
    for c, data, g in zip([c1, c2], [a_inv, a], [g1, g2]):
        c.set_cmap(mpl.cm.get_cmap('jet'))
        c.set_data(data, crs=g)
        c.set_shapefile(fs)
        c.set_plot_params(nlevels=256)
        c.set_lonlat_countours(interval=2)
    ax1 = fig.add_subplot(323)
    ax2 = fig.add_subplot(324)
    c1.visualize(ax1)
    c2.visualize(ax2)

    # Settings
    for c, data in zip([c1, c2], [a_inv, a]):
        c.set_plot_params(nlevels=256, vmax=3)
        c.set_lonlat_countours(interval=1)
        c.set_data(data, interp='linear')
    ax1 = fig.add_subplot(325)
    ax2 = fig.add_subplot(326)
    c1.visualize(ax1)
    c2.visualize(ax2)

    fig.tight_layout()
    if os.path.exists(testdir):
        shutil.rmtree(testdir)


@image_comparison(baseline_images=['test_merca_map'],
                  extensions=['png'])
def test_merca_map():

    grid = salem.grids.local_mercator_grid(center_ll=(11.38, 47.26),
                                               extent=(2000000, 2000000))

    m1 = Map(grid)

    grid = salem.grids.local_mercator_grid(center_ll=(11.38, 47.26),
                                               extent=(2000000, 2000000),
                                               order='ul')
    m2 = Map(grid)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    m1.visualize(ax=ax1, addcbar=False)
    m2.visualize(ax=ax2, addcbar=False)
    plt.tight_layout()


@image_comparison(baseline_images=['test_oceans'],
                  extensions=['png'])
def test_oceans():

    f = os.path.join(get_demo_file('wrf_tip_d1.nc'))
    grid = GeoNetcdf(f).grid
    m = Map(grid, countries=False)
    m.set_shapefile(rivers=True, linewidths=2)
    m.set_shapefile(oceans=True, edgecolor='k', linewidth=3)

    fig, ax = plt.subplots(1, 1)
    m.visualize(ax=ax, addcbar=False)
    plt.tight_layout()


@image_comparison(baseline_images=['test_geometries'],
                  extensions=['png'])
def test_geometries():

    # UL Corner
    g = Grid(nxny=(5, 4), dxdy=(10, 10), ll_corner=(-20, -15), proj=wgs84,
             pixel_ref='corner')
    c = Map(g, ny=4)
    c.set_lonlat_countours(interval=10., colors='crimson')

    c.set_geometry(shpg.Point(10, 10), color='darkred', markersize=60)
    c.set_geometry(shpg.Point(5, 5), s=500, marker='s',
                   facecolor='green', hatch='||||')

    s = np.array([(-5, -10), (0., -5), (-5, 0.), (-10, -5)])
    l1 = shpg.LineString(s)
    l2 = shpg.LinearRing(s+3)
    c.set_geometry(l1)
    c.set_geometry(l2, color='pink', linewidth=3)

    s += 20
    p = shpg.Polygon(shpg.LineString(s), [shpg.LineString(s/4 + 10)])
    c.set_geometry(p, facecolor='red', edgecolor='k', linewidth=3, alpha=0.5)

    p1 = shpg.Point(20, 10)
    p2 = shpg.Point(20, 20)
    p3 = shpg.Point(10, 20)
    mpoints = shpg.MultiPoint([p1, p2, p3])
    c.set_geometry(mpoints, s=250, marker='s',
                   c='purple', hatch='||||')


    c.visualize(addcbar=False)

    c.set_geometry()
    assert_true(len(c._geometries) == 0)


@image_comparison(baseline_images=['test_text'], tol=8,
                  extensions=['png'])
def test_text():

    # UL Corner
    g = Grid(nxny=(5, 4), dxdy=(10, 10), ll_corner=(-20, -15), proj=wgs84,
             pixel_ref='corner')
    c = Map(g, ny=4, countries=False)
    c.set_lonlat_countours(interval=5., colors='crimson')

    c.set_text(-5, -5, 'Less Middle', color='green', style='italic', size=25)
    c.set_geometry(shpg.Point(-10, -10), s=500, marker='o',
                   text='My point', text_delta=[0, 0])

    shape = salem.utils.read_shapefile_to_grid(files['world_borders'], c.grid)
    had_c = set()
    for index, row in shape.iloc[::-1].iterrows():
        if row.CNTRY_NAME in had_c:
            c.set_geometry(row.geometry, crs=c.grid)
        else:
            c.set_geometry(row.geometry, text=row.CNTRY_NAME, crs=c.grid,
                           text_kwargs=dict(horizontalalignment='center',
                           verticalalignment='center', clip_on=True,
                                            color='gray'), text_delta=[0, 0])
        had_c.add(row.CNTRY_NAME)

    c.set_points([20, 20, 10], [10, 20, 20], s=250, marker='s',
                   c='purple', hatch='||||', text='baaaaad', text_delta=[0, 0],
                   text_kwargs=dict(horizontalalignment='center',
                                    verticalalignment='center', color='red'))

    c.visualize(addcbar=False)

    c.set_text()
    assert_true(len(c._text) == 0)

    c.set_geometry()
    assert_true(len(c._geometries) == 0)


@image_comparison(baseline_images=['test_hef_1','test_hef_2','test_hef_3',
                                   'test_hef_5', 'test_hef_array',
                                   'test_hef_topo'],
                  extensions=['png'])
def test_hef():

    grid = salem.grids.local_mercator_grid(center_ll=(10.76, 46.798444),
                                               extent=(10000, 7000))
    c = Map(grid, countries=False)
    c.set_lonlat_countours(interval=10)
    c.set_shapefile(get_demo_file('Hintereisferner_UTM.shp'))
    c.set_topography(get_demo_file('hef_srtm.tif'),
                     interp='linear')
    c.visualize(addcbar=False, title='linear')

    c.set_topography(get_demo_file('hef_srtm.tif'),
                     interp='spline', ks=2)
    c.visualize(addcbar=False, title='spline deg 2')

    c.set_topography(get_demo_file('hef_srtm.tif'))
    c.visualize(addcbar=False, title='Default: spline deg 3')

    h = c.set_topography(get_demo_file('hef_srtm.tif'),
                     interp='spline', ks=5)
    # This test is important, I removed the add_dcbar
    c.visualize(title='spline deg 5')

    dem = salem.GeoTiff(get_demo_file('hef_srtm.tif'))
    mytopo = dem.get_vardata()
    c.set_topography(mytopo, crs=dem.grid, interp='spline')
    c.visualize(addcbar=False, title='From array')

    c.set_lonlat_countours()
    c.set_cmap(cleo.get_cm('topo'))
    c.set_plot_params(nlevels=256)
    c.set_data(h)
    c.visualize()


@image_comparison(baseline_images=['test_googlemap', 'test_googlegrid'],
                  extensions=['png'])
def test_gmap():

    g = GoogleCenterMap(center_ll=(10.762660, 46.794221), zoom=13,
                        size_x=640, size_y=640)

    m = Map(g.grid, countries=False, nx=640)
    m.set_lonlat_countours(interval=0.025)
    m.set_shapefile(get_demo_file('Hintereisferner.shp'),
                    linewidths=2, edgecolor='darkred')
    m.set_rgb(g.get_vardata())
    m.visualize(addcbar=False)

    dem = salem.GeoTiff(get_demo_file('hef_srtm.tif'))
    dem.set_subset(margin=-100)

    dem = salem.grids.local_mercator_grid(center_ll=(10.76, 46.798444),
                                               extent=(10000, 7000))

    i, j = dem.ij_coordinates
    g = GoogleVisibleMap(x=i, y=j, src=dem, size_x=500, size_y=400)
    img = g.get_vardata()

    m = Map(dem, countries=False)

    assert_raises(ValueError, m.set_data, img)

    m.set_lonlat_countours(interval=0.025)
    m.set_shapefile(get_demo_file('Hintereisferner.shp'),
                    linewidths=2, edgecolor='darkred')
    m.set_rgb(img, g.grid)
    m.visualize(addcbar=False)


@image_comparison(baseline_images=['test_googlemap_llconts'],
                  extensions=['png'])
def test_gmap_llconts():

    # This was because some problems were left unnoticed by other tests
    g = salem.GoogleCenterMap(center_ll=(11.38, 47.26), zoom=9)
    m = cleo.Map(g.grid)
    m.set_rgb(g.get_vardata())
    m.set_lonlat_countours(interval=0.2)
    m.visualize(addcbar=False)


if __name__ == '__main__':
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        nose.runmodule(argv=['-s', '-v', '--with-doctest'], exit=False)
