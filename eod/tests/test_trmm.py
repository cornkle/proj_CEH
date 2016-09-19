from eod import trmm
import numpy as np
from numpy.testing import assert_array_equal, assert_allclose
import xarray as xr
import unittest

test_dir = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/eod/tests/test_files/trmm'
test_write = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/eod/tests/test_write'

class TestTrmmRead(unittest.TestCase):

    def test_trmmReadWA(self):

        obj = trmm.ReadWA(test_dir)
        files = obj.fpaths

        assert files == [
            test_dir+'/2007/08/2A25.20070816.55562.7.gra',
            test_dir+'/2007/08/2A25.20070817.55577.7.gra']

    def test_area(self):

        area = [1, 4, 1, 4]
        lat = np.array([[4, 3, 2, 1], ] * 4).transpose()
        lon = np.array([[1, 2, 3, 4], ] * 4)
        lat_result = [3, 3, 2, 2]
        lon_result = [2, 3, 2, 3]

        box = np.where((lon > area[0]) & (lon < area[1]) & (lat > area[2]) & (lat < area[3]))

        assert isinstance(box, type(([], [])))   # np.where is tuple!
        assert len(box[0]) == 4
        assert_array_equal( lat[box], lat_result)
        assert_array_equal(lon[box], lon_result)

    def test_getData_area_cut_date(self):

        # both testfiles have enough rainfall in swath at 0.5 mm threshold and given box
        area=[-15, 5, 15, 20]
        min_rain_swath = 2000
        min_rain_box = 500
        min_tpixel = 2500
        rain_thresh = 0.5

        obj = trmm.ReadWA(test_dir, area=area)

        td=obj.get_data(obj.fpaths[1])
        box = np.where((td['lon'].values > area[0]) & (td['lon'].values < area[2]) & (td['lat'].values > area[1]) & (td['lat'].values < area[3]))
        # files are properly filtered according to the rain/box overlap thresholds
        assert len(box[0]) > min_tpixel
        assert np.sum(td.values[box] > rain_thresh) > min_rain_box

        # use cut on two files
        obj = trmm.ReadWA(test_dir)
        td = obj.get_data(obj.fpaths[0], cut = [8,10])
        assert td['lat'].values[:, 0].max() <= 10
        assert td['lat'].values[:, 0].min() >= 8

        td = obj.get_data(obj.fpaths[1])
        assert td['lat'].values[:, 0].max() >= 10
        assert td['lat'].values[:, 0].min() <= 8

        # get time and via index get same array, cut works the same
        td = obj.get_ddata(2007,8,16,19,16, cut = [8,10])
        assert td['lat'].values[:, 0].max() <= 10
        assert td['lat'].values[:, 0].min() >= 8

        td2 = obj.get_data(obj.fpaths[0], cut = [8, 10])
        assert_array_equal(td2.values, td.values)

    def test_write_netcdf(self):
        obj = trmm.ReadWA(test_dir)

        with self.assertRaises(OSError):
            da = obj.get_data(obj.fpaths[1], netcdf_path='/does/not/exist/test.nc')

        da = obj.get_data(obj.fpaths[1], netcdf_path= test_write+'/test.nc')




