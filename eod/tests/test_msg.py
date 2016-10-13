from eod import msg
import numpy as np
from numpy.testing import assert_array_equal, assert_allclose
import xarray as xr
import unittest
import datetime as dt
import pandas as pd

test_dir = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/eod/tests/test_files/msg'
test_write = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/eod/tests/test_write'
date = '200908010130'

class TestMsgRead(unittest.TestCase):

    def test_read_msg(self):

        obj = msg.ReadMsg(test_dir)

        assert obj.root == test_dir
        assert obj.years == ['2009']

    def test_set_date(self):
        obj = msg.ReadMsg(test_dir)
        obj.set_date(2009,8,1,1,30)

        assert obj.bpath == test_dir+'/cell_blob_files/2009/08/200908010130.gra'
        assert obj.tpath == test_dir+'/bigcell_area_table/rewrite/cell_40c_0130_JJAS.txt'

        obj = msg.ReadMsg(test_dir)
        obj.set_date(2009, 8, 1, 14, 30)

        assert obj.bpath == False
        assert obj.tpath == False

    def test_get_data(self):
        obj = msg.ReadMsg(test_dir)
        obj.set_date(2009, 8, 1, 1, 30)
        d = obj.get_data()

        assert d['t'].values.shape == (580,1640)
        assert d['time.month'].values == 8
        assert d['time.day'].values == 1
        assert d['time.minute'].values == 30

        with self.assertRaises(OSError):
            da = obj.get_data(netcdf_path='/does/not/exist/test.nc')

        da = obj.get_data(netcdf_path= test_write+'/test.nc')

    def test_box(self):
        obj = msg.ReadMsg(test_dir)
        obj.set_date(2009, 8, 1, 1, 30)
        d = obj.get_data(llbox=[-10, 8, 10, 12])

        assert np.isclose( d['lat'].values.max(), 12, atol=0.2)
        assert np.isclose(d['lat'].values.min(), 8, atol=0.2)
        assert np.isclose(d['lon'].values.max(), 10, atol=0.2)
        assert np.isclose(d['lon'].values.min(), -10, atol=0.2)

        obj.set_date(2009, 8, 1, 14, 30)
        d = obj.get_data(llbox=[-5, 3, 2, 7])

        assert np.isclose(d['lat'].values.max(), 7, atol=1)
        assert np.isclose(d['lat'].values.min(), 3, atol=1)
        assert np.isclose(d['lon'].values.max(), 2, atol=1)
        assert np.isclose(d['lon'].values.min(), -5, atol=1)





