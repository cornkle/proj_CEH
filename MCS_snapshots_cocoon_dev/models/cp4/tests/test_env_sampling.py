"""Synthetic regression tests; no real CP4/OBS files are read or written.

Run from models/cp4/:
    python -m unittest discover -s tests -p 'test_env_sampling.py' -v
"""
from pathlib import Path
import importlib.util
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import sample_cp4_environments as env


class TestSpatialAndOutput(unittest.TestCase):
    def test_spatial_box_is_mandatory_0p7_not_nearest_pixel(self):
        # Nonuniform values prove we're averaging a 3x3 support, not one pixel.
        lat = np.array([0.8, 0.6, 0.4, 0.2, 0.0])  # descending
        lon = np.array([-0.4, -0.2, 0.0, 0.2, 0.4])
        field = np.arange(25, dtype=float).reshape(5, 5)
        da = xr.DataArray(field, dims=('lat', 'lon'),
                          coords={'lat': lat, 'lon': lon})
        got = env.sample_box_means(da, lat, lon, [0.4], [0.0])
        expected = field[1:4, 1:4].mean()
        np.testing.assert_allclose(got, [expected])
        self.assertNotEqual(expected, float(field[2, 2]) + 1)  # sanity
        self.assertEqual(env.BOX_RADIUS_DEG, 0.35)

    def test_partial_domain_boxes_are_not_silently_shrunk(self):
        lat = np.linspace(0, 1, 11)
        lon = np.linspace(0, 1, 11)
        da = xr.DataArray(np.ones((11, 11)), dims=('lat', 'lon'))
        got = env.sample_box_means(da, lat, lon, [0.03], [0.5])
        self.assertTrue(np.isnan(got[0]))

    def test_african_grid_longitude_across_360(self):
        lat = np.arange(0, 1.01, 0.1)
        lon = np.arange(359.0, 361.01, 0.1)
        data = np.tile(lon, (len(lat), 1))
        da = xr.DataArray(data, dims=('lat', 'lon'))
        got = env.sample_box_means(da, lat, lon, [0.5], [0.0])
        self.assertAlmostEqual(got[0], 360.0, places=8)

    def test_boxes_are_disjoint_when_extended_to_full_domain(self):
        self.assertEqual(env.locate_regions(-5, 10), ['WAf'])
        self.assertFalse(env.locate_regions(-50, 30))
        root = Path('/tmp/cp4-test-env')
        inside = env.output_path(root, 'A-HIST', 'q850', 2089, 8, 18, scope='boxes')
        outside = env.output_path(root, 'A-HIST', 'q850', 2089, 8, 18, scope='outside')
        self.assertNotEqual(inside, outside)
        self.assertEqual(inside.name, 'part-18utc.parquet')
        self.assertEqual(outside.name, 'part-18utc-outside.parquet')
        self.assertEqual(inside.parent, outside.parent)

    def test_shear_components_magnitude(self):
        sampled = {
            'ua850': np.array([1., 2.]), 'va850': np.array([3., 4.]),
            'ua600': np.array([4., 6.]), 'va600': np.array([7., 7.]),
        }
        np.testing.assert_allclose(env.compute_output('ushear850_600', sampled), [3, 4])
        np.testing.assert_allclose(env.compute_output('vshear850_600', sampled), [4, 3])
        np.testing.assert_allclose(env.compute_output('shear850_600', sampled), [5, 5])

    def test_variable_output_scale_metadata(self):
        req = pd.DataFrame({'storm_id':['x'], 'rel_hour':[-6],
                            'env_time_key':['208908011200']})
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / 'v.parquet'
            # Can't assume pyarrow exists in all local test environments;
            # use parquet test only when installed (always available on JASMIN).
            if not importlib.util.find_spec('pyarrow'):
                self.skipTest('pyarrow unavailable')
            env.save_one('q850', np.array([0.01]), req, p)
            df = pd.read_parquet(p)
            self.assertEqual(df.spatial_sampling.iloc[0], 'box_mean_0p7deg')
            self.assertEqual(df.sampling_box_width_deg.iloc[0], 0.7)


@unittest.skipUnless(importlib.util.find_spec('cftime') is not None,
                     'cftime required for model calendars (available on JASMIN)')
class TestCalendarAndSampling(unittest.TestCase):
    def test_360_day_month_boundary(self):
        storms = pd.DataFrame([dict(
            run='Test-HIST', storm_id='x', year=2089, month=8, day=30,
            hour=17, minute=30, calendar='360_day', clon=5., clat=10.,
            tminlat=10., tminlon=5., storm_hour_end_utc=18.)])
        req = env.make_requests(storms, [18], [-6, 0, 6])
        self.assertEqual(req['era_utc_time'].tolist(), [
            '2089-08-30 12:00:00', '2089-08-30 18:00:00',
            '2089-09-01 00:00:00'])

    def test_native_pressure_field_sampled_at_0p7_deg(self):
        import cftime
        with tempfile.TemporaryDirectory() as tmp:
            run, var = 'Test-HIST', 'ua600'
            p = Path(tmp) / run / var
            p.mkdir(parents=True)
            f = p / 'ua600_GC6hrPtUV_x_208908010600-208909010000.nc'
            tt = [cftime.Datetime360Day(2089, 8, 1, 12),
                  cftime.Datetime360Day(2089, 8, 1, 18)]
            lat = [0., 0.2, 0.4, 0.6, 0.8]
            lon = [-0.4, -0.2, 0., 0.2, 0.4]
            vals = np.stack([np.arange(25).reshape(5, 5),
                             np.arange(25).reshape(5, 5) + 50]).astype('f4')
            ds = xr.Dataset({var: (('time', 'lat', 'lon'), vals)},
                            coords={'time': tt, 'lat': lat, 'lon': lon})
            ds.to_netcdf(f, engine='scipy')
            req = pd.DataFrame(dict(env_time_key=[env.native_key(t) for t in tt],
                                    calendar=['360_day']*2,
                                    tminlat=[0.4]*2, tminlon=[0.0]*2))
            out = env.sample_one_input(run, var, req, Path(tmp))
            np.testing.assert_allclose(out, [12, 62])

    def test_hourly_mean_at_half_hour_is_ending_environment(self):
        import cftime
        with tempfile.TemporaryDirectory() as tmp:
            run, var = 'Test-HIST', 'uas'
            p = Path(tmp) / run / var
            p.mkdir(parents=True)
            f = p / 'uas_GC1hr_x_208908010030-208908302330.nc'
            tt = [cftime.Datetime360Day(2089, 8, 1, 11, 30)]
            lat = [0., 0.2, 0.4, 0.6, 0.8]
            lon = [-0.4, -0.2, 0., 0.2, 0.4]
            ds = xr.Dataset({var: (('time', 'lat', 'lon'),
                                   np.ones((1, 5, 5), dtype='f4') * 5)},
                            coords={'time': tt, 'lat': lat, 'lon': lon})
            ds[var].attrs['cell_methods'] = 'area: time: mean'
            ds.to_netcdf(f, engine='scipy')
            req = pd.DataFrame(dict(env_time_key=['208908011200'],
                                    calendar=['360_day'],
                                    tminlat=[0.4], tminlon=[0.0]))
            out = env.sample_one_input(run, var, req, Path(tmp))
            np.testing.assert_allclose(out, [5.])


if __name__ == '__main__':
    unittest.main()
