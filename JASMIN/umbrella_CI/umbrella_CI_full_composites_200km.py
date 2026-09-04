#!/usr/bin/env python3
"""Build full rotated CI composites from DYAMOND3 5km-RAL3p3-tuned PP files.

The script:
  * reads aligned/right/reverse/left CI lists;
  * retains 14-UTC CIs by default and processes t=-2 and t=-1 separately;
  * excludes water within 200 km and terrain with P95-P5 relief > 700 m;
  * rotates each case so the 600-hPa box-mean wind points in the class direction;
  * calculates SM anomaly, near-surface q, OLR, preceding-hour rain, shear
    magnitude and signed rotated shear components, 10-m divergence, and
    rotated 10-m/600-hPa winds;
  * accumulates valid-pixel sums/counts and writes NetCDF composites and figures;
  * checkpoints by target-time group and resumes automatically.

Run `python umbrella_CI_full_composites.py --help` for options.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import logging
import os
import sys
from functools import lru_cache
from pathlib import Path

import iris
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from scipy.interpolate import RegularGridInterpolator


DIRECTION_ORDER = ["aligned", "right", "reverse", "left"]
TARGET_ANGLE = {"aligned": 90.0, "right": 0.0, "reverse": -90.0, "left": 180.0}
OUTPUT_SCHEMA_VERSION = 2
VARIABLES = {
    "q_low": ("apvera.pp", "m01s03i237"),
    "u10": ("apvera.pp", "m01s03i225"),
    "v10": ("apvera.pp", "m01s03i226"),
    "olr": ("apvera.pp", "m01s02i205"),
    "sm": ("apverb.pp", "m01s08i223"),
    "rain": ("apverb.pp", "m01s04i203"),
    "u600": ("apverc.pp", "m01s15i201"),
    "v600": ("apverc.pp", "m01s15i202"),
}
CUBE_SELECTOR = {"rain": "time_mean"}
TIME_METHOD = {
    "q_low": "exact", "u10": "exact", "v10": "exact", "sm": "exact",
    "rain": "exact", "olr": "nearest", "u600": "nearest", "v600": "nearest",
}
COMPOSITE_FIELDS = [
    "sm_anomaly", "q_low", "olr", "rain", "shear", "shear_u", "shear_v", "div10",
    "u10", "v10", "u600", "v600",
]
FAILURE_COLUMNS = ["category", "time", "sample_time", "lon", "lat", "error"]


def parse_args() -> argparse.Namespace:
    base = Path("/gws/ssde/j25b/kscale/DYAMOND3_reruns")
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--ci-dir", type=Path, default=Path("/home/users/cornkle/pythonWorkspace/proj_CEH/JASMIN/umbrella_CI"))
    parser.add_argument("--model-root", type=Path, default=base / "5km-RAL3p3-tuned/glm/field.pp")
    parser.add_argument("--orography", type=Path, default=base / "orography/orography-n2560e.nc")
    parser.add_argument("--landmask", type=Path, default=base / "landmask/mask-n2560e.nc")
    parser.add_argument("--output-dir", type=Path, default=Path("/gws/ssde/j25b/kscale/USERS/cklein/umbrella/full_composites"))
    parser.add_argument("--year", type=int, default=2020)
    parser.add_argument("--ci-hours", type=int, nargs="+", default=[14])
    parser.add_argument("--lag-hours", type=int, nargs="+", default=[-2, -1],
                        help="One or more sampling times relative to CI")
    parser.add_argument("--midlevel-pressure", type=float, default=600.0, help="hPa")
    parser.add_argument("--water-radius-km", type=float, default=200.0)
    parser.add_argument("--topo-radius-km", type=float, default=200.0)
    parser.add_argument("--topo-relief-max", type=float, default=700.0, help="Maximum P95-P5 orographic relief (m)")
    parser.add_argument("--half-width-km", type=float, default=250.0)
    parser.add_argument("--grid-spacing-km", type=float, default=5.0)
    parser.add_argument("--wind-box-half-width-km", type=float, default=100.0)
    parser.add_argument("--drop-exact-duplicates", action="store_true", help="Remove repeated lon/lat/date/hour rows")
    parser.add_argument("--max-cases-per-category", type=int, default=None, help="Optional reproducible test limit; default processes all")
    parser.add_argument("--random-seed", type=int, default=42)
    parser.add_argument("--checkpoint-every", type=int, default=5, help="Checkpoint after this many target-time groups")
    parser.add_argument("--overwrite", action="store_true", help="Discard an existing checkpoint for this lag")
    parser.add_argument("--lazy-subsets", action="store_true", help="Do not preload full 2-D fields per target time; uses less RAM but repeats PP unpacking")
    parser.add_argument("--no-plots", action="store_true")
    return parser.parse_args()


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stdout,
    )


def lag_tag(lag: int) -> str:
    return f"lag_{lag:+d}h".replace("+", "p").replace("-", "m")


def read_ci_files(ci_dir: Path, year: int, drop_duplicates: bool) -> pd.DataFrame:
    frames = []
    for category in DIRECTION_ORDER:
        path = ci_dir / f"{category}_CIs.txt"
        if not path.exists():
            raise FileNotFoundError(path)
        frame = pd.read_csv(path, sep=r"\s+", header=None, names=["lon", "lat", "month", "day", "hour"])
        frame["category"] = category
        frame["time"] = pd.to_datetime(dict(
            year=np.full(len(frame), year), month=frame["month"], day=frame["day"], hour=frame["hour"]
        ))
        frame["exact_duplicate"] = frame.duplicated(["lon", "lat", "month", "day", "hour"], keep=False)
        if drop_duplicates:
            frame = frame.drop_duplicates(["lon", "lat", "month", "day", "hour"])
        logging.info("%s: %d rows (%d exact duplicate rows)", category, len(frame), int(frame["exact_duplicate"].sum()))
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def wrapped_lon_difference(lon: np.ndarray, lon0: float) -> np.ndarray:
    return (lon - lon0 + 180.0) % 360.0 - 180.0


def calculate_location_metrics(
    cases: pd.DataFrame, topo_file: Path, mask_file: Path,
    water_radius_km: float, topo_radius_km: float,
) -> pd.DataFrame:
    with xr.open_dataset(topo_file) as topo_ds, xr.open_dataset(mask_file) as mask_ds:
        topo = topo_ds["surface_altitude"].load().values.astype(np.float32)
        land = mask_ds["land_binary_mask"].load().values >= 0.5
        lat = topo_ds["latitude"].values.astype(float)
        lon = topo_ds["longitude"].values.astype(float)
        if topo.shape != land.shape or not np.allclose(lat, mask_ds["latitude"].values) or not np.allclose(lon, mask_ds["longitude"].values):
            raise ValueError("Orography and land-mask grids do not match")

    search_radius = max(water_radius_km, topo_radius_km)
    locations = cases[["lon", "lat"]].drop_duplicates().reset_index(drop=True)
    records = []
    for number, row in enumerate(locations.itertuples(index=False), 1):
        lon0, lat0 = float(row.lon), float(row.lat)
        dlat_max = search_radius / 110.574
        dlon_max = search_radius / (111.320 * np.cos(np.deg2rad(lat0)))
        iy = np.where(np.abs(lat - lat0) <= dlat_max)[0]
        dlon = wrapped_lon_difference(lon, lon0)
        ix = np.where(np.abs(dlon) <= dlon_max)[0]
        if not len(iy) or not len(ix):
            records.append({"lon": lon0, "lat": lat0, "central_land": False, "water_distance_km": np.nan,
                            "topo_mean_200km": np.nan, "topo_std_200km": np.nan,
                            "topo_relief_200km": np.nan, "topo_max_200km": np.nan})
            continue

        dy = (lat[iy] - lat0) * 110.574
        dx = dlon[ix] * 111.320 * np.cos(np.deg2rad(lat0))
        distance = np.hypot(dy[:, None], dx[None, :])
        local_topo = topo[np.ix_(iy, ix)]
        local_land = land[np.ix_(iy, ix)]
        cy = int(np.argmin(np.abs(lat - lat0)))
        cx = int(np.argmin(np.abs(wrapped_lon_difference(lon, lon0))))

        water = (distance <= water_radius_km) & ~local_land
        water_distance = float(np.min(distance[water])) if np.any(water) else np.inf
        valid_topo = (distance <= topo_radius_km) & local_land & np.isfinite(local_topo)
        elevation = local_topo[valid_topo]
        if elevation.size:
            p05, p95 = np.percentile(elevation, [5, 95])
            topo_values = (float(np.mean(elevation)), float(np.std(elevation)), float(p95 - p05), float(np.max(elevation)))
        else:
            topo_values = (np.nan, np.nan, np.nan, np.nan)
        records.append({
            "lon": lon0, "lat": lat0, "central_land": bool(land[cy, cx]), "water_distance_km": water_distance,
            "topo_mean_200km": topo_values[0], "topo_std_200km": topo_values[1],
            "topo_relief_200km": topo_values[2], "topo_max_200km": topo_values[3],
        })
        if number % 500 == 0 or number == len(locations):
            logging.info("Ancillary metrics: %d/%d unique locations", number, len(locations))
    return pd.DataFrame(records)


class Field2D:
    """Realized float32 field and its own latitude/longitude coordinates."""
    __slots__ = ("lats", "lons", "data")

    def __init__(self, lats, lons, data):
        self.lats = lats
        self.lons = lons
        self.data = data


class CompositeProcessor:
    def __init__(self, args: argparse.Namespace):
        self.args = args
        self.model_root = args.model_root
        self.run_name = "glm.n2560_RAL3p3_tuned"
        self.xy = np.arange(-args.half_width_km, args.half_width_km + args.grid_spacing_km, args.grid_spacing_km)
        self.X, self.Y = np.meshgrid(self.xy, self.xy)

    def pp_path(self, stream: str, cycle: pd.Timestamp) -> Path:
        name = stream.removesuffix(".pp")
        stamp = pd.Timestamp(cycle).strftime("%Y%m%dT%H%MZ")
        return self.model_root / stream / f"{self.run_name}.{name}_{stamp}.pp"

    @staticmethod
    def take_coord_index(cube, coord_name: str, index: int):
        dim = cube.coord_dims(coord_name)[0]
        selection = [slice(None)] * cube.ndim
        selection[dim] = index
        return cube[tuple(selection)]

    def select_vertical_level(self, cube, variable: str):
        if variable == "sm":
            return self.take_coord_index(cube, "depth", int(np.argmin(cube.coord("depth").points)))
        if variable in ("u600", "v600"):
            pressure = cube.coord("pressure")
            pressure.convert_units("hPa")
            index = int(np.argmin(np.abs(pressure.points - self.args.midlevel_pressure)))
            if not np.isclose(pressure.points[index], self.args.midlevel_pressure):
                raise ValueError(f"{self.args.midlevel_pressure} hPa is unavailable")
            return self.take_coord_index(cube, "pressure", index)
        return cube

    @staticmethod
    def has_time_mean(cube) -> bool:
        return any(method.method == "mean" and "time" in method.coord_names for method in cube.cell_methods)

    @lru_cache(maxsize=8)
    def _load_file_cubes(self, stream: str, cycle_string: str):
        path = self.pp_path(stream, pd.Timestamp(cycle_string))
        return iris.load(str(path))

    def load_base_cube(self, stream: str, stash: str, cycle_string: str, selector: str = "single"):
        cubes = [cube for cube in self._load_file_cubes(stream, cycle_string)
                 if str(cube.attributes.get("STASH", "")) == stash]
        if selector == "time_mean":
            cubes = [cube for cube in cubes if self.has_time_mean(cube)]
        if len(cubes) != 1:
            details = [(cube.shape, cube.cell_methods) for cube in cubes]
            raise ValueError(f"Expected one cube for {stash}, selector={selector}; found {len(cubes)}: {details}")
        return cubes[0]

    @staticmethod
    def as_timestamp(date) -> pd.Timestamp:
        return pd.Timestamp(date.year, date.month, date.day, date.hour, date.minute, date.second)

    def effective_times(self, cube, variable: str) -> list[pd.Timestamp]:
        time = cube.coord("time")
        values = time.bounds[:, 1] if variable == "rain" else time.points
        return [self.as_timestamp(date) for date in time.units.num2date(values)]

    def available_time_slices(self, variable: str, target_time: pd.Timestamp) -> dict[pd.Timestamp, object]:
        stream, stash = VARIABLES[variable]
        selector = CUBE_SELECTOR.get(variable, "single")
        target_time = pd.Timestamp(target_time)
        current_cycle = target_time.floor("12h")
        available = {}
        for cycle in (current_cycle, current_cycle - pd.Timedelta(hours=12)):
            path = self.pp_path(stream, cycle)
            if not path.exists():
                continue
            cube = self.load_base_cube(stream, stash, str(cycle), selector)
            cube = self.select_vertical_level(cube, variable)
            for index, valid_time in enumerate(self.effective_times(cube, variable)):
                available.setdefault(valid_time, self.take_coord_index(cube, "time", index))
        if not available:
            raise FileNotFoundError(f"No {variable} data found around {target_time}")
        return available

    def load_field_at_time(self, variable: str, target_time: pd.Timestamp):
        target_time = pd.Timestamp(target_time)
        available = self.available_time_slices(variable, target_time)
        valid_times = sorted(available)
        method = TIME_METHOD[variable]
        if method == "exact":
            if target_time not in available:
                raise ValueError(f"No exact {variable} field at {target_time}; available={valid_times}")
            used_time = target_time
        else:
            used_time = min(valid_times, key=lambda time: abs(time - target_time))
        return available[used_time], used_time

    @staticmethod
    def consecutive_runs(indices: np.ndarray) -> list[np.ndarray]:
        return np.split(indices, np.where(np.diff(indices) > 1)[0] + 1)

    def extract_local_data(self, cube, lon0: float, lat0: float, half_width_km: float):
        lats = cube.coord("latitude").points
        lons = cube.coord("longitude").points
        lat_pad = half_width_km / 110.574 + 2 * np.median(np.diff(lats))
        lon_pad = half_width_km / (111.320 * np.cos(np.deg2rad(lat0))) + 2 * np.median(np.diff(lons))
        lat_indices = np.where(np.abs(lats - lat0) <= lat_pad)[0]
        lon_difference = wrapped_lon_difference(lons, lon0)
        lon_indices = np.where(np.abs(lon_difference) <= lon_pad)[0]
        if not len(lat_indices) or not len(lon_indices):
            raise ValueError(f"No grid cells found around lon={lon0}, lat={lat0}")
        lat_slice = slice(lat_indices[0], lat_indices[-1] + 1)
        arrays, relative_lons = [], []
        for run in self.consecutive_runs(lon_indices):
            selection = [slice(None)] * cube.ndim
            selection[cube.coord_dims("latitude")[0]] = lat_slice
            selection[cube.coord_dims("longitude")[0]] = slice(run[0], run[-1] + 1)
            subset = cube[tuple(selection)]
            data = np.ma.filled(np.ma.asarray(subset.data, dtype=float), np.nan)
            data = np.moveaxis(data, (subset.coord_dims("latitude")[0], subset.coord_dims("longitude")[0]), (0, 1))
            arrays.append(data)
            relative_lons.append(lon_difference[run])
        data = np.concatenate(arrays, axis=1)
        relative_lons = np.concatenate(relative_lons)
        order = np.argsort(relative_lons)
        local_x = relative_lons[order] * 111.320 * np.cos(np.deg2rad(lat0))
        local_y = (lats[lat_slice] - lat0) * 110.574
        return data[:, order], local_x, local_y

    @staticmethod
    def extract_local_array(field, lon0: float, lat0: float, half_width_km: float):
        lats, lons, data = field.lats, field.lons, field.data
        lat_pad = half_width_km / 110.574 + 2 * np.median(np.diff(lats))
        lon_pad = half_width_km / (111.320 * np.cos(np.deg2rad(lat0))) + 2 * np.median(np.diff(lons))
        lat_indices = np.where(np.abs(lats - lat0) <= lat_pad)[0]
        lon_difference = wrapped_lon_difference(lons, lon0)
        lon_indices = np.where(np.abs(lon_difference) <= lon_pad)[0]
        if not len(lat_indices) or not len(lon_indices):
            raise ValueError(f"No grid cells found around lon={lon0}, lat={lat0}")
        lat_slice = slice(lat_indices[0], lat_indices[-1] + 1)
        arrays, relative_lons = [], []
        for run in CompositeProcessor.consecutive_runs(lon_indices):
            arrays.append(data[lat_slice, run[0]:run[-1] + 1])
            relative_lons.append(lon_difference[run])
        block = np.concatenate(arrays, axis=1)
        relative_lons = np.concatenate(relative_lons)
        order = np.argsort(relative_lons)
        local_x = relative_lons[order] * 111.320 * np.cos(np.deg2rad(lat0))
        local_y = (lats[lat_slice] - lat0) * 110.574
        return block[:, order], local_x, local_y

    def realize_field(self, cube, lat_min: float, lat_max: float):
        lats = cube.coord("latitude").points.astype(np.float64)
        lons = cube.coord("longitude").points.astype(np.float64)
        pad = (np.sqrt(2) * self.args.half_width_km + 2 * self.args.grid_spacing_km) / 110.574 + 0.5
        keep = np.where((lats >= lat_min - pad) & (lats <= lat_max + pad))[0]
        if not len(keep):
            raise ValueError(f"No latitude rows for group range {lat_min} to {lat_max}")
        selection = [slice(None)] * cube.ndim
        selection[cube.coord_dims("latitude")[0]] = slice(keep[0], keep[-1] + 1)
        subset = cube[tuple(selection)]
        data = np.ma.filled(np.ma.asarray(subset.data), np.nan).astype(np.float32, copy=False)
        data = np.moveaxis(data, (subset.coord_dims("latitude")[0], subset.coord_dims("longitude")[0]), (0, 1))
        return Field2D(lats[keep[0]:keep[-1] + 1], lons, data)

    def interpolate_rotated(self, data, source_x, source_y, theta_deg):
        theta = np.deg2rad(theta_deg)
        c, s = np.cos(theta), np.sin(theta)
        earth_x = c * self.X + s * self.Y
        earth_y = -s * self.X + c * self.Y
        interpolator = RegularGridInterpolator((source_y, source_x), data, bounds_error=False, fill_value=np.nan)
        points = np.column_stack([earth_y.ravel(), earth_x.ravel()])
        return interpolator(points).reshape(self.X.shape)

    @staticmethod
    def central_box_mean(data, source_x, source_y, half_width):
        return np.nanmean(data[np.ix_(np.abs(source_y) <= half_width, np.abs(source_x) <= half_width)])

    def process_one_case(self, case: pd.Series, fields: dict):
        lon0, lat0, category = float(case["lon"]), float(case["lat"]), case["category"]
        source_half_width = np.sqrt(2) * self.args.half_width_km + 2 * self.args.grid_spacing_km
        local, coordinates = {}, {}
        for variable, field in fields.items():
            if isinstance(field, Field2D):
                data, source_x, source_y = self.extract_local_array(field, lon0, lat0, source_half_width)
            else:
                data, source_x, source_y = self.extract_local_data(field, lon0, lat0, source_half_width)
            local[variable], coordinates[variable] = data, (source_x, source_y)

        u600_mean = self.central_box_mean(local["u600"], *coordinates["u600"], self.args.wind_box_half_width_km)
        v600_mean = self.central_box_mean(local["v600"], *coordinates["v600"], self.args.wind_box_half_width_km)
        if not np.isfinite(u600_mean) or not np.isfinite(v600_mean) or np.hypot(u600_mean, v600_mean) < 0.1:
            raise ValueError("Invalid or near-zero 600-hPa box-mean wind")
        actual_angle = np.rad2deg(np.arctan2(v600_mean, u600_mean))
        rotation_angle = TARGET_ANGLE[category] - actual_angle

        rotated_scalar = {
            name: self.interpolate_rotated(local[name], *coordinates[name], rotation_angle)
            for name in ("sm", "q_low", "olr", "rain")
        }
        earth_wind = {
            name: self.interpolate_rotated(local[name], *coordinates[name], rotation_angle)
            for name in ("u10", "v10", "u600", "v600")
        }
        theta = np.deg2rad(rotation_angle)
        c, s = np.cos(theta), np.sin(theta)
        u10 = c * earth_wind["u10"] - s * earth_wind["v10"]
        v10 = s * earth_wind["u10"] + c * earth_wind["v10"]
        u600 = c * earth_wind["u600"] - s * earth_wind["v600"]
        v600 = s * earth_wind["u600"] + c * earth_wind["v600"]
        shear_u = u600 - u10
        shear_v = v600 - v10
        spacing_m = self.args.grid_spacing_km * 1000
        divergence = np.gradient(u10, spacing_m, axis=1) + np.gradient(v10, spacing_m, axis=0)
        return {
            "sm_anomaly": rotated_scalar["sm"] - np.nanmean(rotated_scalar["sm"]),
            "q_low": rotated_scalar["q_low"] * 1000,
            "olr": rotated_scalar["olr"],
            "rain": rotated_scalar["rain"] * 3600,
            "shear": np.hypot(shear_u, shear_v),
            "shear_u": shear_u,
            "shear_v": shear_v,
            "div10": divergence,
            "u10": u10, "v10": v10, "u600": u600, "v600": v600,
        }

    def load_fields_for_time(self, target_time: pd.Timestamp, lat_min=-90.0, lat_max=90.0) -> dict:
        fields = {}
        used = {}
        for variable in VARIABLES:
            fields[variable], used[variable] = self.load_field_at_time(variable, target_time)
        logging.info("Field times for %s: %s", target_time, ", ".join(f"{k}={v:%H:%M}" for k, v in used.items()))
        if self.args.lazy_subsets:
            return fields
        logging.info("Realizing latitude-cropped float32 fields for %s", target_time)
        return {variable: self.realize_field(cube, lat_min, lat_max)
                for variable, cube in fields.items()}


def empty_accumulators(shape: tuple[int, int]):
    sums = {cat: {field: np.zeros(shape, dtype=np.float64) for field in COMPOSITE_FIELDS} for cat in DIRECTION_ORDER}
    counts = {cat: {field: np.zeros(shape, dtype=np.int32) for field in COMPOSITE_FIELDS} for cat in DIRECTION_ORDER}
    successful = {cat: 0 for cat in DIRECTION_ORDER}
    attempted = {cat: 0 for cat in DIRECTION_ORDER}
    return sums, counts, successful, attempted


def case_signature(cases: pd.DataFrame, args: argparse.Namespace) -> str:
    columns = ["category", "lon", "lat", "time", "sample_time"]
    config = {
        "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "lag": args.lag_hours, "pressure": args.midlevel_pressure,
        "water_radius": args.water_radius_km, "topo_radius": args.topo_radius_km,
        "relief_max": args.topo_relief_max, "half_width": args.half_width_km,
        "spacing": args.grid_spacing_km, "wind_box": args.wind_box_half_width_km,
    }
    payload = cases[columns].to_csv(index=False).encode() + json.dumps(config, sort_keys=True).encode()
    return hashlib.sha256(payload).hexdigest()


def save_checkpoint(path, signature, next_group, sums, counts, successful, attempted):
    payload = {"signature": np.array(signature), "next_group": np.array(next_group, dtype=np.int32)}
    for category in DIRECTION_ORDER:
        payload[f"successful__{category}"] = np.array(successful[category], dtype=np.int64)
        payload[f"attempted__{category}"] = np.array(attempted[category], dtype=np.int64)
        for field in COMPOSITE_FIELDS:
            payload[f"sum__{category}__{field}"] = sums[category][field]
            payload[f"count__{category}__{field}"] = counts[category][field]
    temporary = path.with_name(path.name + ".tmp.npz")
    np.savez_compressed(temporary, **payload)
    os.replace(temporary, path)


def load_checkpoint(path, signature, shape):
    sums, counts, successful, attempted = empty_accumulators(shape)
    if not path.exists():
        return 0, sums, counts, successful, attempted
    with np.load(path, allow_pickle=False) as data:
        saved_signature = data["signature"].item()
        if saved_signature != signature:
            raise RuntimeError(f"Checkpoint signature mismatch: {path}. Use --overwrite or a different output directory.")
        next_group = int(data["next_group"])
        for category in DIRECTION_ORDER:
            successful[category] = int(data[f"successful__{category}"])
            attempted[category] = int(data[f"attempted__{category}"])
            for field in COMPOSITE_FIELDS:
                sums[category][field][...] = data[f"sum__{category}__{field}"]
                counts[category][field][...] = data[f"count__{category}__{field}"]
    return next_group, sums, counts, successful, attempted


def means_from_accumulators(sums, counts):
    return {
        category: {
            field: np.divide(sums[category][field], counts[category][field],
                             out=np.full_like(sums[category][field], np.nan),
                             where=counts[category][field] > 0)
            for field in COMPOSITE_FIELDS
        }
        for category in DIRECTION_ORDER
    }


def write_netcdf(path: Path, processor: CompositeProcessor, means, counts, successful, args):
    coords = {"category": DIRECTION_ORDER, "y_km": processor.xy, "x_km": processor.xy}
    variables = {}
    for field in COMPOSITE_FIELDS:
        variables[field] = (("category", "y_km", "x_km"), np.stack([means[c][field] for c in DIRECTION_ORDER]))
        variables[f"{field}_count"] = (("category", "y_km", "x_km"), np.stack([counts[c][field] for c in DIRECTION_ORDER]))
    variables["n_cases"] = (("category",), np.array([successful[c] for c in DIRECTION_ORDER], dtype=np.int32))
    ds = xr.Dataset(variables, coords=coords, attrs={
        "description": "Rotated convection-initiation composites",
        "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "lag_hours": args.lag_hours,
        "ci_hours": ",".join(map(str, args.ci_hours)),
        "midlevel_pressure_hpa": args.midlevel_pressure,
        "water_exclusion_radius_km": args.water_radius_km,
        "topographic_radius_km": args.topo_radius_km,
        "topographic_relief_max_m": args.topo_relief_max,
        "sm_anomaly_definition": "top-layer SM minus each rotated case's full plotted-box mean",
        "rain_definition": "one-hour mean stratiform rainfall flux ending at target time, converted to mm h-1",
        "divergence_units": "s-1",
    })
    units = {"sm_anomaly": "kg m-2", "q_low": "g kg-1", "olr": "W m-2", "rain": "mm h-1",
             "shear": "m s-1", "shear_u": "m s-1", "shear_v": "m s-1", "div10": "s-1",
             "u10": "m s-1", "v10": "m s-1", "u600": "m s-1", "v600": "m s-1"}
    for field, unit in units.items():
        ds[field].attrs["units"] = unit
    ds["shear"].attrs["long_name"] = "600 hPa minus 10 m wind shear magnitude"
    ds["shear_u"].attrs.update({
        "long_name": "rotated x-component of 600 hPa minus 10 m wind shear",
        "positive_direction": "right on the composite plot",
    })
    ds["shear_v"].attrs.update({
        "long_name": "rotated y-component of 600 hPa minus 10 m wind shear",
        "positive_direction": "up on the composite plot",
    })
    temporary = path.with_name(path.stem + ".tmp.nc")
    ds.to_netcdf(temporary)
    os.replace(temporary, path)


def plot_composite(path, processor, means, successful, variable, title, cmap, label,
                   symmetric=False, factor=1.0, floor_zero=False):
    plotted = {cat: means[cat][variable] * factor for cat in DIRECTION_ORDER}
    finite_parts = [data[np.isfinite(data)] for data in plotted.values() if np.any(np.isfinite(data))]
    if not finite_parts:
        logging.warning("Skipping empty plot: %s", variable)
        return
    values = np.concatenate(finite_parts)
    if symmetric:
        limit = np.percentile(np.abs(values), 98)
        kwargs = {"norm": TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)}
    else:
        vmin, vmax = np.percentile(values, [2, 98])
        kwargs = {"vmin": max(0, vmin) if floor_zero else vmin, "vmax": vmax}

    fig, axes = plt.subplots(1, 4, figsize=(17, 4.4), sharex=True, sharey=True, constrained_layout=True)
    step = 10
    for ax, category in zip(axes, DIRECTION_ORDER):
        field = means[category]
        shading = ax.pcolormesh(processor.X, processor.Y, plotted[category], cmap=cmap, shading="auto", **kwargs)
        ax.quiver(processor.X[::step, ::step], processor.Y[::step, ::step],
                  field["u600"][::step, ::step], field["v600"][::step, ::step],
                  color="black", scale=180, width=0.003)
        ax.quiver(processor.X[::step, ::step], processor.Y[::step, ::step],
                  field["u10"][::step, ::step], field["v10"][::step, ::step],
                  color="red", scale=180, width=0.003)
        ax.plot(0, 0, "+", color="magenta", ms=11, mew=2)
        ax.axhline(0, color="0.4", lw=0.5)
        ax.axvline(0, color="0.4", lw=0.5)
        ax.set(title=f"{category.capitalize()} (n={successful[category]})", xlabel="x (km)",
               xlim=(-processor.args.half_width_km, processor.args.half_width_km),
               ylim=(-processor.args.half_width_km, processor.args.half_width_km), aspect="equal")
    axes[0].set_ylabel("y (km)")
    fig.suptitle(title)
    fig.legend(handles=[Line2D([0], [0], color="black", lw=2, label="600 hPa wind"),
                        Line2D([0], [0], color="red", lw=2, label="10 m wind")],
               loc="upper center", ncol=2)
    fig.colorbar(shading, ax=axes.tolist(), label=label, shrink=0.82)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def run_one_lag(args: argparse.Namespace, cases_with_filters: pd.DataFrame) -> None:
    run_dir = args.output_dir / lag_tag(args.lag_hours)
    run_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_path = run_dir / "checkpoint.npz"
    failures_path = run_dir / "failures.csv"
    if args.overwrite:
        checkpoint_path.unlink(missing_ok=True)
        failures_path.unlink(missing_ok=True)

    config = {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()}
    (run_dir / "run_config.json").write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")

    cases = cases_with_filters.copy()
    cases["sample_time"] = cases["time"] + pd.to_timedelta(args.lag_hours, unit="h")
    filtered = cases[cases["sample_ok"]].copy()
    filtered = filtered.sort_values(["sample_time", "category", "time", "lat", "lon"], kind="stable").reset_index(drop=True)

    if args.max_cases_per_category is not None:
        filtered = pd.concat([
            group.sample(n=min(args.max_cases_per_category, len(group)), random_state=args.random_seed + i)
            for i, (_, group) in enumerate(filtered.groupby("category", sort=False))
        ], ignore_index=True).sort_values(["sample_time", "category", "time", "lat", "lon"], kind="stable").reset_index(drop=True)
    if filtered.empty:
        raise RuntimeError("No cases remain after filtering")

    cases.to_csv(run_dir / "all_selected_hour_cases_with_filters.csv", index=False)
    filtered.to_csv(run_dir / "filtered_cases.csv", index=False)
    filter_summary = pd.DataFrame({
        "all": cases.groupby("category").size(),
        "inland": cases[cases["coast_ok"]].groupby("category").size(),
        "final": filtered.groupby("category").size(),
    }).fillna(0).astype(int)
    filter_summary.to_csv(run_dir / "filter_summary.csv")
    logging.info("Filter summary:\n%s", filter_summary)

    processor = CompositeProcessor(args)
    signature = case_signature(filtered, args)
    time_groups = list(filtered.groupby("sample_time", sort=True))
    start_group, sums, counts, successful, attempted = load_checkpoint(checkpoint_path, signature, processor.X.shape)
    failures = []
    if failures_path.exists() and start_group and failures_path.stat().st_size:
        try:
            failures = pd.read_csv(failures_path).to_dict("records")
        except pd.errors.EmptyDataError:
            logging.warning("Ignoring empty failure log left by the previous run: %s", failures_path)
    logging.info("Processing %d cases in %d target-time groups; resuming at group %d", len(filtered), len(time_groups), start_group)

    for group_index in range(start_group, len(time_groups)):
        target_time, group = time_groups[group_index]
        logging.info("Time group %d/%d: %s (%d cases)", group_index + 1, len(time_groups), target_time, len(group))
        try:
            fields = processor.load_fields_for_time(target_time, float(group["lat"].min()), float(group["lat"].max()))
        except Exception as exc:
            logging.exception("Unable to load fields for %s", target_time)
            for _, case in group.iterrows():
                category = case["category"]
                attempted[category] += 1
                failures.append({"category": category, "time": case["time"], "sample_time": target_time,
                                 "lon": case["lon"], "lat": case["lat"], "error": repr(exc)})
            fields = None

        if fields is not None:
            for _, case in group.iterrows():
                category = case["category"]
                attempted[category] += 1
                try:
                    rotated = processor.process_one_case(case, fields)
                    for field in COMPOSITE_FIELDS:
                        data = rotated[field]
                        valid = np.isfinite(data)
                        sums[category][field][valid] += data[valid]
                        counts[category][field][valid] += 1
                    successful[category] += 1
                except Exception as exc:
                    failures.append({"category": category, "time": case["time"], "sample_time": target_time,
                                     "lon": case["lon"], "lat": case["lat"], "error": repr(exc)})
                    logging.exception("Case failed: %s %.5f %.5f", category, case["lon"], case["lat"])
        fields = None
        processor._load_file_cubes.cache_clear()
        gc.collect()

        next_group = group_index + 1
        if next_group % args.checkpoint_every == 0 or next_group == len(time_groups):
            save_checkpoint(checkpoint_path, signature, next_group, sums, counts, successful, attempted)
            pd.DataFrame(failures, columns=FAILURE_COLUMNS).to_csv(failures_path, index=False)
            logging.info("Checkpoint saved: group %d; successes=%s; failures=%d", next_group, successful, len(failures))

    means = means_from_accumulators(sums, counts)
    nc_path = run_dir / f"umbrella_CI_composites_{lag_tag(args.lag_hours)}.nc"
    write_netcdf(nc_path, processor, means, counts, successful, args)
    pd.DataFrame({"category": DIRECTION_ORDER,
                  "attempted": [attempted[c] for c in DIRECTION_ORDER],
                  "successful": [successful[c] for c in DIRECTION_ORDER],
                  "failed": [attempted[c] - successful[c] for c in DIRECTION_ORDER]}).to_csv(run_dir / "processing_summary.csv", index=False)

    if not args.no_plots:
        plots = [
            ("sm_anomaly", "Top-layer soil-moisture spatial anomaly", "BrBG", "kg m$^{-2}$", True, 1, False),
            ("olr", "Outgoing longwave radiation", "Greys_r", "W m$^{-2}$", False, 1, False),
            ("q_low", "Near-surface specific humidity", "YlGnBu", "g kg$^{-1}$", False, 1, False),
            ("rain", "Preceding-hour rainfall", "Blues", "mm h$^{-1}$", False, 1, True),
            ("shear", "600 hPa-10 m shear magnitude", "magma", "m s$^{-1}$", False, 1, True),
            ("shear_u", "Rotated zonal shear ($u_{600}-u_{10}$)", "RdBu_r", "m s$^{-1}$", True, 1, False),
            ("shear_v", "Rotated meridional shear ($v_{600}-v_{10}$)", "RdBu_r", "m s$^{-1}$", True, 1, False),
            ("div10", "10 m divergence", "RdBu_r", "$10^{-5}$ s$^{-1}$", True, 1e5, False),
        ]
        for variable, title, cmap, label, symmetric, factor, floor_zero in plots:
            plot_composite(run_dir / f"{variable}_{lag_tag(args.lag_hours)}.png", processor, means, successful,
                           variable, title, cmap, label, symmetric, factor, floor_zero)

    logging.info("Finished. NetCDF: %s", nc_path)


def main() -> None:
    args = parse_args()
    configure_logging()

    cases = read_ci_files(args.ci_dir, args.year, args.drop_exact_duplicates)
    cases = cases[cases["hour"].isin(args.ci_hours)].copy()
    if cases.empty:
        raise RuntimeError(f"No CIs found for hours {args.ci_hours}")

    metrics = calculate_location_metrics(
        cases, args.orography, args.landmask,
        args.water_radius_km, args.topo_radius_km,
    )
    cases = cases.merge(metrics, on=["lon", "lat"], how="left")
    cases["coast_ok"] = cases["central_land"] & (cases["water_distance_km"] >= args.water_radius_km)
    cases["topography_ok"] = (
        np.isfinite(cases["topo_relief_200km"]) &
        (cases["topo_relief_200km"] <= args.topo_relief_max)
    )
    cases["sample_ok"] = cases["coast_ok"] & cases["topography_ok"]

    for lag in args.lag_hours:
        lag_args = argparse.Namespace(**vars(args))
        lag_args.lag_hours = int(lag)
        logging.info("Starting composites for CI hours %s at t=%+d h", args.ci_hours, lag)
        run_one_lag(lag_args, cases)


if __name__ == "__main__":
    main()
