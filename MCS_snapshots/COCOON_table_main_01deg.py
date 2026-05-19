from MCS_snapshots import COCOON_table_create_01deg

import os
import glob
import multiprocessing
import traceback

import numpy as np
import pandas as pd
import xarray as xr
import pyproj

from utils import constants as cnst
from LMCS import glob_util


REGIONS = glob_util.REGIONS

YEARS = range(2000, 2021)
NPROC = 4

LMCS = cnst.lmcs_drive + "/MCS_Feng/global_v2/2d_fields/"
OUTROOT = cnst.lmcs_drive + "/MCS_5000km2_tables_v3/base_tables_tir-prcp/"

T_THRESH = -40
MIN_MCS_SIZE = 5000
DATA_RES_KM = 10

OUT_SUFFIX = "MCS_5000km2_-40C_0.1degTIR-IMERG_hourly.parquet"


def pixel_area_latlon_km2(lat, lon):
    """
    Approximate pixel area for a regular lat/lon grid.

    lat, lon are 1D arrays in degrees.
    Returns 2D area array in km2 with shape (lat, lon).
    """
    R = 6371.0  # km

    dlat = np.abs(np.gradient(lat))
    dlon = np.abs(np.gradient(lon))

    lat_rad = np.deg2rad(lat)
    dlat_rad = np.deg2rad(dlat)
    dlon_rad = np.deg2rad(dlon)

    area = (
        R**2
        * dlat_rad[:, None]
        * dlon_rad[None, :]
        * np.cos(lat_rad[:, None])
    )

    return area


def make_table(reg):
    """
    Create yearly MCS baseline tables for one region.

    Output is one parquet file per region/year.
    """
    print(f"Doing region {reg}")

    box = REGIONS[reg][0]
    outdir = os.path.join(OUTROOT, reg)
    os.makedirs(outdir, exist_ok=True)

    transformer = pyproj.Transformer.from_crs(
        "EPSG:4326",
        "EPSG:6933",  # equal-area, metres
        always_xy=True,
    )

    # Cache this after reading the first valid file.
    # It should remain constant for a given region if all input files share the same grid.
    area_grid_km2 = None
    cached_lat = None
    cached_lon = None

    for yy in YEARS:
        print(f"Doing {reg} {yy}")

        outfile = os.path.join(outdir, f"{yy}_{OUT_SUFFIX}")

        if os.path.isfile(outfile):
            print(outfile, "exists, continue")
            continue

        infiles = sorted(glob.glob(os.path.join(LMCS, f"{yy}*", "*.nc")))

        if len(infiles) == 0:
            print(f"No files found for {reg} {yy}")
            continue

        year_dfs = []

        for infile in infiles:
            print("Doing", infile)

            try:
                with xr.open_dataset(infile) as ds:
                    ds = ds.sel(
                        lon=slice(box[0], box[1]),
                        lat=slice(box[2], box[3]),
                    ).load()

                if ds.lat.size == 0 or ds.lon.size == 0:
                    print("Empty regional subset, continue")
                    continue

                # Build / update area grid only if needed.
                # This avoids recomputing it for every file.
                lat_now = ds.lat.values
                lon_now = ds.lon.values

                if (
                    area_grid_km2 is None
                    or cached_lat is None
                    or cached_lon is None
                    or len(lat_now) != len(cached_lat)
                    or len(lon_now) != len(cached_lon)
                    or not np.allclose(lat_now, cached_lat)
                    or not np.allclose(lon_now, cached_lon)
                ):
                    print("Computing area grid for", reg, yy)
                    area_grid_km2 = pixel_area_latlon_km2(lat_now, lon_now)
                    cached_lat = lat_now.copy()
                    cached_lon = lon_now.copy()

                basic_tab = COCOON_table_create_01deg.process_tir_image(
                    ds,
                    DATA_RES_KM,
                    t_thresh=T_THRESH,
                    min_mcs_size=MIN_MCS_SIZE,
                    area_grid_km2=area_grid_km2,
                    transformer=transformer,
                    rainvar_name="precipitation",
                    min_rain_pmax=1,
                )

                merge_tab = COCOON_table_create_01deg.add_environment_toTable(
                    basic_tab,
                    ds,
                    DATA_RES_KM,
                    envvar_take=[],
                    rainvar_name="precipitation",
                    area_grid_km2=area_grid_km2,
                )

                merge_tab.pop("cloudMask", None)
                merge_tab.pop("tir", None)

                if len(merge_tab["date"]) == 0:
                    continue

                df = pd.DataFrame.from_dict(merge_tab)
                df["region"] = reg
                year_dfs.append(df)

                print("Did", infile)

            except Exception as e:
                print("Failed file:", infile)
                print(e)
                traceback.print_exc()
                continue

        if len(year_dfs) == 0:
            print(f"No storms found for {reg} {yy}")
            continue

        pd_out = pd.concat(year_dfs, ignore_index=True)
        pd_out.to_parquet(outfile, index=False)

        print(f"Wrote {outfile} with {len(pd_out)} storms")


if __name__ == "__main__":
    regs = list(REGIONS.keys())

    with multiprocessing.Pool(processes=NPROC) as pool:
        pool.map(make_table, regs)