import os
import glob
import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr
import multiprocessing as mp

from utils import u_darrays
from GLOBAL import glob_util

# ============================================================
# CONFIG
# ============================================================

REGIONS = glob_util.REGIONS

ERA_ROOT = "/prj/global_water/ERA5_global_0.7/hourly"
BASE_TABLE_DIR = "/prj/global_water/MCS_5000km2_tables_v2/base_tables_tir-prcp"
STATIC_ERA_FILE = ("/prj/global_water/ERA5_global_0.7/static/era5_invariant_07deg_COCOON-study_regrid.nc")


TEST_MODE = True

if TEST_MODE:
    YEARS = [2006]
    TEST_MONTHS = [8]
    TEST_REGIONS = ["GPlains"]
    TEST_N_STORMS_PER_REGION = 1000
    REL_HOURS = np.arange(-6, 4)
    NPROC = 3

    OUT_ROOT = (
        "/prj/global_water/MCS_5000km2_tables_v2/"
        "ERA5_storm_timeseries_TEST"
    )

else:
    YEARS = range(2000, 2021)
    TEST_REGIONS = None
    TEST_N_STORMS_PER_REGION = None
    REL_HOURS = np.arange(-6, 4)
    NPROC = 4

    OUT_ROOT = (
        "/prj/global_water/MCS_5000km2_tables_v2/"
        "ERA5_storm_timeseries"
    )

BOX_RADIUS = 0.35
LOCAL_STORM_TIME = [12,21]


# ============================================================
# VARIABLE CONFIG
# ============================================================

# ============================================================
# ERA INPUT VARIABLES
# ============================================================

ERA_INPUTS = {

    # --------------------------------------------------------
    # SURFACE
    # --------------------------------------------------------

    "tcwv": {
        "kind": "surface",
        "era_name": "total_column_water_vapour",
        "nc_var": "tcwv",
    },

    "cape": {
        "kind": "surface",
        "era_name": "convective_available_potential_energy",
        "nc_var": "cape",
    },

    "mtpr": {
        "kind": "surface",
        "era_name": "mean_total_precipitation_rate",
        "nc_var": "mtpr",
    },
    "qdiv_flux": {
        "kind": "surface",
        "era_name": "vertical_integral_of_divergence_of_moisture_flux",
        "nc_var": "p84.162",
    },

    "tciw": {
        "kind": "surface",
        "era_name": "total_column_cloud_ice_water",
        "nc_var": "tciw",
    },
    "sp": {
        "kind": "surface",
        "era_name": "surface_pressure",
        "nc_var": "sp",
    },

    "cin": {
        "kind": "surface",
        "era_name": "convective_inhibition",
        "nc_var": "cin",
    },

    "t2m": {
        "kind": "surface",
        "era_name": "2m_temperature",
        "nc_var": "t2m",
    },

    "d2m": {
        "kind": "surface",
        "era_name": "2m_dewpoint_temperature",
        "nc_var": "d2m",
    },

    "u100": {
        "kind": "surface",
        "era_name": "100m_u_component_of_wind",
        "nc_var": "u100",
    },

    "v100": {
        "kind": "surface",
        "era_name": "100m_v_component_of_wind",
        "nc_var": "v100",
    },


    # --------------------------------------------------------
    # PRESSURE LEVELS
    # --------------------------------------------------------

    # --------------------
    # Specific humidity
    # --------------------

    "q925": {
        "kind": "pressure_levels",
        "era_name": "specific_humidity",
        "nc_var": "q",
        "level": 925,
    },

    "q850": {
        "kind": "pressure_levels",
        "era_name": "specific_humidity",
        "nc_var": "q",
        "level": 850,
    },

    "q750": {
        "kind": "pressure_levels",
        "era_name": "specific_humidity",
        "nc_var": "q",
        "level": 750,
    },

    "q650": {
        "kind": "pressure_levels",
        "era_name": "specific_humidity",
        "nc_var": "q",
        "level": 650,
    },

    "q500": {
        "kind": "pressure_levels",
        "era_name": "specific_humidity",
        "nc_var": "q",
        "level": 500,
    },


    # --------------------
    # Relative humidity
    # --------------------

    "rh925": {
        "kind": "pressure_levels",
        "era_name": "relative_humidity",
        "nc_var": "r",
        "level": 925,
    },

    "rh850": {
        "kind": "pressure_levels",
        "era_name": "relative_humidity",
        "nc_var": "r",
        "level": 850,
    },

    "rh750": {
        "kind": "pressure_levels",
        "era_name": "relative_humidity",
        "nc_var": "r",
        "level": 750,
    },

    "rh650": {
        "kind": "pressure_levels",
        "era_name": "relative_humidity",
        "nc_var": "r",
        "level": 650,
    },

    "rh500": {
        "kind": "pressure_levels",
        "era_name": "relative_humidity",
        "nc_var": "r",
        "level": 500,
    },


    # --------------------
    # Temperature
    # --------------------

    "t925": {
        "kind": "pressure_levels",
        "era_name": "temperature",
        "nc_var": "t",
        "level": 925,
    },

    "t850": {
        "kind": "pressure_levels",
        "era_name": "temperature",
        "nc_var": "t",
        "level": 850,
    },

    "t750": {
        "kind": "pressure_levels",
        "era_name": "temperature",
        "nc_var": "t",
        "level": 750,
    },

    "t650": {
        "kind": "pressure_levels",
        "era_name": "temperature",
        "nc_var": "t",
        "level": 650,
    },

    "t500": {
        "kind": "pressure_levels",
        "era_name": "temperature",
        "nc_var": "t",
        "level": 500,
    },



    # --------------------
    # w wind
    # --------------------


    "w850": {
        "kind": "pressure_levels",
        "era_name": "vertical_velocity",
        "nc_var": "w",
        "level": 850,
    },


    "w650": {
        "kind": "pressure_levels",
        "era_name": "vertical_velocity",
        "nc_var": "w",
        "level": 650,
    },

    "w500": {
        "kind": "pressure_levels",
        "era_name": "vertical_velocity",
        "nc_var": "w",
        "level": 500,
    },

    "d850": {
        "kind": "pressure_levels",
        "era_name": "divergence",
        "nc_var": "d",
        "level": 850,
    },
    

    # --------------------
    # U wind
    # --------------------

    "u925": {
        "kind": "pressure_levels",
        "era_name": "u_component_of_wind",
        "nc_var": "u",
        "level": 925,
    },

    "u850": {
        "kind": "pressure_levels",
        "era_name": "u_component_of_wind",
        "nc_var": "u",
        "level": 850,
    },

    "u750": {
        "kind": "pressure_levels",
        "era_name": "u_component_of_wind",
        "nc_var": "u",
        "level": 750,
    },

    "u650": {
        "kind": "pressure_levels",
        "era_name": "u_component_of_wind",
        "nc_var": "u",
        "level": 650,
    },


    # --------------------
    # V wind
    # --------------------

    "v925": {
        "kind": "pressure_levels",
        "era_name": "v_component_of_wind",
        "nc_var": "v",
        "level": 925,
    },

    "v850": {
        "kind": "pressure_levels",
        "era_name": "v_component_of_wind",
        "nc_var": "v",
        "level": 850,
    },

    "v750": {
        "kind": "pressure_levels",
        "era_name": "v_component_of_wind",
        "nc_var": "v",
        "level": 750,
    },

    "v650": {
        "kind": "pressure_levels",
        "era_name": "v_component_of_wind",
        "nc_var": "v",
        "level": 650,
    },
}


# ============================================================
# OUTPUT VARIABLES
# ============================================================

OUTPUT_VARIABLES = {

    # --------------------------------------------------------
    # DIRECTLY SAVED VARIABLES
    # --------------------------------------------------------

    "tcwv": {
        "type": "direct",
        "inputs": ["tcwv"],
    },

    "cape": {
        "type": "direct",
        "inputs": ["cape"],
    },
    "era_prcp": {
        "type": "direct",
        "inputs": ["mtpr"],
    },
    "era_ice": {
        "type": "direct",
        "inputs": ["tciw"],
    },
    "qdiv_flux": {
        "type": "direct",
        "inputs": ["qdiv_flux"],
    },

    "div850": {
        "type": "direct",
        "inputs": ["d850"],
    },
    
    "p_srfc": {
        "type": "direct",
        "inputs": ["sp"],
    },

    "cin": {
        "type": "direct",
        "inputs": ["cin"],
    },

    "t2m": {
        "type": "direct",
        "inputs": ["t2m"],
    },

    "d2m": {
        "type": "direct",
        "inputs": ["d2m"],
    },

    "q925": {
        "type": "direct",
        "inputs": ["q925"],
    },

    "q850": {
        "type": "direct",
        "inputs": ["q850"],
    },

    "q750": {
        "type": "direct",
        "inputs": ["q750"],
    },

    "q650": {
        "type": "direct",
        "inputs": ["q650"],
    },

    "q500": {
        "type": "direct",
        "inputs": ["q500"],
    },

    "rh925": {
        "type": "direct",
        "inputs": ["rh925"],
    },

    "rh850": {
        "type": "direct",
        "inputs": ["rh850"],
    },

    "rh750": {
        "type": "direct",
        "inputs": ["rh750"],
    },

    "rh650": {
        "type": "direct",
        "inputs": ["rh650"],
    },

    "rh500": {
        "type": "direct",
        "inputs": ["rh500"],
    },

    "t925": {
        "type": "direct",
        "inputs": ["t925"],
    },

    "t850": {
        "type": "direct",
        "inputs": ["t850"],
    },

    "t750": {
        "type": "direct",
        "inputs": ["t750"],
    },

    "t650": {
        "type": "direct",
        "inputs": ["t650"],
    },

    "t500": {
        "type": "direct",
        "inputs": ["t500"],
    },

    "w850": {
        "type": "direct",
        "inputs": ["w850"],
    },

    "w650": {
        "type": "direct",
        "inputs": ["w650"],
    },

    "w500": {
        "type": "direct",
        "inputs": ["w500"],
    },


    # --------------------------------------------------------
    # DERIVED 2M DIAGNOSTICS
    # --------------------------------------------------------

    "rh2m": {
        "type": "rh2m",
        "inputs": ["t2m", "d2m"],
    },

    "vpd2m": {
        "type": "vpd2m",
        "inputs": ["t2m", "d2m"],
    },


    # --------------------------------------------------------
    # SHEAR VARIABLES
    # --------------------------------------------------------

    "ushear100_650": {
        "type": "ushear100_650",
        "inputs": ["u650", "u100"],
    },

    "vshear100_650": {
        "type": "vshear100_650",
        "inputs": ["v650", "v100"],
    },

    "shear100_650": {
        "type": "shear100_650",
        "inputs": ["u650", "v650", "u100", "v100"],
    },


    # --------------------------------------------------------
    # OPTIONAL: DEEPER SHEAR LAYERS
    # --------------------------------------------------------

    "ushear925_650": {
        "type": "ushear925_650",
        "inputs": ["u650", "u925"],
    },

    "vshear925_650": {
        "type": "vshear925_650",
        "inputs": ["v650", "v925"],
    },

    "shear925_650": {
        "type": "shear925_650",
        "inputs": ["u650", "v650", "u925", "v925"],
    },
}

# ============================================================
# ERA5 FILE HANDLING
# ============================================================

def era_file_path(kind, era_name, edate):
    suffix = "srfc" if kind == "surface" else "pl"

    return (
        f"{ERA_ROOT}/{kind}/{era_name}/"
        f"{era_name}_ERA5_{edate.year}_{edate.month:02d}_{edate.day:02d}_{suffix}.nc"
    )


def open_era_input(input_name, edate):
    spec = ERA_INPUTS[input_name]
    path = era_file_path(spec["kind"], spec["era_name"], edate)

    if not os.path.isfile(path):
        print("Missing ERA5:", path)
        return None

    try:
        ds = xr.open_dataset(path)
        ds = u_darrays.flip_lat(ds)
        ds = ds.sel(time=edate, method="nearest")
        return ds

    except Exception as e:
        print("Could not open/select:", path)
        print(e)
        return None


# ============================================================
# OCEAN MASK
# ============================================================


OCEAN_REGIONS = ["Atl", "Pcf", "InO"]

LAND_THRESHOLD = 0.6


def filter_land_ocean_storms(storm_df):

    before = len(storm_df)

    keep = np.ones(len(storm_df), dtype=bool)

    for i, row in storm_df.iterrows():

        region = row["region"]
        lsm = row["lsm"]

        # ----------------------------------------------------
        # Ocean regions:
        # keep ocean storms only
        # ----------------------------------------------------

        if region in OCEAN_REGIONS:

            keep[i] = lsm < LAND_THRESHOLD

        # ----------------------------------------------------
        # Land regions:
        # keep land storms only
        # ----------------------------------------------------

        else:

            keep[i] = lsm >= LAND_THRESHOLD


    storm_df = storm_df.loc[keep].reset_index(drop=True)

    after = len(storm_df)

    print(
        "Land/ocean filter:",
        before,
        "->",
        after,
    )

    return storm_df

# ============================================================
# STATIC VAR SAMPLING
# ============================================================

def open_static_era():
    ds = xr.open_dataset(STATIC_ERA_FILE)
    ds = u_darrays.flip_lat(ds)

    print("STATIC coords:")
    print("lat min/max:", float(ds.latitude.min()), float(ds.latitude.max()))
    print("lon min/max:", float(ds.longitude.min()), float(ds.longitude.max()))

    return ds


def wrap_lon_to_dataset(lon, ds):
    lon_min = float(ds.longitude.min())
    lon_max = float(ds.longitude.max())

    if lon_min >= 0 and lon < 0:
        return lon % 360
    if lon_max <= 180 and lon > 180:
        return ((lon + 180) % 360) - 180

    return lon


def sample_static_nearest(ds_static, lat, lon):
    out = {}

    # adapt names if your netCDF variables differ slightly
    static_vars = [
        "lsm",
        "topographic_complexity",
        "topographic_complexity_sdfor_slor",
    ]

    for v in static_vars:
        
        if v not in ds_static:
            print("Static variable missing:", v, "available:", list(ds_static.data_vars))
            out[v] = np.nan
            continue

        try:
            out[v] = float(
                ds_static[v].sel(
                    latitude=lat,
                    longitude=lon,
                    method="nearest",
                ).squeeze()
            )
        except Exception:
            out[v] = np.nan

    return out


def add_static_fields_to_storms(storm_df):
    ds_static = open_static_era()

    rows = []

    for _, row in storm_df.iterrows():
        vals = sample_static_nearest(
            ds_static,
            row["tminlat"],
            row["tminlon"],
        )
        rows.append(vals)

    static_df = pd.DataFrame(rows, index=storm_df.index)

    storm_df = pd.concat(
        [storm_df.reset_index(drop=True), static_df.reset_index(drop=True)],
        axis=1,
    )

    return storm_df

# ============================================================
# TIME AND AGGREGATION
# ============================================================


def month_in_season(month, season):
    start, end = season

    if start <= end:
        return start <= month <= end
    else:
        return (month >= start) or (month <= end)


def box_mean(ds, nc_var, lat, lon, level=None):
    if ds is None:
        return np.nan

    if nc_var not in ds:
        print("Variable not found:", nc_var, "available:", list(ds.data_vars))
        return np.nan

    try:
        da = ds[nc_var]

        if level is not None:
            da = da.sel(level=level)

        val = da.sel(
            latitude=slice(lat - BOX_RADIUS, lat + BOX_RADIUS),
            longitude=slice(lon - BOX_RADIUS, lon + BOX_RADIUS),
        ).mean(["latitude", "longitude"])

        return float(val)

    except Exception:
        return np.nan


def nearest_pixel(ds, nc_var, lat, lon, level=None):
    if ds is None:
        return np.nan

    if nc_var not in ds:
        print("Variable not found:", nc_var, "available:", list(ds.data_vars))
        return np.nan

    try:
        da = ds[nc_var]

        if level is not None:
            da = da.sel(level=level)

        val = da.sel(
            latitude=lat,
            longitude=lon,
            method="nearest",
        )

        return float(val)

    except Exception:
        return np.nan


# ============================================================
# ON THE FLY CALCULATIONS
# ============================================================


def calc_rh_from_dewpoint(t_c, td_c):

    es = 6.112 * np.exp((17.67 * t_c) / (t_c + 243.5))
    e = 6.112 * np.exp((17.67 * td_c) / (td_c + 243.5))

    return 100.0 * e / es


def calc_vpd_from_dewpoint(t2m_k, d2m_k):
    t_c = t2m_k - 273.15
    td_c = d2m_k - 273.15

    es = 6.112 * np.exp((17.67 * t_c) / (t_c + 243.5))
    e = 6.112 * np.exp((17.67 * td_c) / (td_c + 243.5))

    # hPa
    return es - e




# ============================================================
# STORM TABLE HANDLING
# ============================================================

def find_base_table(region, year):
    pattern = f"{BASE_TABLE_DIR}/{region}/{year}_MCS_5000km2_*.csv"
    files = glob.glob(pattern)

    if len(files) == 0:
        return None

    if len(files) > 1:
        print("Multiple base files found for", region, year, "using first:")
        for f in files:
            print(" ", f)

    return files[0]


def read_base_table(region, year):

    path = find_base_table(region, year)

    if path is None:
        print("No base table:", region, year)
        return None

    print("Reading:", path)

    df = pd.read_csv(path).reset_index(drop=True)

    df["region"] = region

    df["storm_id"] = (
        region + "_" + str(year) + "_" + df.index.astype(str).str.zfill(7)
    )

    storm_utc = []
    storm_lt = []

    for y, m, d, h in zip(df["year"], df["month"], df["day"], df["hour"]):

        # Base-table storm time is UTC
        utc_dt = pd.Timestamp(
            dt.datetime(
                int(y),
                int(m),
                int(d),
                int(h),
            )
        )

        lt_dt = glob_util.UTC_to_LT_date(utc_dt, region)

        storm_utc.append(utc_dt)
        storm_lt.append(lt_dt)

    df["storm_utc_time"] = pd.to_datetime(storm_utc)
    df["storm_lt_time"] = pd.to_datetime(storm_lt)

    df["storm_utc_hour"] = df["storm_utc_time"].dt.hour
    df["storm_lt_hour"] = df["storm_lt_time"].dt.hour

    # =====================================================
    # Keep only local afternoon/evening storms
    # =====================================================

    before = len(df)

    df = df[
        df["storm_lt_hour"].between(LOCAL_STORM_TIME[0], LOCAL_STORM_TIME[1])
    ].reset_index(drop=True)

    after = len(df)

    print(
        region,
        year,
        "local 15-21 LT storm filter:",
        before,
        "->",
        after,
    )

    ### season filter - only sample months of interest
    season = REGIONS[region][3]
    
    before = len(df)
    
    df = df[
        df["month"].apply(lambda m: month_in_season(int(m), season))
    ].reset_index(drop=True)
    
    after = len(df)
    
    print(
        region,
        year,
        "season filter",
        season,
        ":",
        before,
        "->",
        after,
    )

    
    if TEST_MODE and TEST_MONTHS is not None:
    
        before = len(df)
    
        df = df[
            df["month"].isin(TEST_MONTHS)
        ].reset_index(drop=True)
    
        after = len(df)
    
        print(
            region,
            year,
            "test month filter:",
            before,
            "->",
            after,
        )

    return df


def make_request_table(storm_df, region):
    rows = []

    for _, s in storm_df.iterrows():
        for rel_hour in REL_HOURS:
            era_time = s["storm_utc_time"] + pd.Timedelta(hours=int(rel_hour))

            rows.append({
                "region": region,
                "storm_id": s["storm_id"],
                "storm_utc_time": s["storm_utc_time"],
                "storm_lt_time": s["storm_lt_time"],
                "era_time": era_time.to_pydatetime(),
                "rel_hour": int(rel_hour),
                "year": int(s["year"]),
                "month": int(s["month"]),
                "day": int(s["day"]),
                "hour": int(s["hour"]),
                "tminlat": float(s["tminlat"]),
                "tminlon": float(s["tminlon"]),
            })

    return pd.DataFrame(rows)


def build_request_table_all_regions(year):

    all_requests = []

    if TEST_MODE:
        region_names = TEST_REGIONS
    else:
        region_names = list(REGIONS.keys())

    for region in region_names:

        storm_df = read_base_table(region, year)

        if storm_df is None or len(storm_df) == 0:
            continue

        # ============================================
        # TEST MODE: only keep first few storms
        # ============================================

        if TEST_MODE:
            storm_df = storm_df.head(TEST_N_STORMS_PER_REGION)

            print(
                "TEST MODE:",
                region,
                "keeping",
                len(storm_df),
                "storms"
            )

        req = make_request_table(storm_df, region)

        all_requests.append(req)

    if len(all_requests) == 0:
        return None

    request_df = pd.concat(all_requests, ignore_index=True)

    print("Final request table size:", len(request_df))

    return request_df


def save_storm_metadata(storm_df, region, year):
    outdir = f"{OUT_ROOT}/storm_metadata/region={region}/year={year}"
    os.makedirs(outdir, exist_ok=True)

    outfile = f"{outdir}/part.parquet"
    storm_df.to_parquet(outfile, index=False)

    print("Saved storm metadata:", outfile, "rows:", len(storm_df))


# ============================================================
# EXTRACTION
# ============================================================

def extract_for_one_era_time(group):

    era_time = pd.Timestamp(group["era_time"].iloc[0]).to_pydatetime()

    print("Extracting:", era_time, "requests:", len(group))

    needed_inputs = sorted({
        inp
        for outspec in OUTPUT_VARIABLES.values()
        for inp in outspec["inputs"]
    })

    opened = {
        inp: open_era_input(inp, era_time)
        for inp in needed_inputs
    }

    output_rows = {outvar: [] for outvar in OUTPUT_VARIABLES}

    for _, req in group.iterrows():

        lat = req["tminlat"]
        lon = req["tminlon"]

        sampled = {}

        for inp in needed_inputs:
            spec = ERA_INPUTS[inp]

            sampled[inp] = nearest_pixel(
                opened[inp],
                spec["nc_var"],
                lat,
                lon,
                level=spec.get("level", None),
            )

        for outvar, outspec in OUTPUT_VARIABLES.items():

            otype = outspec["type"]

            if otype == "direct":
                value = sampled[outspec["inputs"][0]]

            elif otype == "rh2m":
                t2m = sampled["t2m"]
                d2m = sampled["d2m"]

                if np.isfinite(t2m) and np.isfinite(d2m):
                    value = calc_rh_from_dewpoint(t2m - 273.15, d2m - 273.15)
                else:
                    value = np.nan

            elif otype == "vpd2m":
                t2m = sampled["t2m"]
                d2m = sampled["d2m"]

                if np.isfinite(t2m) and np.isfinite(d2m):
                    value = calc_vpd_from_dewpoint(t2m, d2m)
                else:
                    value = np.nan

            elif otype == "ushear100_650":
                value = (
                    sampled["u650"] - sampled["u100"]
                    if np.isfinite(sampled["u650"]) and np.isfinite(sampled["u100"])
                    else np.nan
                )

            elif otype == "vshear100_650":
                value = (
                    sampled["v650"] - sampled["v100"]
                    if np.isfinite(sampled["v650"]) and np.isfinite(sampled["v100"])
                    else np.nan
                )

            elif otype == "shear100_650":
                if all(np.isfinite(sampled[x]) for x in ["u650", "v650", "u100", "v100"]):
                    ushear = sampled["u650"] - sampled["u100"]
                    vshear = sampled["v650"] - sampled["v100"]
                    value = np.sqrt(ushear ** 2 + vshear ** 2)
                else:
                    value = np.nan

            elif otype == "ushear925_650":
                value = (
                    sampled["u650"] - sampled["u925"]
                    if np.isfinite(sampled["u650"]) and np.isfinite(sampled["u925"])
                    else np.nan
                )

            elif otype == "vshear925_650":
                value = (
                    sampled["v650"] - sampled["v925"]
                    if np.isfinite(sampled["v650"]) and np.isfinite(sampled["v925"])
                    else np.nan
                )

            elif otype == "shear925_650":
                if all(np.isfinite(sampled[x]) for x in ["u650", "v650", "u925", "v925"]):
                    ushear = sampled["u650"] - sampled["u925"]
                    vshear = sampled["v650"] - sampled["v925"]
                    value = np.sqrt(ushear ** 2 + vshear ** 2)
                else:
                    value = np.nan

            else:
                raise ValueError(f"Unknown output variable type: {otype}")

            storm_utc_time = pd.Timestamp(req["storm_utc_time"])
            storm_lt_time = pd.Timestamp(req["storm_lt_time"])
            era_utc_time = pd.Timestamp(req["era_time"])
            era_lt_time = glob_util.UTC_to_LT_date(era_utc_time, req["region"])

            output_rows[outvar].append({
                "region": req["region"],
                "storm_id": req["storm_id"],

                "storm_utc_time": storm_utc_time,
                "storm_lt_time": storm_lt_time,
                "era_utc_time": era_utc_time,
                "era_lt_time": era_lt_time,

                "storm_utc_hour": storm_utc_time.hour,
                "storm_lt_hour": storm_lt_time.hour,
                "era_utc_hour": era_utc_time.hour,
                "era_lt_hour": era_lt_time.hour,

                "rel_hour": req["rel_hour"],

                "year": req["year"],
                "month": req["month"],
                "day": req["day"],
                "hour": req["hour"],

                "tminlat": req["tminlat"],
                "tminlon": req["tminlon"],

                "variable": outvar,
                "value": value,
            })

    return {
        outvar: pd.DataFrame(rows)
        for outvar, rows in output_rows.items()
    }


# ============================================================
# SAVING
# ============================================================

def variable_output_file(varname, region, year):

    outdir = (
        f"{OUT_ROOT}/{varname}/"
        f"region={region}/year={year}"
    )

    os.makedirs(outdir, exist_ok=True)

    return f"{outdir}/part.parquet"


def save_variable_region_year(varname, df, region, year):

    outfile = variable_output_file(varname, region, year)

    df.to_parquet(outfile, index=False)

    print("Saved:", outfile, "rows:", len(df))


def output_exists(varname, region, year):

    outfile = variable_output_file(varname, region, year)

    return os.path.isfile(outfile)


def save_all_outputs(collected, year):
    for varname, frames in collected.items():
        if len(frames) == 0:
            continue

        df = pd.concat(frames, ignore_index=True)

        for region, rchunk in df.groupby("region"):
            save_variable_region_year(varname, rchunk, region, year)


# ============================================================
# YEAR DRIVER
# ============================================================

def process_year_all_regions(year):

    print("\n==============================")
    print("Processing year:", year)
    print("==============================")

    # =====================================================
    # Determine which outputs are missing
    # =====================================================

    missing = []

    region_names = (
        TEST_REGIONS if TEST_MODE
        else list(REGIONS.keys())
    )

    for varname in OUTPUT_VARIABLES:

        for region in region_names:

            if not output_exists(varname, region, year):

                missing.append((varname, region))

    # -----------------------------------------------------

    if len(missing) == 0:

        print("All outputs already exist for year", year)
        return

    print("\nMissing outputs:\n")

    for m in missing:
        print(m)

    needed_regions = sorted(set(r for _, r in missing))

    print("\nNeed to process regions:", needed_regions)

    # =====================================================
    # Build request table only for needed regions
    # =====================================================

    all_requests = []

    for region in needed_regions:

        storm_df = read_base_table(region, year)

        if storm_df is None or len(storm_df) == 0:
            continue

        if TEST_MODE:
            storm_df = storm_df.head(TEST_N_STORMS_PER_REGION)

            print(
                "TEST MODE:",
                region,
                "keeping",
                len(storm_df),
                "storms"
            )

        storm_df = add_static_fields_to_storms(storm_df)
        
        storm_df = filter_land_ocean_storms(storm_df)
        
        save_storm_metadata(storm_df, region, year)
        
        req = make_request_table(storm_df, region)

        all_requests.append(req)

    if len(all_requests) == 0:

        print("No requests for year:", year)
        return

    request_df = pd.concat(all_requests, ignore_index=True)

    print("\nRequest table preview:\n")

    print(
        request_df[
            [
                "region",
                "storm_id",
                "storm_utc_time",
                "era_time",
                "rel_hour",
            ]
        ].head(20)
    )

    print("\nTotal requests:", len(request_df))
    print("Unique ERA times:", request_df["era_time"].nunique())

    # =====================================================
    # Group by ERA time
    # =====================================================

    grouped = [
        g.copy()
        for _, g in request_df.groupby("era_time")
    ]

    print("Number of ERA groups:", len(grouped))

    # =====================================================
    # Run extraction
    # =====================================================

    if NPROC == 1:
        results = [extract_for_one_era_time(g) for g in grouped]

    else:
        with mp.Pool(processes=NPROC) as pool:
            results = pool.map(extract_for_one_era_time, grouped)

    # =====================================================
    # Collect results
    # =====================================================

    collected = {var: [] for var in OUTPUT_VARIABLES}

    for res in results:

        for var in OUTPUT_VARIABLES:

            collected[var].append(res[var])

    # =====================================================
    # Save only missing outputs
    # =====================================================

    for varname, frames in collected.items():

        if len(frames) == 0:
            continue

        df = pd.concat(frames, ignore_index=True)

        for region, rchunk in df.groupby("region"):

            if (varname, region) not in missing:

                print(
                    "Output already exists, skipping:",
                    varname,
                    region,
                    year
                )

                continue

            save_variable_region_year(
                varname,
                rchunk,
                region,
                year,
            )

    print("Finished year:", year)


def main():
    for year in YEARS:
        process_year_all_regions(year)


if __name__ == "__main__":
    main()