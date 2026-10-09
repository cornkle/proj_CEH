#!/usr/bin/env python3
"""
Compare occurrence of strongly negative SM anomalies on MCS days with their
background occurrence at MCS locations.

Design
------
1. Use the same eligible MCS sample as the current DRY/WET analysis:
      |v_mid_fullBox| < 50 m s-1
      storm area < 500,000 km2
      pmax >= 1 mm h-1

2. Map all eligible storms to their CP4 environmental-grid cell and retain
   each spatial location only once for the background sampling.

3. For EVERY JAS day, sample the 10-12 UTC mean SM anomaly at every unique
   MCS location, using the same 2-degree box geometry and noon-rain validity
   mask used for storm SM_2deg.

4. Save those daily values YEAR BY YEAR:
      background_SM_daily_hist_1998.nc
      ...
      background_SM_daily_fut_2006.nc

   Completed years are automatically reused on subsequent runs.

5. Obtain each storm's independently sampled SM value by looking up its
   date + location in the SAME daily cache.

6. Compare:
      P(DRY | storm)
      P(DRY | all location-days)
   in HIST and FUT.

Changing DRY_THRESHOLD does NOT require rebuilding the daily caches.
"""

from pathlib import Path
import glob
import hashlib
import os
import re
import sys

import numpy as np
import pandas as pd
import xarray as xr


# =============================================================================
# IMPORTS / PATHS
# =============================================================================

REPO_ROOT = Path("/home/users/cornkle/pythonWorkspace/proj_CEH")
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from JASMIN import MetUM_variables as mu
from shared.paths import load_paths
from shared.utils import u_interpolate as u_int

paths = load_paths()
time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)


# =============================================================================
# SETTINGS
# =============================================================================

# Historical storm-associated P25 threshold from current analysis
DRY_THRESHOLD = -3.069373148619687

YEARS = tuple(range(1998, 2007))
MONTHS = (7, 8, 9)
STORM_HOURS = (17, 18)

# Current storm-sample filters
VMID_MIN, VMID_MAX = -50.0, 50.0
AREA_MAX_KM2 = 500000.0
GRID_KM = 4.4
PMAX_MIN = 1.0

# Environmental sampling
ENV_HOUR = 12                      # 10-12 UTC mean
HALF_BOX = 22                      # ±22 px = SM_2deg definition
MIN_VALID = 0.5 * (44 ** 2)       # same criterion as table code
DOMAIN = (-18.0, 25.0, 5.0, 25.0)

# Increment if the sampling algorithm itself changes.
# Existing caches with an older version will then be rebuilt automatically.
CACHE_VERSION = "SMbackground_yearly_v1"

# Normally False. Completed valid yearly caches are reused automatically.
FORCE_REBUILD_DAILY = False

RAW_ROOT = Path("/home/users/cornkle/linked_CP4")

CLIM_ROOT = Path(
    "/gws/ssde/j25b/lmcs/cklein/CP_models/MCS_files/WAf/climatology"
)

ANCIL_ROOT = Path(
    "/gws/ssde/j25b/impala/shared/CP4A/ncfiles/4km/ANCILS"
)

TABLE_DIR = (
    Path(paths.cp4_2d)
    / "CP4_box_JASMIN/mean3h_v2/tables"
)

OUT_DIR = (
    Path(paths.cp4_2d)
    / "CP4_box_JASMIN/mean3h_v2/SM_background_frequency"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

YEARLY_CACHE_DIR = OUT_DIR / "yearly_daily_cache"
YEARLY_CACHE_DIR.mkdir(parents=True, exist_ok=True)

LOCATION_FILE = OUT_DIR / "unique_storm_locations.csv"

STORM_AUDIT_FILES = {
    "hist": OUT_DIR / "storm_SM_audit_hist.csv",
    "fut": OUT_DIR / "storm_SM_audit_fut.csv",
}

SUMMARY_FILE = OUT_DIR / "SM_dry_background_vs_storm_summary.csv"


def yearly_cache_file(climate, year):
    return (
        YEARLY_CACHE_DIR
        / f"background_SM_daily_{climate}_{year}.nc"
    )


# =============================================================================
# TABLE HELPERS
# =============================================================================

def read_numeric_table(path):
    df = pd.read_csv(
        path,
        index_col=0,
        low_memory=False
    ).T

    for col in df.columns:
        try:
            df[col] = pd.to_numeric(df[col])
        except (ValueError, TypeError):
            pass

    return df


def first_existing(paths_):
    for p in paths_:
        if p.exists():
            return p

    raise FileNotFoundError(
        "None of these files exists:\n"
        + "\n".join(str(p) for p in paths_)
    )


def load_storm_sample(tag):
    """
    Load and align mean/anomaly tables, then retain the eligible sample used
    for the current SM DRY/WET analysis.
    """

    mean_path = first_existing([
        TABLE_DIR
        / f"{tag}_mean_table_JASMIN_3hmeansVersion_rainMask_fullBox_ICAPE.csv",

        TABLE_DIR
        / f"{tag}_mean_table_JASMIN_3hmeansVersion_rainMask_fullBox.csv",
    ])

    anom_path = first_existing([
        TABLE_DIR
        / f"{tag}_anom_table_JASMIN_3hmeansVersion_rainMask_fullBox_ICAPE.csv",

        TABLE_DIR
        / f"{tag}_anom_table_JASMIN_3hmeansVersion_rainMask_fullBox.csv",
    ])

    mean = read_numeric_table(mean_path)
    anom = read_numeric_table(anom_path)

    # Same common-case logic as existing table analysis
    keys = ["pmax", "area", "area_1mm"]

    common = pd.merge(
        mean[keys].drop_duplicates(),
        anom[keys].drop_duplicates(),
        on=keys,
        how="inner",
    )

    mean = mean.merge(
        common,
        on=keys,
        how="inner",
    )

    anom = anom.merge(
        common,
        on=keys,
        how="inner",
    )

    # Explicitly put both tables into the same deterministic ordering
    sort_cols = [
        c for c in [
            "year", "month", "day", "hour",
            "lon", "lat",
            "pmax", "area", "area_1mm",
        ]
        if c in mean.columns and c in anom.columns
    ]

    mean = (
        mean
        .sort_values(sort_cols)
        .reset_index(drop=True)
    )

    anom = (
        anom
        .sort_values(sort_cols)
        .reset_index(drop=True)
    )

    if len(mean) != len(anom):
        raise ValueError(
            f"{tag}: mean/anomaly lengths differ after matching: "
            f"{len(mean)} vs {len(anom)}"
        )

    # Check event identity
    for col in [
        "year", "month", "day", "hour",
        "lon", "lat",
        "pmax", "area", "area_1mm",
    ]:

        if col not in mean or col not in anom:
            continue

        a = pd.to_numeric(
            mean[col],
            errors="coerce"
        ).to_numpy(float)

        b = pd.to_numeric(
            anom[col],
            errors="coerce"
        ).to_numpy(float)

        good = np.isfinite(a) & np.isfinite(b)

        if good.any():
            if not np.allclose(
                a[good],
                b[good],
                rtol=0,
                atol=1e-6,
            ):
                raise ValueError(
                    f"{tag}: mean/anomaly mismatch in {col}"
                )

    # Eligible sample
    keep = (
        mean["pmax"].ge(PMAX_MIN)
        & mean["month"].isin(MONTHS)
        & mean["hour"].isin(STORM_HOURS)
        & mean["year"].isin(YEARS)
        & mean["v_mid_fullBox"].between(
            VMID_MIN,
            VMID_MAX,
            inclusive="neither",
        )
        & (
            mean["area"] * GRID_KM**2
            < AREA_MAX_KM2
        )
        & np.isfinite(anom["SM_2deg"])
        & np.isfinite(anom["lon"])
        & np.isfinite(anom["lat"])
    )

    mean = (
        mean.loc[keep]
        .reset_index(drop=True)
    )

    anom = (
        anom.loc[keep]
        .reset_index(drop=True)
    )

    # Keep useful absolute-table diagnostics
    anom["v_mid_fullBox_abs"] = (
        mean["v_mid_fullBox"].to_numpy()
    )

    anom["area_abs"] = (
        mean["area"].to_numpy()
    )

    anom["pmax_abs"] = (
        mean["pmax"].to_numpy()
    )

    print(f"\n{tag.upper()} storm sample")
    print(
        f"  eligible storms: "
        f"{len(anom):,}"
    )
    print(
        f"  table SM P25:    "
        f"{np.nanpercentile(anom['SM_2deg'], 25):.12f}"
    )
    print(
        f"  fixed threshold: "
        f"{DRY_THRESHOLD:.12f}"
    )

    return anom


# =============================================================================
# RAW FILE HELPERS
# =============================================================================

def shift_longitude(obj):
    if (
        "longitude" in obj.coords
        and float(obj.longitude.max()) > 180
    ):
        obj = obj.assign_coords(
            longitude=obj.longitude - 360
        ).sortby("longitude")

    return obj


def file_date(path):
    """
    Parse YYYYMMDD from native CP4 filename.
    """

    m = re.search(
        r"_(\d{4})(\d{2})(\d{2})\d{4}-",
        Path(path).name,
    )

    if m is None:
        return None

    return tuple(
        map(int, m.groups())
    )


def index_raw_files(climate, var):

    raw_climate = {
        "hist": "hist",
        "fut": "future",
    }[climate]

    files = sorted(
        glob.glob(
            str(
                RAW_ROOT
                / raw_climate
                / var
                / "*.nc"
            )
        )
    )

    lookup = {}

    for path in files:
        date = file_date(path)

        if date is not None:
            lookup[date] = path

    if not lookup:
        raise FileNotFoundError(
            f"No dated files found in "
            f"{RAW_ROOT / raw_climate / var}"
        )

    return lookup


# =============================================================================
# TARGET GRID / LAND MASK / SURFACE REGRIDDING
# =============================================================================

def load_target_grid_and_weights():

    plfile = sorted(
        glob.glob(
            str(
                RAW_ROOT
                / "hist/q_pl/*.nc"
            )
        )
    )[0]

    srfcfile = sorted(
        glob.glob(
            str(
                RAW_ROOT
                / "hist/lw_out_PBLtop/*.nc"
            )
        )
    )[0]

    with xr.open_dataset(
        plfile,
        decode_times=False
    ) as ds:

        pl = shift_longitude(
            ds[mu.create_CP4_filename("q_pl")]
        )

        pl = pl.sel(
            longitude=slice(
                DOMAIN[0],
                DOMAIN[1]
            ),
            latitude=slice(
                DOMAIN[2],
                DOMAIN[3]
            ),
        )

        pl_lon = (
            pl.longitude.values.copy()
        )

        pl_lat = (
            pl.latitude.values.copy()
        )

    with xr.open_dataset(
        srfcfile,
        decode_times=False
    ) as ds:

        srfc = shift_longitude(
            ds[
                mu.create_CP4_filename(
                    "lw_out_PBLtop"
                )
            ]
        )

        srfc = srfc.sel(
            longitude=slice(
                DOMAIN[0],
                DOMAIN[1]
            ),
            latitude=slice(
                DOMAIN[2],
                DOMAIN[3]
            ),
        )

        srfc_lon = (
            srfc.longitude.values.copy()
        )

        srfc_lat = (
            srfc.latitude.values.copy()
        )

    inds, weights, shape = (
        u_int.interpolation_weights(
            srfc_lon,
            srfc_lat,
            pl_lon,
            pl_lat,
        )
    )

    return (
        pl_lat,
        pl_lon,
        inds,
        weights,
        shape,
    )


(
    PL_LAT,
    PL_LON,
    INTERP_INDS,
    INTERP_WEIGHTS,
    INTERP_SHAPE,
) = load_target_grid_and_weights()


def load_sea_mask():

    f = sorted(
        glob.glob(
            str(
                ANCIL_ROOT
                / "landseamask*.nc"
            )
        )
    )[0]

    with xr.open_dataset(
        f,
        decode_times=False
    ) as ds:

        lsm = ds["lsm"]

        if "rlon" in lsm.coords:

            lsm = lsm.assign_coords(
                rlon=lsm.rlon.values - 360
            )

            lsm = lsm.sel(
                rlon=slice(
                    DOMAIN[0],
                    DOMAIN[1]
                ),
                rlat=slice(
                    DOMAIN[2],
                    DOMAIN[3]
                ),
            )

        arr = np.asarray(
            lsm.squeeze().values
        )

    target_shape = (
        len(PL_LAT),
        len(PL_LON)
    )

    if arr.shape != target_shape:
        raise ValueError(
            f"Land-mask shape {arr.shape} "
            f"does not match target grid "
            f"{target_shape}"
        )

    return arr == 0


SEA_MASK = load_sea_mask()


def regrid_surface(da):

    vals = u_int.interpolate_data(
        np.asarray(da.values),
        INTERP_INDS,
        INTERP_WEIGHTS,
        INTERP_SHAPE,
    )

    vals = np.asarray(
        vals,
        dtype=float
    )

    vals[SEA_MASK] = np.nan

    return vals


# =============================================================================
# DAILY 10-12 UTC RAW / CLIMATOLOGY FIELDS
# =============================================================================

def raw_noon_field(path, var):
    """
    Reproduce environmental-time sampling.

    SM:
        10-12 UTC mean

    lsRain:
        10-12 UTC accumulation corresponding to lsRain_noon
    """

    vn = mu.create_CP4_filename(var)

    with xr.open_dataset(
        path,
        decode_times=time_coder,
    ) as ds:

        da = shift_longitude(
            ds[vn]
        )

        da = da.sel(
            longitude=slice(
                DOMAIN[0],
                DOMAIN[1]
            ),
            latitude=slice(
                DOMAIN[2],
                DOMAIN[3]
            ),
        )

        if var == "SM":

            try:
                da = da.sel(
                    depth=0.05
                )
            except Exception:
                pass

            da = da.where(
                da < 500
            )

        hh = da.time.dt.hour

        da = da.where(
            (
                hh >= ENV_HOUR - 2
            )
            &
            (
                hh <= ENV_HOUR
            ),
            drop=True,
        )

        if (
            da.sizes.get("time", 0)
            == 0
        ):
            raise ValueError(
                f"No 10-12 UTC data in "
                f"{path}"
            )

        if var == "SM":

            out = da.mean(
                "time",
                skipna=True
            )

        elif var == "lsRain":

            out = (
                da.sum(
                    "time",
                    skipna=True
                )
                * 3600.0
            )

        else:
            raise ValueError(var)

        return out.load()


def climatology_field(
    climate,
    month,
    day,
    var,
):

    ctag = (
        "historical"
        if climate == "hist"
        else "future"
    )

    if var == "SM":

        fname = (
            CLIM_ROOT
            / ctag
            / (
                f"SM_{month:02d}-{day:02d}_"
                f"12_{ctag}_3HourMean.nc"
            )
        )

    elif var == "lsRain":

        fname = (
            CLIM_ROOT
            / ctag
            / (
                f"lsRain_{month:02d}-{day:02d}_"
                f"12_{ctag}.nc"
            )
        )

    else:
        raise ValueError(var)

    if not fname.exists():
        raise FileNotFoundError(
            fname
        )

    vn = mu.create_CP4_filename(
        var
    )

    with xr.open_dataset(
        fname
    ) as ds:

        da = shift_longitude(
            ds[vn]
        )

        da = da.sel(
            longitude=slice(
                DOMAIN[0],
                DOMAIN[1]
            ),
            latitude=slice(
                DOMAIN[2],
                DOMAIN[3]
            ),
        ).squeeze()

        if var == "SM":

            try:
                da = da.sel(
                    depth=0.05
                )
            except Exception:
                pass

            da = da.where(
                da < 500
            )

        elif var == "lsRain":

            da = (
                da * 3600.0
            )

        return da.load()


def make_daily_sm_anomaly(
    sm_path,
    rain_path,
    clim_sm,
    clim_rain,
):
    """
    Build daily SM anomaly on the environmental target grid and apply
    the same noon-rain validity mask used for storm SM_2deg.
    """

    raw_sm = raw_noon_field(
        sm_path,
        "SM"
    )

    raw_rain = raw_noon_field(
        rain_path,
        "lsRain"
    )

    # Fail explicitly if the raw/climatology grids differ unexpectedly
    raw_sm, clim_sm_here = xr.align(
        raw_sm,
        clim_sm,
        join="exact",
    )

    raw_rain, clim_rain_here = xr.align(
        raw_rain,
        clim_rain,
        join="exact",
    )

    # Anomaly before interpolation, as in the original sampling workflow
    sm_anom = regrid_surface(
        raw_sm
        - clim_sm_here
    )

    rain_anom = regrid_surface(
        raw_rain
        - clim_rain_here
    )

    # Same validity criterion used for SM_2deg
    valid = (
        np.isfinite(rain_anom)
        & (rain_anom <= 0)
    )

    sm_anom[~valid] = np.nan

    return sm_anom


# =============================================================================
# MAP STORMS TO UNIQUE TARGET-GRID LOCATIONS
# =============================================================================

def nearest_indices(
    values,
    grid,
):
    return np.asarray(
        [
            int(
                np.argmin(
                    np.abs(
                        grid - val
                    )
                )
            )
            for val in values
        ],
        dtype=int,
    )


def add_grid_cells(df):

    out = df.copy()

    out["iy"] = nearest_indices(
        out["lat"].to_numpy(float),
        PL_LAT,
    )

    out["ix"] = nearest_indices(
        out["lon"].to_numpy(float),
        PL_LON,
    )

    out["grid_lat"] = (
        PL_LAT[out["iy"]]
    )

    out["grid_lon"] = (
        PL_LON[out["ix"]]
    )

    out["date_id"] = (
        out["year"].astype(int) * 10000
        + out["month"].astype(int) * 100
        + out["day"].astype(int)
    )

    return out


def make_location_table(
    hist,
    fut,
):

    h = (
        hist[
            ["iy", "ix"]
        ]
        .drop_duplicates()
    )

    f = (
        fut[
            ["iy", "ix"]
        ]
        .drop_duplicates()
    )

    hset = set(
        map(
            tuple,
            h.to_numpy()
        )
    )

    fset = set(
        map(
            tuple,
            f.to_numpy()
        )
    )

    union = sorted(
        hset | fset
    )

    locs = pd.DataFrame(
        union,
        columns=[
            "iy",
            "ix",
        ],
    )

    locs.insert(
        0,
        "location_id",
        np.arange(
            len(locs),
            dtype=int,
        ),
    )

    locs["lat"] = (
        PL_LAT[locs["iy"]]
    )

    locs["lon"] = (
        PL_LON[locs["ix"]]
    )

    locs["hist_location"] = [
        (iy, ix) in hset
        for iy, ix in zip(
            locs["iy"],
            locs["ix"],
        )
    ]

    locs["fut_location"] = [
        (iy, ix) in fset
        for iy, ix in zip(
            locs["iy"],
            locs["ix"],
        )
    ]

    return locs


def add_location_ids(
    df,
    locs,
):

    lookup = {
        (
            int(r.iy),
            int(r.ix),
        ):
        int(r.location_id)

        for r
        in locs.itertuples()
    }

    out = df.copy()

    out["location_id"] = [
        lookup[
            (
                int(iy),
                int(ix),
            )
        ]
        for iy, ix
        in zip(
            out["iy"],
            out["ix"],
        )
    ]

    return out


def location_signature(locs):
    """
    Fingerprint of the exact set/order of sampled locations.
    """

    arr = (
        locs[
            ["iy", "ix"]
        ]
        .to_numpy(
            dtype=np.int32
        )
    )

    return hashlib.sha1(
        arr.tobytes()
    ).hexdigest()


# =============================================================================
# FAST ±22-PIXEL BOX SAMPLING
# =============================================================================

def integral_image(a):

    return (
        np.pad(
            a,
            (
                (1, 0),
                (1, 0),
            ),
        )
        .cumsum(axis=0)
        .cumsum(axis=1)
    )


def rectangle_sum(
    ii,
    y0,
    y1,
    x0,
    x1,
):

    return (
        ii[y1, x1]
        - ii[y0, x1]
        - ii[y1, x0]
        + ii[y0, x0]
    )


def sample_sm2deg(
    field,
    iy,
    ix,
):
    """
    Sample ±22-gridpoint SM box around every requested location.
    """

    iy = np.asarray(
        iy,
        dtype=int,
    )

    ix = np.asarray(
        ix,
        dtype=int,
    )

    ny, nx = field.shape

    y0 = iy - HALF_BOX
    y1 = iy + HALF_BOX + 1

    x0 = ix - HALF_BOX
    x1 = ix + HALF_BOX + 1

    inside = (
        (y0 >= 0)
        & (x0 >= 0)
        & (y1 <= ny)
        & (x1 <= nx)
    )

    result = np.full(
        len(iy),
        np.nan,
        dtype=float,
    )

    if not inside.any():
        return result

    finite = np.isfinite(
        field
    )

    sums = integral_image(
        np.where(
            finite,
            field,
            0.0,
        )
    )

    counts = integral_image(
        finite.astype(
            np.int32
        )
    )

    ii = np.where(
        inside
    )[0]

    ss = rectangle_sum(
        sums,
        y0[ii],
        y1[ii],
        x0[ii],
        x1[ii],
    )

    nn = rectangle_sum(
        counts,
        y0[ii],
        y1[ii],
        x0[ii],
        x1[ii],
    )

    good = (
        nn > MIN_VALID
    )

    result[
        ii[good]
    ] = (
        ss[good]
        / nn[good]
    )

    return result


# =============================================================================
# YEARLY CACHE VALIDATION
# =============================================================================

def yearly_cache_is_valid(
    path,
    climate,
    year,
    loc_signature,
    n_locations,
):
    """
    Only reuse a cache if its key metadata agree with this run.
    """

    if not path.exists():
        return False

    try:

        with xr.open_dataset(
            path
        ) as ds:

            if (
                ds.attrs.get(
                    "cache_version"
                )
                != CACHE_VERSION
            ):
                return False

            if (
                ds.attrs.get(
                    "climate"
                )
                != climate
            ):
                return False

            if (
                int(
                    ds.attrs.get(
                        "year",
                        -999,
                    )
                )
                != int(year)
            ):
                return False

            if (
                ds.attrs.get(
                    "location_signature"
                )
                != loc_signature
            ):
                return False

            if (
                ds.sizes.get(
                    "location"
                )
                != n_locations
            ):
                return False

            if (
                ds.sizes.get(
                    "date",
                    0
                )
                == 0
            ):
                return False

            if (
                "SM_2deg_anom"
                not in ds
            ):
                return False

            years_here = (
                np.asarray(
                    ds["year"].values
                )
            )

            if not np.all(
                years_here == year
            ):
                return False

        return True

    except Exception as exc:

        print(
            f"Cache validation failed for "
            f"{path}: {exc}"
        )

        return False


# =============================================================================
# BUILD ONE YEAR OF DAILY DATA
# =============================================================================

def build_yearly_cache(
    climate,
    year,
    locs,
    outfile,
    loc_signature,
    sm_files,
    rain_files,
):

    print(
        f"\nBuilding "
        f"{climate.upper()} {year} cache"
    )

    print(
        f"Output: {outfile}"
    )

    dates = sorted(
        d
        for d
        in (
            set(sm_files)
            & set(rain_files)
        )
        if (
            d[0] == year
            and d[1] in MONTHS
        )
    )

    if not dates:
        raise ValueError(
            f"No usable {climate} "
            f"dates for {year}"
        )

    values = np.full(
        (
            len(dates),
            len(locs),
        ),
        np.nan,
        dtype=np.float32,
    )

    loc_iy = (
        locs["iy"]
        .to_numpy(int)
    )

    loc_ix = (
        locs["ix"]
        .to_numpy(int)
    )

    # For one model year there is only one entry for each calendar day,
    # but this structure keeps climatology loading explicit and simple.
    calendar_groups = {}

    for pos, (
        yy,
        month,
        day,
    ) in enumerate(dates):

        calendar_groups.setdefault(
            (
                month,
                day,
            ),
            [],
        ).append(
            (
                pos,
                yy,
                month,
                day,
            )
        )

    n_calendar = len(
        calendar_groups
    )

    for cc, (
        (
            month,
            day,
        ),
        entries,
    ) in enumerate(
        sorted(
            calendar_groups.items()
        ),
        start=1,
    ):

        print(
            f"{climate.upper()} "
            f"{year} "
            f"{month:02d}-{day:02d} "
            f"({cc}/{n_calendar})"
        )

        clim_sm = climatology_field(
            climate,
            month,
            day,
            "SM",
        )

        clim_rain = climatology_field(
            climate,
            month,
            day,
            "lsRain",
        )

        for (
            pos,
            yy,
            month,
            day,
        ) in entries:

            key = (
                yy,
                month,
                day,
            )

            field = make_daily_sm_anomaly(
                sm_files[key],
                rain_files[key],
                clim_sm,
                clim_rain,
            )

            values[
                pos,
                :
            ] = sample_sm2deg(
                field,
                loc_iy,
                loc_ix,
            ).astype(
                np.float32
            )

    date_ids = np.asarray(
        [
            yy * 10000
            + month * 100
            + day

            for (
                yy,
                month,
                day,
            )
            in dates
        ],
        dtype=np.int32,
    )

    ds = xr.Dataset(
        data_vars={
            "SM_2deg_anom": (
                (
                    "date",
                    "location",
                ),
                values,
            ),
        },

        coords={
            "date": np.arange(
                len(dates),
                dtype=np.int32,
            ),

            "location": (
                locs[
                    "location_id"
                ]
                .to_numpy(
                    np.int32
                )
            ),

            "date_id": (
                "date",
                date_ids,
            ),

            "year": (
                "date",
                np.asarray(
                    [
                        d[0]
                        for d
                        in dates
                    ],
                    dtype=np.int16,
                ),
            ),

            "month": (
                "date",
                np.asarray(
                    [
                        d[1]
                        for d
                        in dates
                    ],
                    dtype=np.int8,
                ),
            ),

            "day": (
                "date",
                np.asarray(
                    [
                        d[2]
                        for d
                        in dates
                    ],
                    dtype=np.int8,
                ),
            ),

            "iy": (
                "location",
                locs[
                    "iy"
                ].to_numpy(
                    np.int32
                ),
            ),

            "ix": (
                "location",
                locs[
                    "ix"
                ].to_numpy(
                    np.int32
                ),
            ),

            "lat": (
                "location",
                locs[
                    "lat"
                ].to_numpy(
                    float
                ),
            ),

            "lon": (
                "location",
                locs[
                    "lon"
                ].to_numpy(
                    float
                ),
            ),

            "hist_location": (
                "location",
                locs[
                    "hist_location"
                ].to_numpy(
                    np.int8
                ),
            ),

            "fut_location": (
                "location",
                locs[
                    "fut_location"
                ].to_numpy(
                    np.int8
                ),
            ),
        },

        attrs={
            "climate": climate,
            "year": int(year),

            "cache_version": (
                CACHE_VERSION
            ),

            "location_signature": (
                loc_signature
            ),

            "description": (
                "Daily 10-12 UTC "
                "SM_2deg anomalies "
                "sampled once per "
                "unique eligible MCS "
                "location."
            ),

            "threshold_not_applied": (
                "Raw SM_2deg anomalies "
                "are stored so DRY "
                "thresholds can be "
                "changed without "
                "rebuilding."
            ),

            "rain_mask": (
                "SM pixels retained "
                "where noon rainfall "
                "anomaly <= 0."
            ),

            "half_box_pixels": (
                HALF_BOX
            ),

            "environment_hour": (
                ENV_HOUR
            ),
        },
    )

    encoding = {
        "SM_2deg_anom": {
            "zlib": True,
            "complevel": 4,
            "dtype": "float32",
        }
    }

    # Atomic-ish write:
    # interrupted write remains *.tmp rather than being accepted as cache
    tmpfile = Path(
        str(outfile)
        + ".tmp"
    )

    if tmpfile.exists():
        tmpfile.unlink()

    ds.to_netcdf(
        tmpfile,
        encoding=encoding,
    )

    ds.close()

    os.replace(
        tmpfile,
        outfile,
    )

    print(
        f"Saved completed cache: "
        f"{outfile}"
    )


# =============================================================================
# BUILD / REUSE ALL YEARLY CACHES
# =============================================================================

def ensure_yearly_caches(
    climate,
    locs,
    loc_signature,
):

    print(
        "\n"
        + "=" * 70
    )

    print(
        f"{climate.upper()} "
        f"YEARLY DAILY CACHE"
    )

    print(
        "=" * 70
    )

    sm_files = index_raw_files(
        climate,
        "SM",
    )

    rain_files = index_raw_files(
        climate,
        "lsRain",
    )

    usable = []

    for year in YEARS:

        outfile = yearly_cache_file(
            climate,
            year,
        )

        valid = yearly_cache_is_valid(
            outfile,
            climate,
            year,
            loc_signature,
            len(locs),
        )

        if (
            FORCE_REBUILD_DAILY
            or not valid
        ):

            build_yearly_cache(
                climate,
                year,
                locs,
                outfile,
                loc_signature,
                sm_files,
                rain_files,
            )

        else:

            print(
                f"Using existing "
                f"{climate.upper()} "
                f"{year} cache: "
                f"{outfile.name}"
            )

        # Validate once more after creation/reuse
        if not yearly_cache_is_valid(
            outfile,
            climate,
            year,
            loc_signature,
            len(locs),
        ):
            raise RuntimeError(
                f"Cache is still invalid "
                f"after processing: "
                f"{outfile}"
            )

        usable.append(
            outfile
        )

    return usable


# =============================================================================
# STORM VALUES FROM EXACT SAME DAILY CACHE
# =============================================================================

def storm_values_from_yearly_caches(
    storms,
    climate,
):

    result = storms.copy()

    sampled = np.full(
        len(result),
        np.nan,
        dtype=float,
    )

    for year in YEARS:

        rows = (
            result["year"].astype(int)
            == year
        )

        if not rows.any():
            continue

        cache_file = yearly_cache_file(
            climate,
            year,
        )

        with xr.open_dataset(
            cache_file
        ) as ds:

            date_ids = (
                ds[
                    "date_id"
                ]
                .values
                .astype(int)
            )

            date_lookup = {
                int(date_id): i
                for i, date_id
                in enumerate(date_ids)
            }

            # Small cache array:
            # ~92 days x N unique locations
            sm = (
                ds[
                    "SM_2deg_anom"
                ]
                .values
            )

        inds = np.where(
            rows.to_numpy()
        )[0]

        for i in inds:

            row = result.iloc[i]

            date_pos = (
                date_lookup.get(
                    int(
                        row["date_id"]
                    )
                )
            )

            if date_pos is None:
                continue

            sampled[i] = sm[
                date_pos,
                int(
                    row[
                        "location_id"
                    ]
                ),
            ]

    result[
        "SM_2deg_resampled"
    ] = sampled

    return result


# =============================================================================
# STATISTICS
# =============================================================================

def fraction_below(
    values,
    threshold,
):

    values = np.asarray(
        values,
        dtype=float,
    )

    good = np.isfinite(
        values
    )

    if not good.any():
        return (
            np.nan,
            0,
            0,
        )

    n = int(
        good.sum()
    )

    k = int(
        np.sum(
            values[good]
            < threshold
        )
    )

    return (
        k / n,
        k,
        n,
    )


def background_fraction(
    climate,
    location_mask,
):
    """
    Pool unique location-days across all yearly caches.

    Also returns the mean location-wise DRY fraction as a sensitivity.
    """

    location_mask = np.asarray(
        location_mask,
        dtype=bool,
    )

    total_valid = 0
    total_dry = 0

    loc_valid = np.zeros(
        location_mask.sum(),
        dtype=np.int64,
    )

    loc_dry = np.zeros(
        location_mask.sum(),
        dtype=np.int64,
    )

    for year in YEARS:

        cache_file = yearly_cache_file(
            climate,
            year,
        )

        with xr.open_dataset(
            cache_file
        ) as ds:

            arr = (
                ds[
                    "SM_2deg_anom"
                ]
                .values[
                    :,
                    location_mask,
                ]
            )

        good = np.isfinite(
            arr
        )

        dry = (
            good
            & (
                arr
                < DRY_THRESHOLD
            )
        )

        total_valid += int(
            good.sum()
        )

        total_dry += int(
            dry.sum()
        )

        loc_valid += (
            good.sum(
                axis=0
            )
        )

        loc_dry += (
            dry.sum(
                axis=0
            )
        )

    pooled = (
        total_dry
        / total_valid
    )

    loc_frac = np.divide(
        loc_dry,
        loc_valid,

        out=np.full(
            loc_valid.shape,
            np.nan,
            dtype=float,
        ),

        where=(
            loc_valid > 0
        ),
    )

    location_equal = (
        np.nanmean(
            loc_frac
        )
    )

    return (
        pooled,
        total_dry,
        total_valid,
        location_equal,
    )


# =============================================================================
# MAIN
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Eligible storm samples
# -----------------------------------------------------------------------------

lahist = add_grid_cells(
    load_storm_sample(
        "hist"
    )
)

lafut = add_grid_cells(
    load_storm_sample(
        "fut"
    )
)


# -----------------------------------------------------------------------------
# 2. Unique storm-location grid cells
# -----------------------------------------------------------------------------

locs = make_location_table(
    lahist,
    lafut,
)

lahist = add_location_ids(
    lahist,
    locs,
)

lafut = add_location_ids(
    lafut,
    locs,
)

locs.to_csv(
    LOCATION_FILE,
    index=False,
)

loc_sig = location_signature(
    locs
)

print(
    "\nUnique target-grid "
    "MCS locations"
)

print(
    f"  HIST:  "
    f"{locs['hist_location'].sum():,}"
)

print(
    f"  FUT:   "
    f"{locs['fut_location'].sum():,}"
)

print(
    f"  Union: "
    f"{len(locs):,}"
)


# -----------------------------------------------------------------------------
# 3. Build or reuse YEARLY daily caches
# -----------------------------------------------------------------------------

ensure_yearly_caches(
    "hist",
    locs,
    loc_sig,
)

ensure_yearly_caches(
    "fut",
    locs,
    loc_sig,
)


# -----------------------------------------------------------------------------
# 4. Independently sampled storm values
#
# IMPORTANT:
# These are not separately calculated using a second storm sampler.
# They are looked up directly from the same daily location x date caches
# used for the background statistics.
# -----------------------------------------------------------------------------

lahist_audit = (
    storm_values_from_yearly_caches(
        lahist,
        "hist",
    )
)

lafut_audit = (
    storm_values_from_yearly_caches(
        lafut,
        "fut",
    )
)

lahist_audit.to_csv(
    STORM_AUDIT_FILES[
        "hist"
    ],
    index=False,
)

lafut_audit.to_csv(
    STORM_AUDIT_FILES[
        "fut"
    ],
    index=False,
)


# -----------------------------------------------------------------------------
# 5. Existing table cross-check
# -----------------------------------------------------------------------------

(
    table_h,
    table_h_n,
    table_h_tot,
) = fraction_below(
    lahist_audit[
        "SM_2deg"
    ],
    DRY_THRESHOLD,
)

(
    table_f,
    table_f_n,
    table_f_tot,
) = fraction_below(
    lafut_audit[
        "SM_2deg"
    ],
    DRY_THRESHOLD,
)


# -----------------------------------------------------------------------------
# 6. Self-consistent storm fractions from daily cache
# -----------------------------------------------------------------------------

(
    storm_h,
    storm_h_n,
    storm_h_tot,
) = fraction_below(
    lahist_audit[
        "SM_2deg_resampled"
    ],
    DRY_THRESHOLD,
)

(
    storm_f,
    storm_f_n,
    storm_f_tot,
) = fraction_below(
    lafut_audit[
        "SM_2deg_resampled"
    ],
    DRY_THRESHOLD,
)


# -----------------------------------------------------------------------------
# 7. Background occurrence
#
# Primary:
# each climate's own unique MCS-location set
# -----------------------------------------------------------------------------

hist_location_mask = (
    locs[
        "hist_location"
    ]
    .to_numpy(bool)
)

fut_location_mask = (
    locs[
        "fut_location"
    ]
    .to_numpy(bool)
)

all_location_mask = np.ones(
    len(locs),
    dtype=bool,
)


(
    bg_h,
    bg_h_n,
    bg_h_tot,
    bg_h_equal,
) = background_fraction(
    "hist",
    hist_location_mask,
)


(
    bg_f,
    bg_f_n,
    bg_f_tot,
    bg_f_equal,
) = background_fraction(
    "fut",
    fut_location_mask,
)


# Same spatial sample in both climates:
# union of locations occupied by storms in either climate
(
    bg_h_union,
    _,
    _,
    _,
) = background_fraction(
    "hist",
    all_location_mask,
)


(
    bg_f_union,
    _,
    _,
    _,
) = background_fraction(
    "fut",
    all_location_mask,
)


# -----------------------------------------------------------------------------
# 8. Relative climate changes / enrichment
# -----------------------------------------------------------------------------

storm_ratio = (
    storm_f
    / storm_h
)

background_ratio = (
    bg_f
    / bg_h
)

enrich_h = (
    storm_h
    / bg_h
)

enrich_f = (
    storm_f
    / bg_f
)

# Equivalent to enrich_f / enrich_h
excess_ratio = (
    storm_ratio
    / background_ratio
)


# Descriptive counterfactual:
#
# What FUT storm DRY fraction would occur if the HIST storm/background
# enrichment stayed fixed and only the background DRY frequency changed?
storm_f_background_only = (
    storm_h
    * background_ratio
)

extra_storm_pp = (
    100
    * (
        storm_f
        - storm_f_background_only
    )
)


# Absolute storm shift decomposition
storm_change = (
    storm_f
    - storm_h
)

background_component = (
    storm_f_background_only
    - storm_h
)

association_component = (
    storm_f
    - storm_f_background_only
)

if not np.isclose(
    storm_change,
    0,
):

    background_share = (
        100
        * background_component
        / storm_change
    )

    association_share = (
        100
        * association_component
        / storm_change
    )

else:

    background_share = np.nan
    association_share = np.nan


# =============================================================================
# PRINT RESULTS
# =============================================================================

print(
    "\n"
    + "=" * 72
)

print(
    f"DRY threshold: "
    f"SM_2deg anomaly < "
    f"{DRY_THRESHOLD:.6f} "
    f"kg m-2"
)

print(
    "=" * 72
)


print(
    "\nTABLE CROSS-CHECK"
)

print(
    f"HIST table storms DRY:       "
    f"{100*table_h:6.2f}%  "
    f"({table_h_n:,}/"
    f"{table_h_tot:,})"
)

print(
    f"FUT  table storms DRY:       "
    f"{100*table_f:6.2f}%  "
    f"({table_f_n:,}/"
    f"{table_f_tot:,})"
)


print(
    "\nSELF-CONSISTENT DAILY "
    "SAMPLER — STORM DAYS"
)

print(
    f"HIST storms DRY:             "
    f"{100*storm_h:6.2f}%  "
    f"({storm_h_n:,}/"
    f"{storm_h_tot:,})"
)

print(
    f"FUT  storms DRY:             "
    f"{100*storm_f:6.2f}%  "
    f"({storm_f_n:,}/"
    f"{storm_f_tot:,})"
)


print(
    "\nDifference from existing "
    "table result"
)

print(
    f"HIST: "
    f"{100*(storm_h-table_h):+.2f} "
    f"percentage points"
)

print(
    f"FUT : "
    f"{100*(storm_f-table_f):+.2f} "
    f"percentage points"
)


print(
    "\nBACKGROUND — UNIQUE MCS "
    "LOCATIONS, ONE VALUE PER DAY"
)

print(
    f"HIST DRY location-days:      "
    f"{100*bg_h:6.2f}%  "
    f"({bg_h_n:,}/"
    f"{bg_h_tot:,})"
)

print(
    f"FUT  DRY location-days:      "
    f"{100*bg_f:6.2f}%  "
    f"({bg_f_n:,}/"
    f"{bg_f_tot:,})"
)


print(
    "\nRelative HIST -> FUT increase"
)

print(
    f"Storm DRY occurrence:        "
    f"x{storm_ratio:.3f}  "
    f"({100*(storm_ratio-1):+.1f}%)"
)

print(
    f"Background DRY occurrence:   "
    f"x{background_ratio:.3f}  "
    f"({100*(background_ratio-1):+.1f}%)"
)


print(
    "\nStorm-day enrichment "
    "relative to background"
)

print(
    f"HIST storm/background:       "
    f"{enrich_h:.3f}"
)

print(
    f"FUT  storm/background:       "
    f"{enrich_f:.3f}"
)

print(
    f"FUT/HIST enrichment change:  "
    f"x{excess_ratio:.3f}"
)

print(
    f"Change beyond background:     "
    f"{100*(excess_ratio-1):+.1f}%"
)


print(
    "\nBACKGROUND-SHIFT "
    "COUNTERFACTUAL"
)

print(
    f"Expected FUT storm DRY fraction "
    f"if only background occurrence "
    f"changed: "
    f"{100*storm_f_background_only:.2f}%"
)

print(
    f"Observed FUT storm DRY fraction: "
    f"{100*storm_f:.2f}%"
)

print(
    f"Extra storm association:         "
    f"{extra_storm_pp:+.2f} "
    f"percentage points"
)


print(
    "\nFraction of total "
    "storm-day DRY shift"
)

print(
    f"Background-frequency "
    f"contribution: "
    f"{background_share:.1f}%"
)

print(
    f"Additional storm association: "
    f"{association_share:.1f}%"
)


print(
    "\nCOMMON-SPATIAL-SAMPLE "
    "SENSITIVITY"
)

print(
    f"HIST union-location "
    f"background: "
    f"{100*bg_h_union:.2f}%"
)

print(
    f"FUT  union-location "
    f"background: "
    f"{100*bg_f_union:.2f}%"
)


print(
    "\nLOCATION-EQUAL SENSITIVITY"
)

print(
    f"HIST mean location "
    f"DRY fraction: "
    f"{100*bg_h_equal:.2f}%"
)

print(
    f"FUT  mean location "
    f"DRY fraction: "
    f"{100*bg_f_equal:.2f}%"
)


# =============================================================================
# SAVE SUMMARY
# =============================================================================

summary = pd.DataFrame(
    {
        "metric": [
            "dry_threshold",

            "table_storm_hist_fraction",
            "table_storm_fut_fraction",

            "resampled_storm_hist_fraction",
            "resampled_storm_fut_fraction",

            "background_hist_fraction",
            "background_fut_fraction",

            "storm_fut_over_hist",
            "background_fut_over_hist",

            "enrichment_hist",
            "enrichment_fut",

            "enrichment_fut_over_hist",
            "enrichment_change_percent",

            "future_storm_fraction_if_background_only",
            "extra_future_storm_association_percentage_points",

            "background_share_of_storm_shift_percent",
            "association_share_of_storm_shift_percent",

            "background_hist_union_locations",
            "background_fut_union_locations",

            "background_hist_location_equal",
            "background_fut_location_equal",
        ],

        "value": [
            DRY_THRESHOLD,

            table_h,
            table_f,

            storm_h,
            storm_f,

            bg_h,
            bg_f,

            storm_ratio,
            background_ratio,

            enrich_h,
            enrich_f,

            excess_ratio,
            100 * (
                excess_ratio
                - 1
            ),

            storm_f_background_only,
            extra_storm_pp,

            background_share,
            association_share,

            bg_h_union,
            bg_f_union,

            bg_h_equal,
            bg_f_equal,
        ],
    }
)

summary.to_csv(
    SUMMARY_FILE,
    index=False,
)


# =============================================================================
# FINAL FILE SUMMARY
# =============================================================================

print(
    "\nSaved / reused yearly caches:"
)

for climate in (
    "hist",
    "fut",
):
    for year in YEARS:

        print(
            "  "
            + str(
                yearly_cache_file(
                    climate,
                    year,
                )
            )
        )


print(
    "\nSaved analysis files:"
)

print(
    f"  {LOCATION_FILE}"
)

print(
    f"  {STORM_AUDIT_FILES['hist']}"
)

print(
    f"  {STORM_AUDIT_FILES['fut']}"
)

print(
    f"  {SUMMARY_FILE}"
)