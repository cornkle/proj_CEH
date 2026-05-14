import xarray as xr
import numpy as np
import xesmf as xe


hi_file = '/prj/global_water/ERA5_global_0.7/static/era5_invariant_025deg.nc'
target_file = '/prj/global_water/ERA5_global_0.7/hourly/surface/convective_inhibition/convective_inhibition_ERA5_2000_01_02_srfc.nc'
out_file = '/prj/global_water/ERA5_global_0.7/static/era5_invariant_07deg_COCOON-study_regrid.nc'

G = 9.80665

def standardise_lon(ds, lon_name="longitude"):
    """Convert longitudes to -180..180 if needed and sort."""
    if ds[lon_name].max() > 180:
        ds = ds.assign_coords({
            lon_name: ((ds[lon_name] + 180) % 360) - 180
        })
        ds = ds.sortby(lon_name)
    return ds


def ensure_lat_ascending(ds, lat_name="latitude"):
    """xESMF is happier with ascending latitude."""
    if ds[lat_name][0] > ds[lat_name][-1]:
        ds = ds.sortby(lat_name)
    return ds


def add_bounds(ds, lat_name="latitude", lon_name="longitude"):
    """
    Add explicit 1D lat/lon bounds for conservative regridding.
    Clips latitude bounds to [-90, 90] to avoid xESMF warnings.
    """

    lat = ds[lat_name].values
    lon = ds[lon_name].values

    def bounds_from_centres(x, is_lat=False):
        mid = (x[:-1] + x[1:]) / 2

        first = x[0] - (mid[0] - x[0])
        last = x[-1] + (x[-1] - mid[-1])

        bounds = np.concatenate([[first], mid, [last]])

        if is_lat:
            bounds = np.clip(bounds, -90, 90)

        return bounds

    lat_b = bounds_from_centres(lat, is_lat=True)
    lon_b = bounds_from_centres(lon, is_lat=False)

    grid = xr.Dataset(
        {
            "lat": ds[lat_name],
            "lon": ds[lon_name],
            "lat_b": xr.DataArray(lat_b, dims=("lat_b",)),
            "lon_b": xr.DataArray(lon_b, dims=("lon_b",)),
        }
    )

    return grid


# -------------------------------------------------------------------
# 1. Open files
# -------------------------------------------------------------------

hi = xr.open_dataset(hi_file)
target = xr.open_dataset(target_file)

# If your files use lat/lon instead of latitude/longitude, rename them
rename_hi = {}
rename_target = {}

if "lat" in hi.coords and "latitude" not in hi.coords:
    rename_hi["lat"] = "latitude"
if "lon" in hi.coords and "longitude" not in hi.coords:
    rename_hi["lon"] = "longitude"

if "lat" in target.coords and "latitude" not in target.coords:
    rename_target["lat"] = "latitude"
if "lon" in target.coords and "longitude" not in target.coords:
    rename_target["lon"] = "longitude"

if rename_hi:
    hi = hi.rename(rename_hi)
if rename_target:
    target = target.rename(rename_target)


# -------------------------------------------------------------------
# 2. Standardise coordinates
# -------------------------------------------------------------------

hi = standardise_lon(hi)
target = standardise_lon(target)

hi = ensure_lat_ascending(hi)
target = ensure_lat_ascending(target)


# 3. Derive useful high-resolution variables before regridding

if "z" in hi:
    hi["elevation"] = hi["z"] / G
    hi["elevation"].attrs.update(
        units="m",
        long_name="surface elevation from ERA5 geopotential"
    )

# ERA5 name in your file: sdfor
if "sdfor" in hi:
    hi["topographic_complexity"] = hi["sdfor"]
    hi["topographic_complexity"].attrs.update(
        units=hi["sdfor"].attrs.get("units", "m"),
        long_name="standard deviation of filtered subgrid orography"
    )

if "sdfor" in hi and "slor" in hi:
    hi["topographic_complexity_sdfor_slor"] = hi["sdfor"] * hi["slor"]
    hi["topographic_complexity_sdfor_slor"].attrs.update(
        long_name="standard deviation of filtered subgrid orography multiplied by subgrid slope"
    )

# -------------------------------------------------------------------
# 4. Build source and target grids with explicit bounds
# -------------------------------------------------------------------

source_grid = add_bounds(hi)
target_grid = add_bounds(target)


# -------------------------------------------------------------------
# 5. Regrid to your existing ERA5 subset grid
# -------------------------------------------------------------------

regridder = xe.Regridder(
    source_grid,
    target_grid,
    method="conservative",
    periodic=False,
    reuse_weights=False,
)

coarse = regridder(hi)


# -------------------------------------------------------------------
# 6. Keep useful variables only
# -------------------------------------------------------------------

keep = [
    "lsm",
    "z",
    "elevation",
    "sdfor",
    "slor",
    "topographic_complexity",
    "topographic_complexity_sdfor_slor",
]

keep = [v for v in keep if v in coarse.data_vars]
coarse = coarse[keep]


# -------------------------------------------------------------------
# 7. Rename output dimensions back to ERA5-style names if needed
# -------------------------------------------------------------------

if "lat" in coarse.dims:
    coarse = coarse.rename({"lat": "latitude"})
if "lon" in coarse.dims:
    coarse = coarse.rename({"lon": "longitude"})


# -------------------------------------------------------------------
# 8. Add metadata and save to disk
# -------------------------------------------------------------------

coarse.attrs.update(
    description="ERA5 invariant fields regridded to existing 0.7 degree ERA5 subset grid",
    source_file=hi_file,
    target_grid_file=target_file,
)

for v in coarse.data_vars:
    coarse[v].encoding.clear()
    coarse[v].attrs.pop("_FillValue", None)
    coarse[v].attrs.pop("missing_value", None)

encoding = {
    v: {
        "_FillValue": np.float32(-9999.0),
        "dtype": "float32",
    }
    for v in coarse.data_vars
}

coarse.to_netcdf(
    out_file,
    encoding=encoding,
    format="NETCDF3_64BIT",
)

print(f"Saved: {out_file}")
print(coarse)