import numpy as np
import pandas as pd

from scipy.ndimage import label
from shapely.geometry import Polygon
from skimage import measure


def dictionary():
    vars = [
        "date", "month", "hour", "minute", "year", "day",
        "area", "area70", "area60", "area50",
        "minlon", "minlat", "maxlon", "maxlat",
        "clon", "clat", "tminlon", "tminlat",
        "tmin", "tmean", "tp1", "tp99",
        "stormID", "cloudMask", "tir",
        "shape_complexity", "major_axis_km", "minor_axis_km", "aspect_ratio",
    ]

    return {v: [] for v in vars}


def mask_to_projected_polygon(mask, lat, lon, transformer):
    """
    Convert a boolean storm mask to a projected shapely polygon.
    lat/lon are 1D arrays corresponding to mask rows/columns.
    """
    contours = measure.find_contours(mask.astype(float), level=0.5)
    polygons = []

    for contour in contours:
        rows = contour[:, 0]
        cols = contour[:, 1]

        contour_lat = np.interp(rows, np.arange(len(lat)), lat)
        contour_lon = np.interp(cols, np.arange(len(lon)), lon)

        x, y = transformer.transform(contour_lon, contour_lat)
        poly = Polygon(zip(x, y))

        if poly.is_valid and poly.area > 0:
            polygons.append(poly)

    if len(polygons) == 0:
        return None

    return max(polygons, key=lambda p: p.area)


def shape_complexity_index(poly):
    """
    1 - polygon area / convex hull area.
    Higher = more ragged / concave / irregular.
    """
    if poly is None or poly.is_empty:
        return np.nan

    hull_area = poly.convex_hull.area

    if hull_area == 0:
        return np.nan

    return 1 - poly.area / hull_area


def major_minor_axis_lengths(poly):
    """
    Major/minor axis estimate from minimum rotated rectangle.
    Returns metres.
    """
    if poly is None or poly.is_empty:
        return np.nan, np.nan

    rect = poly.minimum_rotated_rectangle
    coords = np.asarray(rect.exterior.coords)

    side_lengths = np.sqrt(np.sum(np.diff(coords, axis=0) ** 2, axis=1))

    major = np.max(side_lengths)
    minor = np.min(side_lengths)

    return major, minor


def mcs_define_area(array, thresh, area_grid_km2, min_area_km2=None, max_area_km2=None):
    """
    Label contiguous cold-cloud objects and filter them by true km2 area.

    Parameters
    ----------
    array : 2D np.ndarray
        Brightness temperature in degC.
    thresh : float
        Cold-cloud threshold, e.g. -40.
    area_grid_km2 : 2D np.ndarray
        Pixel area grid in km2, same shape as array.
    min_area_km2 : float
        Minimum object area in km2.
    max_area_km2 : float
        Maximum object area in km2.

    Returns
    -------
    labels : 2D np.ndarray
        Labelled object array with rejected objects set to 0.
    goodinds : np.ndarray
        Label IDs that passed filtering.
    object_areas : dict
        Mapping from label ID to true km2 area.
    """
    cold_mask = np.isfinite(array) & (array < thresh)

    labels, _ = label(cold_mask)

    goodinds = []
    object_areas = {}

    for lab in np.unique(labels):
        if lab == 0:
            continue

        storm_mask = labels == lab
        area = float(np.nansum(area_grid_km2[storm_mask]))
        object_areas[lab] = area

        if min_area_km2 is not None and area < min_area_km2:
            labels[storm_mask] = 0
            continue

        if max_area_km2 is not None and area > max_area_km2:
            labels[storm_mask] = 0
            continue

        goodinds.append(lab)

    return labels, np.asarray(goodinds), object_areas


def process_tir_image(ds, data_res, t_thresh=-40, min_mcs_size=5000, max_mcs_size=300000, area_grid_km2=None,
    transformer=None, rainvar_name=None, min_rain_pmax=None):
    """
    Identify MCS-like contiguous cold-cloud objects.

    The minimum-size filter is applied using true pixel area in km2 if
    area_grid_km2 is provided. This is preferred for lat/lon grids.

    MCS definition here:
    - contiguous cloud area colder than t_thresh, usually <= -40C
    - true area >= min_mcs_size km2
    - minimum brightness temperature <= -50C
    """
    ctt = ds["tb"].squeeze() - 273.15

    if area_grid_km2 is None:
        # Fallback to old approximate behaviour.
        area_grid_km2 = np.ones_like(ctt.values, dtype=float) * data_res**2

    labels, goodinds, object_areas = mcs_define_area(ctt.values, t_thresh, area_grid_km2=area_grid_km2,
        min_area_km2=min_mcs_size, max_area_km2=max_mcs_size)

    dic = dictionary()

    for g in goodinds:
        
        storm_mask = labels == g
        
        ###rainfall at least 1mm test - only minimally rainy storms!
        if rainvar_name is not None and min_rain_pmax is not None:
            rain = ds[rainvar_name].squeeze().where(storm_mask)
        
            if np.all(np.isnan(rain.values)):
                continue
        
            rain_max = float(np.nanmax(rain.values))
        
            if rain_max < min_rain_pmax:
                continue

        storm = ctt.where(storm_mask)

        if np.all(np.isnan(storm.values)):
            continue

        tmin = float(np.nanmin(storm.values))

        # Additional core-intensity requirement.
        if tmin > -50:
            continue

        tmin_pos = np.nanargmin(storm.values)
        tpos_2d = np.unravel_index(tmin_pos, storm.shape)

        pos = np.where(storm_mask)

        latmin = float(np.nanmin(ctt.lat.values[pos[0]]))
        latmax = float(np.nanmax(ctt.lat.values[pos[0]]))
        lonmin = float(np.nanmin(ctt.lon.values[pos[1]]))
        lonmax = float(np.nanmax(ctt.lon.values[pos[1]]))

        # Shape metrics
        if transformer is not None:
            try:
                poly = mask_to_projected_polygon(storm_mask, ctt.lat.values,
                    ctt.lon.values, transformer)

                sci = shape_complexity_index(poly)
                major_m, minor_m = major_minor_axis_lengths(poly)

                dic["shape_complexity"].append(sci)
                dic["major_axis_km"].append(major_m / 1000)
                dic["minor_axis_km"].append(minor_m / 1000)
                dic["aspect_ratio"].append(major_m / minor_m if minor_m > 0 else np.nan)

            except Exception as e:
                print("Storm shape error, setting shape metrics to NaN:", e)
                dic["shape_complexity"].append(np.nan)
                dic["major_axis_km"].append(np.nan)
                dic["minor_axis_km"].append(np.nan)
                dic["aspect_ratio"].append(np.nan)

        else:
            dic["shape_complexity"].append(np.nan)
            dic["major_axis_km"].append(np.nan)
            dic["minor_axis_km"].append(np.nan)
            dic["aspect_ratio"].append(np.nan)

        # Area metrics
        storm_values = storm.values

        area = float(np.nansum(area_grid_km2[storm_mask]))
        area70 = float(np.nansum(area_grid_km2[storm_values <= -70]))
        area60 = float(np.nansum(area_grid_km2[storm_values <= -60]))
        area50 = float(np.nansum(area_grid_km2[storm_values <= -50]))

        dic["area"].append(area)
        dic["area70"].append(area70)
        dic["area60"].append(area60)
        dic["area50"].append(area50)

        datestr = (
            str(int(ctt["time.year"].values))
            + "-"
            + str(int(ctt["time.month"].values)).zfill(2)
            + "-"
            + str(int(ctt["time.day"].values)).zfill(2)
            + "_"
            + str(int(ctt["time.hour"].values)).zfill(2)
            + ":"
            + str(int(ctt["time.minute"].values)).zfill(2)
        )

        dic["date"].append(datestr)
        dic["month"].append(int(ctt["time.month"]))
        dic["hour"].append(int(ctt["time.hour"]))
        dic["year"].append(int(ctt["time.year"]))
        dic["day"].append(int(ctt["time.day"]))
        dic["minute"].append(int(ctt["time.minute"]))

        dic["minlon"].append(lonmin)
        dic["minlat"].append(latmin)
        dic["maxlon"].append(lonmax)
        dic["maxlat"].append(latmax)
        dic["clon"].append(lonmin + (lonmax - lonmin) / 2)
        dic["clat"].append(latmin + (latmax - latmin) / 2)

        dic["tmin"].append(tmin)
        dic["tminlat"].append(float(ctt.lat[tpos_2d[0]].values))
        dic["tminlon"].append(float(ctt.lon[tpos_2d[1]].values))
        dic["tmean"].append(float(np.nanmean(storm_values)))
        dic["tp1"].append(float(np.nanpercentile(storm_values, 1)))
        dic["tp99"].append(float(np.nanpercentile(storm_values, 99)))

        dic["stormID"].append(datestr + "_" + str(g))
        dic["cloudMask"].append(storm_mask)
        dic["tir"].append(storm_values)

    return dic


def add_environment_toTable(
    tab,
    in_ds,
    data_res,
    envvar_take=[],
    tabvar_skip=[],
    rainvar_name=None,
    env_tformat="%Y-%m-%d %H:%M:%S",
    env_hour=12,
    area_grid_km2=None,
):
    """
    Add rainfall and optional environmental variables to storm table.

    Rainfall is sampled at storm time across the cold-cloud mask.
    Optional environment variables are sampled around the minimum-Tb location.
    """
    tab_tformat = "%Y-%m-%d_%H:%M"

    dic = {}
    for k in tab.keys():
        if k in tabvar_skip:
            continue
        dic[k] = tab[k]

    ds = in_ds
    envdates = pd.to_datetime(ds.time.dt.floor("min"), format=env_tformat)

    for tlat, tlon, date, mask, tir in zip(
        dic["tminlat"],
        dic["tminlon"],
        dic["date"],
        dic["cloudMask"],
        dic["tir"],
    ):

        if rainvar_name is not None:
            tabdate = pd.to_datetime(date, format=tab_tformat)
            pos = envdates == tabdate

            if np.sum(pos) == 0:
                rain = None
            else:
                rain = ds[rainvar_name].isel(time=pos).where(mask).squeeze()

            if (rainvar_name + "_mean") not in dic.keys():
                for tag in [
                    "_mean", "_max", "_p80", "_p90", "_p95", "_p99",
                    "_area1mm", "_area3mm", "_area8mm",
                ]:
                    dic[rainvar_name + tag] = []

            if rain is None or np.all(np.isnan(rain.values)):
                dic[rainvar_name + "_mean"].append(np.nan)
                dic[rainvar_name + "_max"].append(np.nan)
                dic[rainvar_name + "_p80"].append(np.nan)
                dic[rainvar_name + "_p90"].append(np.nan)
                dic[rainvar_name + "_p95"].append(np.nan)
                dic[rainvar_name + "_p99"].append(np.nan)
                dic[rainvar_name + "_area1mm"].append(np.nan)
                dic[rainvar_name + "_area3mm"].append(np.nan)
                dic[rainvar_name + "_area8mm"].append(np.nan)

            else:
                pmax_pos = np.nanargmax(rain.values)
                ppos_2d = np.unravel_index(pmax_pos, rain.shape)

                pmax_lon = float(rain.lon[ppos_2d[1]].values)
                pmax_lat = float(rain.lat[ppos_2d[0]].values)

                pmax = rain.sel(
                    lon=slice(pmax_lon - 0.075, pmax_lon + 0.075),
                    lat=slice(pmax_lat - 0.075, pmax_lat + 0.075),
                ).mean().values

                dic[rainvar_name + "_mean"].append(float(rain.mean().values))
                dic[rainvar_name + "_max"].append(float(pmax))
                dic[rainvar_name + "_p80"].append(float(rain.quantile(0.80).values))
                dic[rainvar_name + "_p90"].append(float(rain.quantile(0.90).values))
                dic[rainvar_name + "_p95"].append(float(rain.quantile(0.95).values))
                dic[rainvar_name + "_p99"].append(float(rain.quantile(0.99).values))

                if area_grid_km2 is None:
                    parea = float(np.sum(rain.values >= 1) * data_res**2)
                    strati = float(np.sum(rain.values >= 3) * data_res**2)
                    conv = float(np.sum(rain.values >= 8) * data_res**2)
                else:
                    parea = float(np.nansum(area_grid_km2[rain.values >= 1]))
                    strati = float(np.nansum(area_grid_km2[rain.values >= 3]))
                    conv = float(np.nansum(area_grid_km2[rain.values >= 8]))

                dic[rainvar_name + "_area1mm"].append(parea)
                dic[rainvar_name + "_area3mm"].append(strati)
                dic[rainvar_name + "_area8mm"].append(conv)

        if len(envvar_take) > 0:
            tabdate = pd.to_datetime(date, format=tab_tformat).replace(
                hour=env_hour,
                minute=0,
            )

            pos = envdates == tabdate

            if np.sum(pos) == 0:
                for vt in envvar_take:
                    if vt not in dic:
                        dic[vt] = []
                    dic[vt].append(np.nan)
                continue

            single = ds.isel(time=pos).sel(
                longitude=slice(tlon - 0.375, tlon + 0.375),
                latitude=slice(tlat - 0.375, tlat + 0.375),
            )

            single = single.mean()

            for vt in envvar_take:
                if vt not in dic:
                    dic[vt] = []
                dic[vt].append(single[vt].values)

    return dic