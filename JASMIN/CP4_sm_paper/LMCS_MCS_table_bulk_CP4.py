import glob
import multiprocessing
import os

import numpy as np
import pandas as pd
import xarray as xr
from metpy import calc
from metpy.units import units

from shared.utils import u_arrays as ua


TIMETAGS = ("hist", "fut")
ANOMTAGS = ("anom", "mean")
MAIN_PATH = "/gws/ssde/j25b/lmcs/cklein/CP_models/MCS_files/WAf/CP4_box_JASMIN/mean3h_v2/"
SKIP_VARS = ["lw_out_PBLtop", "lsRain"]
HOURS = (17, 18)


def dictionary(dummy):
    dic = {}

    for vn in dummy.data_vars:
        if vn in SKIP_VARS:
            continue
        dic[f"{vn}_0.25deg"] = []
        dic[f"{vn}_1deg"] = []
        dic[f"{vn}_2deg"] = []

    for vn in dummy.data_vars:
        dic[f"{vn}_Smean"] = []

    output_vars = [
        "hour", "month", "year", "day", "index", "indx", "indy", "lon", "lat", "tmin", "pmax", "pmax_native",
        "pgt30", "pgt1", "isvalid", "tgrad", "SMgrad", "SHgrad", "LHgrad", "tcwv_cl", "sh_cl", "lh_cl",
        "lsRain_noon_cl", "area", "area-70", "area-80", "area-60", "area_1mm",
        "area_5mm", "area_10mm","lsRain_noon_pix_1deg" 
    ]

    for vn in output_vars:
        dic[vn] = []

    return dic


def add_derived_fields(ds, anomtag, dx_m=4400, p_low_hpa=925, p_mid_hpa=650):
    """Add divergence to all files and thermodynamics to absolute mean files."""

    ds = ds.copy()
    u = ds["u_srfc"].values * units("m/s")
    v = ds["v_srfc"].values * units("m/s")
    spacing = dx_m * units.m
    div = calc.divergence(u, v, dx=spacing, dy=spacing).to("1/s")
    ds["div"] = (ds["u_srfc"].dims, div.magnitude)
    ds["div"].attrs.update(units="s-1", long_name="925-hPa horizontal wind divergence")

    # Theta-e, RH and the proxies are nonlinear and require absolute T and q.
    # Do not interpret anomaly T (K anomaly) and q (kg kg-1 anomaly) as a state.
    if anomtag != "mean":
        return ds

    # Assumptions:
    # t_srfc and t_mid: K
    # q_srfc and q_mid: kg kg-1
    # u_srfc and v_srfc: m s-1
    t_low = ds["t_srfc"].values * units.K
    t_mid = ds["t_mid"].values * units.K
    q_low = ds["q_srfc"].values * units("kg/kg")
    q_mid = ds["q_mid"].values * units("kg/kg")
    p_low = p_low_hpa * units.hPa
    p_mid = p_mid_hpa * units.hPa

    # Relative humidity as a fraction; multiply by 100 below for percentage.
    rh_low = calc.relative_humidity_from_specific_humidity(p_low, t_low, q_low)
    rh_mid = calc.relative_humidity_from_specific_humidity(p_mid, t_mid, q_mid)

    # Equivalent potential temperatures.
    td_low = calc.dewpoint_from_specific_humidity(p_low, t_low, q_low)
    td_mid = calc.dewpoint_from_specific_humidity(p_mid, t_mid, q_mid)

    thetae_low = calc.equivalent_potential_temperature(p_low, t_low, td_low)
    thetae_mid = calc.equivalent_potential_temperature(p_mid, t_mid, td_mid)
    thetaes_mid = calc.saturation_equivalent_potential_temperature(p_mid, t_mid)

    # Positive values indicate greater potential for a low-level parcel to
    # remain buoyant after reaching the mid-level environment.
    cape_proxy = (thetae_low - thetaes_mid).to("K")

    # This is not true CIN. It measures how far the mid-level environment is
    # from saturation in moist-entropy space. Larger values imply a drier,
    # potentially more entraining/suppressive mid-level environment.
    cin_proxy = (thetaes_mid - thetae_mid).to("K")

    ds["rh_srfc"] = (ds["q_srfc"].dims, rh_low.to("dimensionless").magnitude * 100)
    ds["rh_mid"] = (ds["q_mid"].dims, rh_mid.to("dimensionless").magnitude * 100)
    ds["cape_proxy"] = (ds["t_srfc"].dims, cape_proxy.magnitude)
    ds["cin_proxy"] = (ds["t_mid"].dims, cin_proxy.magnitude)
    ds["thetae_srfc"] = (ds["t_srfc"].dims, thetae_low.to("K").magnitude)
    ds["thetae_mid"] = (ds["t_mid"].dims, thetae_mid.to("K").magnitude)
    ds["thetaes_mid"] = (ds["t_mid"].dims, thetaes_mid.to("K").magnitude)

    ds["rh_srfc"].attrs.update(units="%", long_name="925-hPa relative humidity")
    ds["rh_mid"].attrs.update(units="%", long_name="650-hPa relative humidity")
    ds["cape_proxy"].attrs.update(
        units="K",
        long_name="two-level instability proxy: theta-e 925 minus saturated theta-e 650"
    )
    ds["cin_proxy"].attrs.update(
        units="K",
        long_name="650-hPa saturation deficit: saturated theta-e minus theta-e"
    )

    return ds


def file_loop(args):
  
    f, timetag, anomtag = args
    print(f"Doing file: {f}")

    ds = xr.open_dataset(f)
    ds = add_derived_fields(ds, anomtag)
    out = dictionary(ds)

    bname = os.path.basename(f)
    parts = bname.split("_")
    index = int(parts[2])
    ilon = float(parts[4].strip("[]"))
    ilat = float(parts[5].strip("[]"))

    outt = ds["lw_out_PBLtop"].values
    outp = ds["lsRain"].values
    outp_noon = ds["lsRain_noon"].values

    # if np.sum(np.isnan(outp_noon[54:62, 54:62])) > 0:
    #     print("Noon rainfall, continue")
    #     ds.close()
    #     return None

    if np.nanmax(outp_noon[54:62, 54:62]) > 0.25:
        print("Noon rainfall, continue")
        ds.close()
        return None
    print("RAINFALL FREE!!!!!!")

    u = units.Quantity(ds["u_srfc"].values, "m/s")
    v = units.Quantity(ds["v_srfc"].values, "m/s")
    dx = units.Quantity(4400, "m")
    div = calc.divergence(u, v, dx=dx, dy=dx)

    out["lon"] = ds.minlon
    out["lat"] = ds.minlat
    out["hour"] = ds["time.hour"].item()
    out["month"] = ds["time.month"].item()
    out["year"] = ds["time.year"].item()
    out["day"] = ds["time.day"].item()

    # Unique identification of MCS files from the source table
    out["index"] = index
    out["indx"] = ilon
    out["indy"] = ilat

    t_thresh = -50
    mask = outt <= t_thresh
    varmask = np.isfinite(outp_noon) & (outp_noon <= 0)

    if np.sum(mask) <= 3:
        ds.close()
        return None

    rainfield = outp[mask]

    out["area"] = np.sum(outt <= t_thresh)
    out["area-60"] = np.sum(outt <= -60)
    out["area-70"] = np.sum(outt <= -70)
    out["area-80"] = np.sum(outt <= -80)
    out["area_1mm"] = np.sum(outp >= 1)
    out["area_5mm"] = np.sum(outp >= 5)
    out["area_10mm"] = np.sum(outp >= 10)

    if not np.any(np.isfinite(outp)):
        print(f"No finite storm-time rainfall: {f}")
        ds.close()
        return None

    maxpos = np.unravel_index(np.nanargmax(outp), outp.shape)
    minpos = np.unravel_index(np.nanargmin(outt), outt.shape)

    out["tmin"] = np.nanmin(outt)
    out["pmax"] = np.nanmean(ua.cut_kernel(outp, maxpos[1], maxpos[0], 1))
    out["pmax_native"] = np.nanmax(outp[mask])

    if anomtag == "anom":
        basefile = os.path.basename(f)
        dpath = os.path.join(MAIN_PATH, f"mean_{timetag}")

        try:
            with xr.open_dataset(os.path.join(dpath, basefile)) as dcl:
                out["tcwv_cl"] = np.nanmean(ua.cut_kernel(dcl["tcwv"].values, minpos[1], minpos[0], 11))
                out["sh_cl"] = np.nanmean(ua.cut_kernel(dcl["sh"].values, minpos[1], minpos[0], 11))
                out["lh_cl"] = np.nanmean(ua.cut_kernel(dcl["lh"].values, minpos[1], minpos[0], 11))
                out["lsRain_noon_cl"] = np.nanmean(
                    ua.cut_kernel(dcl["lsRain_noon"].values, minpos[1], minpos[0], 11)
                )
                
        except Exception as exc:
            print(f"Could not open matching mean file for {basefile}: {exc}")
            ds.close()
            return None

    boxes = ((11, "_1deg"), (22, "_2deg"), (3, "_0.25deg"), (mask, "_Smean"))

    for dist, tag in boxes:
        for vn in ds.data_vars:

            if isinstance(dist, int) & (vn not in SKIP_VARS):
                cnt = ds[vn].where(varmask).sel(latitude=slice(dist*-1, dist), longitude=slice(dist*-1, dist)).count(['longitude', 'latitude'])
                if cnt > 0.5 * ((dist*2)**2):
                    data = float(ds[vn].where(varmask).sel(latitude=slice(dist*-1, dist), longitude=slice(dist*-1, dist)).mean(['longitude', 'latitude']).values)  # box of 1deg or 0.25deg
                else:
                    data = np.nan
            else:
                data = float(ds[vn].where((mask) & (varmask)).mean(['longitude', 'latitude']).values)  # storm mask

            out[f"{vn}{tag}"] = float(data)
          
            if '1deg' in tag:
                out["lsRain_noon_pix"+tag] = np.nansum(ds["lsRain_noon"].sel(latitude=slice(dist*-1, dist), longitude=slice(dist*-1, dist)) > 1)

    

    tgrad_lat = ds.isel(longitude=slice(minpos[1] - 11, minpos[1] + 11)).mean("longitude").squeeze()
    tgrad = tgrad_lat.polyfit(dim="latitude", deg=1)

    out["tgrad"] = float(tgrad["t2_polyfit_coefficients"].sel(degree=1))
    out["SMgrad"] = float(tgrad["SM_polyfit_coefficients"].sel(degree=1))
    out["SHgrad"] = float(tgrad["sh_polyfit_coefficients"].sel(degree=1))
    out["LHgrad"] = float(tgrad["lh_polyfit_coefficients"].sel(degree=1))

    out["pgt30"] = np.sum(rainfield > 30)
    out["isvalid"] = np.sum(mask)
    out["pgt1"] = np.sum(rainfield > 1)

    ds.close()
    return out


def run_case(timetag, anomtag, processes=4):
    input_dir = os.path.join(MAIN_PATH, f"{anomtag}_{timetag}")
    other_dir = os.path.join(MAIN_PATH, f"{'mean' if anomtag == 'anom' else 'anom'}_{timetag}")

    # Only use 17/18 UTC cases whose complete filename exists in both folders.
    files_here = {os.path.basename(f): f for f in glob.glob(os.path.join(input_dir, "*.nc"))
                  if int(os.path.basename(f).split("_")[1].split(":")[0]) in HOURS}
    names_other = {os.path.basename(f) for f in glob.glob(os.path.join(other_dir, "*.nc"))
                   if int(os.path.basename(f).split("_")[1].split(":")[0]) in HOURS}
    matched_names = sorted(set(files_here) & names_other)
    files = [files_here[name] for name in matched_names]

    print(f"\nProcessing timetag={timetag}, anomtag={anomtag}")
    print(f"Nb matched 17/18 UTC files: {len(files)}")
    print(f"Excluded because counterpart is missing: {len(set(files_here) - names_other)}")

    if not files:
        print(f"No NetCDF files found in {input_dir}; skipping.")
        return

    with xr.open_dataset(files[0]) as dummy:
        dummy = add_derived_fields(dummy, anomtag)
        mdic = dictionary(dummy)

    tasks = [(f, timetag, anomtag) for f in files]

    with multiprocessing.Pool(processes=processes) as pool:
        res = pool.map(file_loop, tasks)

    for result in res:
        if result is None:
            continue
        for key in mdic:
            mdic[key].append(result[key])

    df = pd.DataFrame.from_dict(mdic, orient="index")
    output_file = os.path.join(
        MAIN_PATH, "tables", f"{timetag}_{anomtag}_table_JASMIN_3hmeansVersion_rainMask.csv"
    )
    df.to_csv(output_file)
    print(f"Saved: {output_file}")


def main():
    os.makedirs(os.path.join(MAIN_PATH, "tables"), exist_ok=True)

    for timetag in TIMETAGS:
        for anomtag in ANOMTAGS:
            run_case(timetag, anomtag, processes=4)


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
