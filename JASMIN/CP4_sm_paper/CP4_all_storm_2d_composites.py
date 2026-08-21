import glob
import multiprocessing as mp
import os

import numpy as np
import pandas as pd
import xarray as xr
from metpy import calc
from metpy.units import units


MAIN_PATH = "/gws/ssde/j25b/lmcs/cklein/CP_models/MCS_files/WAf/CP4_box_JASMIN/mean3h_v2"
HOURS = (17, 18)
PROCESSES = 4
CHUNK_SIZE = 250
DX_M = 4400
COMPOSITE_VARS = None  # None = every 2-D variable


def filename_hour(path):
    return int(os.path.basename(path).split("_")[1].split(":")[0])


def add_derived_fields(ds, dataset_type):
    """Divergence for both datasets; nonlinear thermodynamics for absolute means only."""
    ds = ds.copy()
    if {"u_srfc", "v_srfc"}.issubset(ds.data_vars):
        u, v = ds["u_srfc"].values * units("m/s"), ds["v_srfc"].values * units("m/s")
        div = calc.divergence(u, v, dx=DX_M * units.m, dy=DX_M * units.m).to("1/s")
        ds["div"] = (ds["u_srfc"].dims, div.magnitude)
        ds["div"].attrs.update(units="s-1", long_name="925-hPa horizontal wind divergence")

    if dataset_type != "mean" or not {"t_srfc", "q_srfc", "t_mid", "q_mid"}.issubset(ds.data_vars):
        return ds

    # File units: temperature K and specific humidity kg kg-1.
    tl, tm = ds["t_srfc"].values * units.K, ds["t_mid"].values * units.K
    ql, qm = ds["q_srfc"].values * units("kg/kg"), ds["q_mid"].values * units("kg/kg")
    pl, pm = 925 * units.hPa, 650 * units.hPa
    rhl = calc.relative_humidity_from_specific_humidity(pl, tl, ql)
    rhm = calc.relative_humidity_from_specific_humidity(pm, tm, qm)
    tdl = calc.dewpoint_from_specific_humidity(pl, tl, ql)
    tdm = calc.dewpoint_from_specific_humidity(pm, tm, qm)
    thel = calc.equivalent_potential_temperature(pl, tl, tdl).to("K")
    them = calc.equivalent_potential_temperature(pm, tm, tdm).to("K")
    thes_m = calc.saturation_equivalent_potential_temperature(pm, tm).to("K")

    derived = {
        "rh_srfc": (ds["q_srfc"].dims, rhl.to("dimensionless").magnitude * 100, "%"),
        "rh_mid": (ds["q_mid"].dims, rhm.to("dimensionless").magnitude * 100, "%"),
        "thetae_srfc": (ds["t_srfc"].dims, thel.magnitude, "K"),
        "thetae_mid": (ds["t_mid"].dims, them.magnitude, "K"),
        "thetaes_mid": (ds["t_mid"].dims, thes_m.magnitude, "K"),
        "cape_proxy": (ds["t_srfc"].dims, (thel - thes_m).magnitude, "K"),
        "cin_proxy": (ds["t_mid"].dims, (thes_m - them).magnitude, "K"),
    }
    for vn, (dims, values, unit) in derived.items():
        ds[vn] = (dims, values)
        ds[vn].attrs["units"] = unit
    return ds


def select_vars(ds):
    names = [vn for vn in ds.data_vars if ds[vn].dims == ("latitude", "longitude")]
    return names if COMPOSITE_VARS is None else [vn for vn in COMPOSITE_VARS if vn in names]


def new_acc(shape, varnames):
    return {kind: {"sum": {vn: np.zeros(shape) for vn in varnames[kind]},
                   "count": {vn: np.zeros(shape, dtype=np.int64) for vn in varnames[kind]}, "n": 0}
            for kind in ("anom", "mean")}


def add_case(acc, ds, kind, varnames):
    acc[kind]["n"] += 1
    for vn in varnames[kind]:
        values = np.asarray(ds[vn].values, dtype=np.float64)
        valid = np.isfinite(values)
        acc[kind]["sum"][vn][valid] += values[valid]
        acc[kind]["count"][vn][valid] += 1


def process_chunk(args):
    pairs, shape, varnames = args
    acc, records = new_acc(shape, varnames), []
    for afile, mfile in pairs:
        name = os.path.basename(afile)
        try:
            with xr.open_dataset(afile) as ar, xr.open_dataset(mfile) as mr:
                # Use the anomaly file for the same case-selection filters as the DRY/WET script.
                noon, rain, olr = ar["lsRain_noon"], ar["lsRain"], ar["lw_out_PBLtop"]
                central = noon.isel(latitude=slice(54, 62), longitude=slice(54, 62)).values
                if np.any(np.isfinite(central)) and np.nanmax(central) > 0.25:
                    records.append((name, "rejected_noon_rain")); continue
                if int((olr <= -50).sum()) <= 3:
                    records.append((name, "rejected_small_storm")); continue
                if not np.any(np.isfinite(rain.values)):
                    records.append((name, "rejected_no_storm_rain")); continue

                anom, mean = add_derived_fields(ar, "anom"), add_derived_fields(mr, "mean")
                add_case(acc, anom, "anom", varnames)
                add_case(acc, mean, "mean", varnames)
                records.append((name, "kept"))
        except Exception as exc:
            records.append((name, f"error: {type(exc).__name__}: {exc}"))
    return acc, records


def matched_pairs(climate):
    adir, mdir = os.path.join(MAIN_PATH, f"anom_{climate}"), os.path.join(MAIN_PATH, f"mean_{climate}")
    amap = {os.path.basename(f): f for f in glob.glob(os.path.join(adir, "*.nc")) if filename_hour(f) in HOURS}
    mmap = {os.path.basename(f): f for f in glob.glob(os.path.join(mdir, "*.nc")) if filename_hour(f) in HOURS}
    common = sorted(set(amap) & set(mmap))
    return [(amap[n], mmap[n]) for n in common], sorted(set(amap) - set(mmap)), sorted(set(mmap) - set(amap))


def run_climate(climate):
    pairs, anom_only, mean_only = matched_pairs(climate)
    if not pairs:
        raise FileNotFoundError(f"No matched 17/18 UTC files for {climate}")
    print(f"{climate}: matched={len(pairs)}, anomaly-only={len(anom_only)}, mean-only={len(mean_only)}")

    with xr.open_dataset(pairs[0][0]) as ar, xr.open_dataset(pairs[0][1]) as mr:
        templates = {"anom": add_derived_fields(ar, "anom").load(), "mean": add_derived_fields(mr, "mean").load()}
    varnames = {kind: select_vars(templates[kind]) for kind in templates}
    shape = (templates["anom"].sizes["latitude"], templates["anom"].sizes["longitude"])
    tasks = [(pairs[i:i + CHUNK_SIZE], shape, varnames) for i in range(0, len(pairs), CHUNK_SIZE)]
    with mp.Pool(PROCESSES) as pool:
        results = pool.map(process_chunk, tasks)

    total, records = new_acc(shape, varnames), []
    for acc, rec in results:
        records.extend(rec)
        for kind in ("anom", "mean"):
            total[kind]["n"] += acc[kind]["n"]
            for vn in varnames[kind]:
                total[kind]["sum"][vn] += acc[kind]["sum"][vn]
                total[kind]["count"][vn] += acc[kind]["count"][vn]

    outdir = os.path.join(MAIN_PATH, "composites")
    os.makedirs(outdir, exist_ok=True)
    for kind in ("anom", "mean"):
        data_vars = {}
        for vn in varnames[kind]:
            sums, counts = total[kind]["sum"][vn], total[kind]["count"][vn]
            avg = np.full(shape, np.nan); np.divide(sums, counts, out=avg, where=counts > 0)
            data_vars[vn] = (("latitude", "longitude"), avg, dict(templates[kind][vn].attrs))
            data_vars[f"{vn}_n_valid"] = (("latitude", "longitude"), counts)
        out = xr.Dataset(data_vars, coords={"latitude": templates[kind].latitude, "longitude": templates[kind].longitude})
        out.attrs.update(climate=climate, dataset_type=kind, hours="17,18", n_storms=total[kind]["n"],
                         selection="Complete filename present in mean and anomaly directories; validity filters from anomaly file")
        path = os.path.join(outdir, f"{climate}_{kind}_all_2d_composite_17-18UTC.nc")
        out.to_netcdf(path, encoding={vn: {"zlib": True, "complevel": 4} for vn in out.data_vars})
        print(f"Saved {path} ({total[kind]['n']} storms)")

    pd.DataFrame(records, columns=["file", "classification"]).to_csv(
        os.path.join(outdir, f"{climate}_all_composite_audit_17-18UTC.csv"), index=False)
    mismatch = [(n, "anomaly_without_mean") for n in anom_only] + [(n, "mean_without_anomaly") for n in mean_only]
    pd.DataFrame(mismatch, columns=["file", "status"]).to_csv(
        os.path.join(outdir, f"{climate}_mean_anomaly_mismatches_17-18UTC.csv"), index=False)


def main():
    for climate in ("hist", "fut"):
        run_climate(climate)


if __name__ == "__main__":
    mp.freeze_support()
    main()
