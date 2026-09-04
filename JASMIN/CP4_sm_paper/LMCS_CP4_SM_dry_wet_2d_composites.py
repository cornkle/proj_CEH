"""DRY/WET composites. Keep this beside CP4_all_storm_2d_composites.py."""
import multiprocessing as mp
import os

import numpy as np
import pandas as pd
import xarray as xr

from CP4_all_storm_2d_composites import (MAIN_PATH, PROCESSES, CHUNK_SIZE, add_derived_fields,
                                         select_vars, matched_pairs)

THRESHOLDS = {"hist": {"dry": -3.0608628459778258, "wet": 1.2792528717828768},
              "fut": {"dry": -4.436951070853239, "wet": 0.08351219035011467}}

SELECTIONS = ("all", "area_lt10000", "vmid_2deg_negative", "vmid_full_negative",
              "vmid_full_negative_area_lt10000")
SELECTION_DESCRIPTIONS = {
    "all": "no additional wind or storm-area filter",
    "area_lt10000": "cold-cloud area < 10000 km2",
    "vmid_2deg_negative": "absolute mean-file v_mid_2deg < 0 m s-1",
    "vmid_full_negative": "absolute mean-file full-box mean v_mid < 0 m s-1",
    "vmid_full_negative_area_lt10000": "full-box mean v_mid < 0 m s-1 and cold-cloud area < 10000 km2",
}


def new_acc(shape, varnames):
    return {selection: {kind: {group: {"sum": {vn: np.zeros(shape) for vn in varnames[kind]},
                                      "count": {vn: np.zeros(shape, dtype=np.int64) for vn in varnames[kind]}, "n": 0}
                              for group in ("dry", "wet")} for kind in ("anom", "mean")}
            for selection in SELECTIONS}


def add_case(acc, ds, selection, kind, group, varnames):
    acc[selection][kind][group]["n"] += 1
    for vn in varnames[kind]:
        values = np.asarray(ds[vn].values, dtype=float); valid = np.isfinite(values)
        acc[selection][kind][group]["sum"][vn][valid] += values[valid]
        acc[selection][kind][group]["count"][vn][valid] += 1


def process_chunk(args):
    pairs, climate, shape, varnames = args
    acc, records = new_acc(shape, varnames), []
    for afile, mfile in pairs:
        name = os.path.basename(afile)
        try:
            with xr.open_dataset(afile) as ar, xr.open_dataset(mfile) as mr:
                noon, rain, olr = ar.lsRain_noon, ar.lsRain, ar.lw_out_PBLtop
                central = noon.isel(latitude=slice(54, 62), longitude=slice(54, 62)).values
                if np.any(np.isfinite(central)) and np.nanmax(central) > 0.25:
                    records.append((name, np.nan, np.nan, np.nan, np.nan, "rejected_noon_rain", False, False, False, False)); continue
                if int((olr <= -50).sum()) <= 3:
                    records.append((name, np.nan, np.nan, np.nan, np.nan, "rejected_small_storm", False, False, False, False)); continue
                if not np.any(np.isfinite(rain.values)):
                    records.append((name, np.nan, np.nan, np.nan, np.nan, "rejected_no_storm_rain", False, False, False, False)); continue

                varmask = np.isfinite(noon) & (noon <= 0)
                smbox = ar.SM.where(varmask).sel(latitude=slice(-22, 22), longitude=slice(-22, 22))
                if int(smbox.count()) <= 0.5 * (44 ** 2):
                    records.append((name, np.nan, np.nan, np.nan, np.nan, "rejected_insufficient_SM", False, False, False, False)); continue
                sm2 = float(smbox.mean(skipna=True))
                group = "dry" if sm2 <= THRESHOLDS[climate]["dry"] else "wet" if sm2 >= THRESHOLDS[climate]["wet"] else None
                if group is None:
                    records.append((name, sm2, np.nan, np.nan, np.nan, "middle", False, False, False, False)); continue

                # Reproduce v_mid_2deg from the absolute mean file. This uses
                # the mean file's rainfall-free mask, as in the mean table.
                mean_varmask = np.isfinite(mr.lsRain_noon) & (mr.lsRain_noon <= 0)
                vbox = mr.v_mid.where(mean_varmask).sel(latitude=slice(-22, 22), longitude=slice(-22, 22))
                vmid2 = float(vbox.mean(skipna=True)) if int(vbox.count()) > 0.5 * (44 ** 2) else np.nan
                vfull = mr.v_mid.where(mean_varmask)
                vmid_full = float(vfull.mean(skipna=True)) if int(vfull.count()) > 0.5 * vfull.size else np.nan
                storm_area_km2 = float((olr <= -50).sum()) * (4400 ** 2) / 1e6
                flags = {"all": True, "area_lt10000": storm_area_km2 < 10000,
                         "vmid_2deg_negative": np.isfinite(vmid2) and vmid2 < 0,
                         "vmid_full_negative": np.isfinite(vmid_full) and vmid_full < 0,
                         "vmid_full_negative_area_lt10000": np.isfinite(vmid_full) and vmid_full < 0 and storm_area_km2 < 10000}

                # Imported function uses K and kg kg-1, and calculates nonlinear
                # thermodynamics only from the absolute mean file.
                anom, mean = add_derived_fields(ar, "anom"), add_derived_fields(mr, "mean")
                for selection, keep in flags.items():
                    if keep:
                        add_case(acc, anom, selection, "anom", group, varnames)
                        add_case(acc, mean, selection, "mean", group, varnames)
                records.append((name, sm2, vmid2, vmid_full, storm_area_km2, group, flags["area_lt10000"],
                                flags["vmid_2deg_negative"], flags["vmid_full_negative"], flags["vmid_full_negative_area_lt10000"]))
        except Exception as exc:
            records.append((name, np.nan, np.nan, np.nan, np.nan, f"error: {type(exc).__name__}: {exc}", False, False, False, False))
    return acc, records


def run_climate(climate):
    pairs, anom_only, mean_only = matched_pairs(climate)  # already restricted to matched 17/18 UTC files
    if not pairs:
        raise FileNotFoundError(f"No matched 17/18 UTC files for {climate}")
    print(f"{climate}: matched={len(pairs)}, anomaly-only={len(anom_only)}, mean-only={len(mean_only)}")
    with xr.open_dataset(pairs[0][0]) as ar, xr.open_dataset(pairs[0][1]) as mr:
        templates = {"anom": add_derived_fields(ar, "anom").load(), "mean": add_derived_fields(mr, "mean").load()}
    varnames = {kind: select_vars(templates[kind]) for kind in templates}
    shape = (templates["anom"].sizes["latitude"], templates["anom"].sizes["longitude"])
    tasks = [(pairs[i:i + CHUNK_SIZE], climate, shape, varnames) for i in range(0, len(pairs), CHUNK_SIZE)]
    with mp.Pool(PROCESSES) as pool:
        results = pool.map(process_chunk, tasks)

    total, records = new_acc(shape, varnames), []
    for acc, rec in results:
        records.extend(rec)
        for selection in SELECTIONS:
            for kind in ("anom", "mean"):
                for group in ("dry", "wet"):
                    total[selection][kind][group]["n"] += acc[selection][kind][group]["n"]
                    for vn in varnames[kind]:
                        total[selection][kind][group]["sum"][vn] += acc[selection][kind][group]["sum"][vn]
                        total[selection][kind][group]["count"][vn] += acc[selection][kind][group]["count"][vn]

    outdir = os.path.join(MAIN_PATH, "composites"); os.makedirs(outdir, exist_ok=True)
    for selection in SELECTIONS:
        for kind in ("anom", "mean"):
            for group in ("dry", "wet"):
                data_vars = {}
                for vn in varnames[kind]:
                    sums = total[selection][kind][group]["sum"][vn]
                    counts = total[selection][kind][group]["count"][vn]
                    avg = np.full(shape, np.nan); np.divide(sums, counts, out=avg, where=counts > 0)
                    data_vars[vn] = (("latitude", "longitude"), avg, dict(templates[kind][vn].attrs))
                    data_vars[f"{vn}_n_valid"] = (("latitude", "longitude"), counts)
                out = xr.Dataset(data_vars, coords={"latitude": templates[kind].latitude, "longitude": templates[kind].longitude})
                out.attrs.update(climate=climate, dataset_type=kind, soil_moisture_class=group.upper(), hours="17,18",
                                 n_storms=total[selection][kind][group]["n"], sm_threshold=THRESHOLDS[climate][group],
                                 sm_sampling="Anomaly-file SM_2deg over -22:22 after lsRain_noon <= 0 mask",
                                 wind_selection=selection,
                                 additional_selection=SELECTION_DESCRIPTIONS[selection])
                suffix = "" if selection == "all" else f"_{selection}"
                path = os.path.join(outdir, f"{climate}_{kind}_{group}_2d_composite_17-18UTC{suffix}.nc")
                out.to_netcdf(path, encoding={vn: {"zlib": True, "complevel": 4} for vn in out.data_vars})
                print(f"Saved {path} ({total[selection][kind][group]['n']} storms)")

    audit_columns = ["file", "SM_2deg", "v_mid_2deg", "v_mid_full", "storm_area_km2", "classification",
                     "area_lt10000", "vmid_2deg_negative", "vmid_full_negative", "vmid_full_negative_area_lt10000"]
    pd.DataFrame(records, columns=audit_columns).to_csv(
        os.path.join(outdir, f"{climate}_SM_composite_classification_17-18UTC.csv"), index=False)
    mismatch = [(n, "anomaly_without_mean") for n in anom_only] + [(n, "mean_without_anomaly") for n in mean_only]
    pd.DataFrame(mismatch, columns=["file", "status"]).to_csv(
        os.path.join(outdir, f"{climate}_mean_anomaly_mismatches_17-18UTC.csv"), index=False)


def main():
    for climate in ("hist", "fut"):
        run_climate(climate)


if __name__ == "__main__":
    mp.freeze_support(); main()