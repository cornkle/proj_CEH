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


def new_acc(shape, varnames):
    return {kind: {group: {"sum": {vn: np.zeros(shape) for vn in varnames[kind]},
                           "count": {vn: np.zeros(shape, dtype=np.int64) for vn in varnames[kind]}, "n": 0}
                   for group in ("dry", "wet")} for kind in ("anom", "mean")}


def add_case(acc, ds, kind, group, varnames):
    acc[kind][group]["n"] += 1
    for vn in varnames[kind]:
        values = np.asarray(ds[vn].values, dtype=float); valid = np.isfinite(values)
        acc[kind][group]["sum"][vn][valid] += values[valid]
        acc[kind][group]["count"][vn][valid] += 1


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
                    records.append((name, np.nan, "rejected_noon_rain")); continue
                if int((olr <= -50).sum()) <= 3:
                    records.append((name, np.nan, "rejected_small_storm")); continue
                if not np.any(np.isfinite(rain.values)):
                    records.append((name, np.nan, "rejected_no_storm_rain")); continue

                varmask = np.isfinite(noon) & (noon <= 0)
                smbox = ar.SM.where(varmask).sel(latitude=slice(-22, 22), longitude=slice(-22, 22))
                if int(smbox.count()) <= 0.5 * (44 ** 2):
                    records.append((name, np.nan, "rejected_insufficient_SM")); continue
                sm2 = float(smbox.mean(skipna=True))
                group = "dry" if sm2 <= THRESHOLDS[climate]["dry"] else "wet" if sm2 >= THRESHOLDS[climate]["wet"] else None
                if group is None:
                    records.append((name, sm2, "middle")); continue

                # Imported function uses K and kg kg-1, and calculates nonlinear
                # thermodynamics only from the absolute mean file.
                anom, mean = add_derived_fields(ar, "anom"), add_derived_fields(mr, "mean")
                add_case(acc, anom, "anom", group, varnames)
                add_case(acc, mean, "mean", group, varnames)
                records.append((name, sm2, group))
        except Exception as exc:
            records.append((name, np.nan, f"error: {type(exc).__name__}: {exc}"))
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
        for kind in ("anom", "mean"):
            for group in ("dry", "wet"):
                total[kind][group]["n"] += acc[kind][group]["n"]
                for vn in varnames[kind]:
                    total[kind][group]["sum"][vn] += acc[kind][group]["sum"][vn]
                    total[kind][group]["count"][vn] += acc[kind][group]["count"][vn]

    outdir = os.path.join(MAIN_PATH, "composites"); os.makedirs(outdir, exist_ok=True)
    for kind in ("anom", "mean"):
        for group in ("dry", "wet"):
            data_vars = {}
            for vn in varnames[kind]:
                sums, counts = total[kind][group]["sum"][vn], total[kind][group]["count"][vn]
                avg = np.full(shape, np.nan); np.divide(sums, counts, out=avg, where=counts > 0)
                data_vars[vn] = (("latitude", "longitude"), avg, dict(templates[kind][vn].attrs))
                data_vars[f"{vn}_n_valid"] = (("latitude", "longitude"), counts)
            out = xr.Dataset(data_vars, coords={"latitude": templates[kind].latitude, "longitude": templates[kind].longitude})
            out.attrs.update(climate=climate, dataset_type=kind, soil_moisture_class=group.upper(), hours="17,18",
                             n_storms=total[kind][group]["n"], sm_threshold=THRESHOLDS[climate][group],
                             sm_sampling="Anomaly-file SM_2deg over -22:22 after lsRain_noon <= 0 mask")
            path = os.path.join(outdir, f"{climate}_{kind}_{group}_2d_composite_17-18UTC.nc")
            out.to_netcdf(path, encoding={vn: {"zlib": True, "complevel": 4} for vn in out.data_vars})
            print(f"Saved {path} ({total[kind][group]['n']} storms)")

    pd.DataFrame(records, columns=["file", "SM_2deg", "classification"]).to_csv(
        os.path.join(outdir, f"{climate}_SM_composite_classification_17-18UTC.csv"), index=False)
    mismatch = [(n, "anomaly_without_mean") for n in anom_only] + [(n, "mean_without_anomaly") for n in mean_only]
    pd.DataFrame(mismatch, columns=["file", "status"]).to_csv(
        os.path.join(outdir, f"{climate}_mean_anomaly_mismatches_17-18UTC.csv"), index=False)


def main():
    for climate in ("hist", "fut"):
        run_climate(climate)


if __name__ == "__main__":
    mp.freeze_support(); main()