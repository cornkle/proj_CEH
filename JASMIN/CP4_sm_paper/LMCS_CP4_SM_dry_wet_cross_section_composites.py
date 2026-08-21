import glob
import multiprocessing as mp
import os

import numpy as np
import pandas as pd
import xarray as xr
from metpy import calc
from metpy.units import units


MAIN_PATH = "/gws/ssde/j25b/lmcs/cklein/CP_models/MCS_files/WAf/CP4_box_JASMIN/mean3h_v2"
OUTDIR = os.path.join(MAIN_PATH, "composites")
HOURS, DIRECTIONS = (17, 18), ("XDIR", "YDIR")
PROCESSES, CHUNK_SIZE, DX_M = 4, 250, 4400


def parse_key(path_or_name):
    """Return a key shared by 2-D and XDIR/YDIR filenames."""
    p = os.path.basename(path_or_name).split("_")
    if p[2] in DIRECTIONS:  # date_time_XDIR_ID_lonXlat_lon_lat.nc
        date, time, storm_id, lon, lat = p[0], p[1], p[3], p[5], p[6]
    else:                   # date_time_ID_lonXlat_lon_lat_3Hmeans.nc
        date, time, storm_id, lon, lat = p[0], p[1], p[2], p[4], p[5]
    return date, time, int(storm_id), float(lon.strip("[]")), float(lat.replace(".nc", "").strip("[]"))


def file_direction(path):
    return os.path.basename(path).split("_")[2]


def file_hour(path):
    return int(os.path.basename(path).split("_")[1].split(":")[0])


def read_groups(climate):
    """Read the exact DRY/WET membership produced by the 2-D script."""
    path = os.path.join(OUTDIR, f"{climate}_SM_composite_classification_17-18UTC.csv")
    table = pd.read_csv(path)
    table = table[table["classification"].isin(["dry", "wet"])].copy()
    groups = {parse_key(row.file): row.classification for row in table.itertuples()}
    print(f"{climate}: loaded {len(groups)} retained DRY/WET classifications from {path}")
    return groups


def add_derived_fields(ds, kind, direction):
    ds = ds.copy()
    hdim = "longitude" if direction == "XDIR" else "latitude"
    wind = "u_cross" if direction == "XDIR" else "v_cross"

    # Along-section horizontal divergence component only, not full divergence.
    if wind in ds:
        deriv = np.gradient(ds[wind].values, DX_M, axis=ds[wind].get_axis_num(hdim))
        ds["div_horizontal_cross"] = (ds[wind].dims, deriv)
        ds["div_horizontal_cross"].attrs.update(
            units="s-1", long_name=f"along-section horizontal divergence component ({'du/dx' if direction == 'XDIR' else 'dv/dy'})")

    # Nonlinear thermodynamics require absolute T and q, so use mean files only.
    if kind != "mean" or not {"t_cross", "q_cross", "pressure"}.issubset(ds.variables):
        return ds

    pressure = ds["pressure"].values[:, None] * units.hPa
    temperature = ds["t_cross"].values * units.K
    specific_humidity = ds["q_cross"].values * units("kg/kg")
    rh = calc.relative_humidity_from_specific_humidity(pressure, temperature, specific_humidity)
    dewpoint = calc.dewpoint_from_specific_humidity(pressure, temperature, specific_humidity)
    thetae = calc.equivalent_potential_temperature(pressure, temperature, dewpoint).to("K")
    thetaes = calc.saturation_equivalent_potential_temperature(pressure, temperature).to("K")

    ds["rh_cross"] = (ds["t_cross"].dims, rh.to("dimensionless").magnitude * 100)
    ds["thetae_cross"] = (ds["t_cross"].dims, thetae.magnitude)
    ds["thetaes_cross"] = (ds["t_cross"].dims, thetaes.magnitude)
    ds["rh_cross"].attrs.update(units="%", long_name="relative humidity")
    ds["thetae_cross"].attrs.update(units="K", long_name="equivalent potential temperature")
    ds["thetaes_cross"].attrs.update(units="K", long_name="saturation equivalent potential temperature")
    return ds


def select_vars(ds, direction):
    allowed = {"pressure", "longitude" if direction == "XDIR" else "latitude"}
    return [vn for vn in ds.data_vars if ds[vn].ndim > 0 and set(ds[vn].dims).issubset(allowed)]


def new_acc(templates, varnames):
    return {kind: {group: {"sum": {vn: np.zeros(templates[kind][vn].shape) for vn in varnames[kind]},
                           "count": {vn: np.zeros(templates[kind][vn].shape, dtype=np.int64) for vn in varnames[kind]}, "n": 0}
                   for group in ("dry", "wet")} for kind in ("anom", "mean")}


def add_case(acc, ds, kind, group, varnames):
    acc[kind][group]["n"] += 1
    for vn in varnames[kind]:
        values = np.asarray(ds[vn].values, dtype=float)
        valid = np.isfinite(values)
        acc[kind][group]["sum"][vn][valid] += values[valid]
        acc[kind][group]["count"][vn][valid] += 1


def process_chunk(args):
    pairs, key_groups, direction, templates, varnames = args
    acc, records = new_acc(templates, varnames), []
    for key, afile, mfile in pairs:
        group, name = key_groups[key], os.path.basename(afile)
        try:
            with xr.open_dataset(afile) as ar, xr.open_dataset(mfile) as mr:
                anom = add_derived_fields(ar, "anom", direction)
                mean = add_derived_fields(mr, "mean", direction)
                add_case(acc, anom, "anom", group, varnames)
                add_case(acc, mean, "mean", group, varnames)
                records.append((name, group, "kept"))
        except Exception as exc:
            records.append((name, group, f"error: {type(exc).__name__}: {exc}"))
    return acc, records


def matched_files(climate, direction, groups):
    adir, mdir = os.path.join(MAIN_PATH, f"pl_anom_{climate}"), os.path.join(MAIN_PATH, f"pl_mean_{climate}")
    afiles = [f for f in glob.glob(os.path.join(adir, "*.nc")) if file_direction(f) == direction and file_hour(f) in HOURS]
    mfiles = [f for f in glob.glob(os.path.join(mdir, "*.nc")) if file_direction(f) == direction and file_hour(f) in HOURS]
    amap, mmap = {parse_key(f): f for f in afiles}, {parse_key(f): f for f in mfiles}
    selected = sorted(set(amap) & set(mmap) & set(groups))
    pairs = [(key, amap[key], mmap[key]) for key in selected]
    stats = {"anom_total": len(amap), "mean_total": len(mmap), "matched_classified": len(pairs),
             "classified_without_anom_cross": len(set(groups) - set(amap)),
             "classified_without_mean_cross": len((set(groups) & set(amap)) - set(mmap))}
    return pairs, stats


def run_direction(climate, direction, groups):
    pairs, stats = matched_files(climate, direction, groups)
    print(f"{climate} {direction}: {stats}")
    if not pairs:
        print(f"No matched classified files for {climate} {direction}; skipping")
        return

    with xr.open_dataset(pairs[0][1]) as ar, xr.open_dataset(pairs[0][2]) as mr:
        templates = {"anom": add_derived_fields(ar, "anom", direction).load(),
                     "mean": add_derived_fields(mr, "mean", direction).load()}
    varnames = {kind: select_vars(templates[kind], direction) for kind in ("anom", "mean")}
    tasks = [(pairs[i:i + CHUNK_SIZE], groups, direction, templates, varnames) for i in range(0, len(pairs), CHUNK_SIZE)]
    with mp.Pool(PROCESSES) as pool:
        results = pool.map(process_chunk, tasks)

    total, records = new_acc(templates, varnames), []
    for acc, rec in results:
        records.extend(rec)
        for kind in ("anom", "mean"):
            for group in ("dry", "wet"):
                total[kind][group]["n"] += acc[kind][group]["n"]
                for vn in varnames[kind]:
                    total[kind][group]["sum"][vn] += acc[kind][group]["sum"][vn]
                    total[kind][group]["count"][vn] += acc[kind][group]["count"][vn]

    hdim = "longitude" if direction == "XDIR" else "latitude"
    coords = {"pressure": templates["anom"]["pressure"].values, hdim: templates["anom"][hdim].values}
    for kind in ("anom", "mean"):
        for group in ("dry", "wet"):
            data_vars = {}
            for vn in varnames[kind]:
                sums, counts = total[kind][group]["sum"][vn], total[kind][group]["count"][vn]
                avg = np.full(sums.shape, np.nan); np.divide(sums, counts, out=avg, where=counts > 0)
                dims = templates[kind][vn].dims
                data_vars[vn] = (dims, avg, dict(templates[kind][vn].attrs))
                data_vars[f"{vn}_n_valid"] = (dims, counts)
            out = xr.Dataset(data_vars, coords=coords)
            out.attrs.update(climate=climate, dataset_type=kind, direction=direction, soil_moisture_class=group.upper(),
                             hours="17,18", n_storms=total[kind][group]["n"],
                             classification_source=f"{climate}_SM_composite_classification_17-18UTC.csv")
            path = os.path.join(OUTDIR, f"{climate}_{kind}_{group}_{direction}_cross_section_17-18UTC.nc")
            out.to_netcdf(path, encoding={vn: {"zlib": True, "complevel": 4} for vn in out.data_vars})
            print(f"Saved {path} ({total[kind][group]['n']} storms)")

    pd.DataFrame(records, columns=["cross_section_file", "classification", "status"]).to_csv(
        os.path.join(OUTDIR, f"{climate}_{direction}_cross_section_audit_17-18UTC.csv"), index=False)


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    for climate in ("hist", "fut"):
        groups = read_groups(climate)
        for direction in DIRECTIONS:
            run_direction(climate, direction, groups)


if __name__ == "__main__":
    mp.freeze_support()
    main()
