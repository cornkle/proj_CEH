#!/usr/bin/env python3
"""Add lower-tropospheric and full ICAPE to the final CP4 MCS tables.

The existing hist/fut mean tables define the cases and storm centres. Raw CP4
T/q profiles plus 2-m T/q and surface pressure are sampled at 12 UTC. The
matching 2-D storm file supplies the same no-noon-rain mask used by the table.

For every fixed spatial scale (0.25deg, 1deg, 2deg and fullBox), the script:
  1. forms an area-mean environmental sounding;
  2. calculates CAPE for parcels originating from the surface through 800 hPa;
  3. integrates parcel CAPE over source pressure to give J m-2.

ICAPE500 limits each parcel's buoyancy integral at 500 hPa. ICAPEfull uses all
available pressure levels. Work is checkpointed by date and safely resumable.
Original tables are never overwritten unless an explicit output path is used.
"""

from __future__ import annotations

import argparse
import glob
import multiprocessing
import os
import re
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from metpy import calc
from metpy.units import units

from JASMIN import MetUM_variables as mu
from shared.utils import u_interpolate as u_int


RAW_ROOT = Path("/home/users/cornkle/linked_CP4")
BOX_ROOT = Path("/gws/ssde/j25b/lmcs/cklein/CP_models/MCS_files/WAf/CP4_box_JASMIN/mean3h_v2")
TABLE_DIR = BOX_ROOT / "tables"
DOMAIN = (-18.0, 25.0, 5.0, 25.0)
ENV_HOUR = 12
SOURCE_TOP_HPA = 800.0
PARTIAL_TOP_HPA = 500.0
HALF_FULL_BOX = 57
MIN_VALID_FRACTION = 0.5
BOXES = {"0.25deg": 3, "1deg": 11, "2deg": 22, "fullBox": None}
OUTPUT_COLUMNS = [f"{kind}_{scale}" for kind in ("ICAPE500", "ICAPEfull") for scale in BOXES]
NUMBER_RE = re.compile(r"[-+]?\d+(?:\.\d+)?")
LONLAT_RE = re.compile(r"_lonXlat_\[?([-+]?\d+(?:\.\d+)?)\]?_\[?([-+]?\d+(?:\.\d+)?)\]?")
G = 9.80665


@dataclass
class DayTask:
    tag: str
    date: str
    rows: list[dict]
    raw_root: str
    box_root: str
    daily_file: str
    surface_time_mode: str
    overwrite: bool


def number(value, name: str) -> float:
    """Extract one numeric scalar from normal or stringified table values."""
    if pd.isna(value):
        raise ValueError(f"Missing {name}")
    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value)
    match = NUMBER_RE.search(str(value))
    if not match:
        raise ValueError(f"Cannot parse {name}={value!r}")
    return float(match.group())


def read_table(path: Path) -> tuple[pd.DataFrame, bool]:
    """Read either variables-as-rows or cases-as-rows table orientation."""
    raw = pd.read_csv(path, index_col=0, low_memory=False)
    required = {"year", "month", "day", "hour", "lon", "lat"}
    if required.issubset(raw.index.astype(str)):
        return raw.T.copy(), True
    if required.issubset(raw.columns.astype(str)):
        return raw.copy(), False
    raise ValueError(f"Cannot identify table orientation or required fields in {path}")


def resolve_table(table_dir: Path, tag: str, explicit: Path | None) -> Path:
    if explicit is not None:
        if not explicit.exists():
            raise FileNotFoundError(explicit)
        return explicit
    exact = table_dir / f"{tag}_mean_table_JASMIN_3hmeansVersion_rainMask_fullBox.csv"
    if exact.exists():
        return exact
    candidates = sorted(p for p in table_dir.glob(
        f"{tag}_mean_table_JASMIN_3hmeansVersion_rainMask_fullBox*.csv"
    ) if "ICAPE" not in p.stem and "failure" not in p.stem)
    if len(candidates) != 1:
        raise FileNotFoundError(f"Expected one {tag} mean fullBox table, found {len(candidates)}: {candidates}")
    return candidates[0]


def table_records(frame: pd.DataFrame) -> list[dict]:
    records = []
    for case_id, (_, row) in enumerate(frame.iterrows()):
        rec = row.to_dict()
        rec["case_id"] = case_id
        for field in ("year", "month", "day", "hour", "lon", "lat"):
            rec[field] = number(rec[field], field)
        for field in ("index", "indx", "indy"):
            if field in rec and not pd.isna(rec[field]):
                rec[field] = number(rec[field], field)
        rec["date"] = f"{int(rec['year']):04d}{int(rec['month']):02d}{int(rec['day']):02d}"
        records.append(rec)
    return records


def raw_file(raw_dir: Path, var: str, date: str) -> Path:
    matches = sorted(Path(p) for p in glob.glob(str(raw_dir / var / f"*_{date}*.nc")))
    if not matches:
        raise FileNotFoundError(f"No {var} file for {date} in {raw_dir / var}")
    if len(matches) > 1:
        raise RuntimeError(f"Multiple {var} files for {date}: {matches}")
    return matches[0]


def select_time(da: xr.DataArray, hour: int, mode: str) -> xr.DataArray:
    if "time" not in da.dims:
        return da.squeeze()
    hours = da["time"].dt.hour
    if mode == "3hmean":
        chosen = da.where((hours >= hour - 2) & (hours <= hour), drop=True)
        if chosen.sizes.get("time", 0) == 0:
            raise ValueError(f"No {hour-2:02d}-{hour:02d} UTC data for {da.name}")
        return chosen.mean("time").squeeze()
    chosen = da.where(hours == hour, drop=True)
    if chosen.sizes.get("time", 0) == 0:
        raise ValueError(f"No {hour:02d} UTC data for {da.name}")
    if chosen.sizes.get("time", 1) != 1:
        raise ValueError(f"Expected one {hour:02d} UTC field for {da.name}, found {chosen.sizes['time']}")
    return chosen.squeeze("time", drop=True)


def load_raw(raw_dir: Path, var: str, date: str, mode: str) -> xr.DataArray:
    path = raw_file(raw_dir, var, date)
    with xr.open_dataset(path) as ds:
        name = mu.create_CP4_filename(var)
        if name not in ds:
            raise KeyError(f"{name!r} not found in {path}; variables={list(ds.data_vars)}")
        da = ds[name]
        if "longitude" in da.coords and float(da.longitude.max()) > 180:
            da = da.assign_coords(longitude=da.longitude - 360)
        da = da.sel(longitude=slice(DOMAIN[0], DOMAIN[1]), latitude=slice(DOMAIN[2], DOMAIN[3]))
        da = select_time(da, ENV_HOUR, mode).load()
    return da


def normalise_temperature(values: np.ndarray, label: str) -> np.ndarray:
    out = np.asarray(values, dtype=float)
    median = float(np.nanmedian(out))
    if median > 1000:
        out = out / 100.0
        median /= 100.0
    if not 150 < median < 350:
        raise ValueError(f"Implausible {label} median temperature: {median:g} K")
    return out


def normalise_humidity(values: np.ndarray, label: str) -> np.ndarray:
    out = np.asarray(values, dtype=float)
    positive = out[np.isfinite(out) & (out > 0)]
    median = float(np.nanmedian(positive)) if positive.size else np.nan
    if np.isfinite(median) and median > 0.2:
        out = out / 1000.0
        median /= 1000.0
    if not np.isfinite(median) or not 0 < median < 0.1:
        raise ValueError(f"Implausible {label} median specific humidity: {median:g} kg kg-1")
    return out


def normalise_surface_pressure(values: np.ndarray) -> np.ndarray:
    out = np.asarray(values, dtype=float)
    median = float(np.nanmedian(out))
    if median > 2000:
        out = out / 100.0
        median /= 100.0
    if not 700 < median < 1100:
        raise ValueError(f"Implausible surface-pressure median: {median:g} hPa")
    return out


def regrid_surface(da: xr.DataArray, target: xr.DataArray) -> np.ndarray:
    if (np.array_equal(da.longitude.values, target.longitude.values)
            and np.array_equal(da.latitude.values, target.latitude.values)):
        return np.asarray(da.values, dtype=float)
    inds, weights, shape = u_int.interpolation_weights(
        da.longitude, da.latitude, target.longitude, target.latitude
    )
    return np.asarray(u_int.interpolate_data(da.values, inds, weights, shape), dtype=float)


def load_day_fields(raw_root: Path, tag: str, date: str, surface_time_mode: str) -> dict:
    raw_dir = raw_root / ("hist" if tag == "hist" else "future")
    t_da = load_raw(raw_dir, "t_pl", date, "instant").sortby("pressure", ascending=False)
    q_da = load_raw(raw_dir, "q_pl", date, "instant").sortby("pressure", ascending=False)
    t_da = t_da.transpose("pressure", "latitude", "longitude")
    q_da = q_da.transpose("pressure", "latitude", "longitude")
    if not (np.array_equal(t_da.pressure.values, q_da.pressure.values)
            and np.array_equal(t_da.latitude.values, q_da.latitude.values)
            and np.array_equal(t_da.longitude.values, q_da.longitude.values)):
        raise ValueError(f"T/q pressure grids differ on {date}")

    t2_da = load_raw(raw_dir, "t2", date, surface_time_mode)
    q2_da = load_raw(raw_dir, "q2", date, surface_time_mode)
    ps_da = load_raw(raw_dir, "p_srfc", date, surface_time_mode)
    pressure = np.asarray(t_da.pressure.values, dtype=float)
    if float(np.nanmedian(pressure)) > 2000:
        pressure = pressure / 100.0

    t = normalise_temperature(t_da.values, "pressure-level")
    q = normalise_humidity(q_da.values, "pressure-level")
    t[t == 0], q[q == 0] = np.nan, np.nan
    return {"pressure": pressure, "temperature": t, "humidity": q,
            "t2": normalise_temperature(regrid_surface(t2_da, t_da), "2-m"),
            "q2": normalise_humidity(regrid_surface(q2_da, t_da), "2-m"),
            "ps": normalise_surface_pressure(regrid_surface(ps_da, t_da)),
            "longitude": np.asarray(t_da.longitude.values), "latitude": np.asarray(t_da.latitude.values)}


def candidate_location(path: Path) -> tuple[float, float] | None:
    match = LONLAT_RE.search(path.name)
    return (float(match.group(1)), float(match.group(2))) if match else None


def find_2d_file(row: dict, box_root: Path, tag: str) -> Path:
    date = f"{int(row['year']):04d}-{int(row['month']):02d}-{int(row['day']):02d}"
    stamp = f"{date}_{int(row['hour']):02d}:00:00"
    folder = box_root / f"mean_{tag}"
    if "index" in row and np.isfinite(row["index"]):
        candidates = sorted(folder.glob(f"{stamp}_{int(row['index'])}_lonXlat_*.nc"))
    else:
        candidates = []
    if not candidates:
        candidates = sorted(folder.glob(f"{stamp}_*_lonXlat_*.nc"))
    if not candidates:
        raise FileNotFoundError(f"No matching 2-D file for case {row['case_id']} ({stamp})")

    target_lon = row.get("indx", np.nan)
    target_lat = row.get("indy", np.nan)
    target_lon = target_lon if np.isfinite(target_lon) else round(row["lon"], 1)
    target_lat = target_lat if np.isfinite(target_lat) else round(row["lat"], 1)
    scored = []
    for path in candidates:
        location = candidate_location(path)
        if location is None:
            continue
        scored.append(((location[0] - target_lon) ** 2 + (location[1] - target_lat) ** 2, path))
    if not scored:
        if len(candidates) == 1:
            return candidates[0]
        raise RuntimeError(f"Cannot distinguish 2-D candidates for case {row['case_id']}: {candidates}")
    scored.sort(key=lambda item: item[0])
    if len(scored) > 1 and np.isclose(scored[0][0], scored[1][0]):
        raise RuntimeError(f"Ambiguous 2-D candidates for case {row['case_id']}: {scored[:2]}")
    if scored[0][0] > 0.05 ** 2 * 2:
        raise RuntimeError(f"2-D location mismatch for case {row['case_id']}: {scored[0]}")
    return scored[0][1]


def noon_valid_mask(path: Path) -> np.ndarray:
    with xr.open_dataset(path) as ds:
        if "lsRain_noon" not in ds:
            raise KeyError(f"lsRain_noon missing from {path}")
        rain = np.asarray(ds["lsRain_noon"].squeeze().values, dtype=float)
    if rain.ndim != 2:
        raise ValueError(f"Expected 2-D lsRain_noon in {path}, got {rain.shape}")
    return np.isfinite(rain) & (rain <= 0)


def centred_cut(fields: dict, lon: float, lat: float) -> dict:
    xpos = int(np.nanargmin(np.abs(fields["longitude"] - lon)))
    ypos = int(np.nanargmin(np.abs(fields["latitude"] - lat)))
    ys = slice(ypos - HALF_FULL_BOX, ypos + HALF_FULL_BOX + 1)
    xs = slice(xpos - HALF_FULL_BOX, xpos + HALF_FULL_BOX + 1)
    if ypos - HALF_FULL_BOX < 0 or xpos - HALF_FULL_BOX < 0:
        raise IndexError(f"Storm centre ({lon}, {lat}) too close to lower domain boundary")
    out = {"pressure": fields["pressure"]}
    out["temperature"] = fields["temperature"][:, ys, xs]
    out["humidity"] = fields["humidity"][:, ys, xs]
    for name in ("t2", "q2", "ps"):
        out[name] = fields[name][ys, xs]
    expected = HALF_FULL_BOX * 2 + 1
    if out["temperature"].shape[1:] != (expected, expected):
        raise IndexError(f"Incomplete raw fullBox around ({lon}, {lat}): {out['temperature'].shape}")
    return out


def area_mean_profile(cut: dict, valid: np.ndarray, distance: int | None) -> tuple:
    size = valid.shape[0]
    if valid.shape != (size, size) or size != HALF_FULL_BOX * 2 + 1:
        raise ValueError(f"2-D mask/raw cut shape mismatch: {valid.shape}")
    centre = size // 2
    sl = slice(None) if distance is None else slice(centre - distance, centre + distance + 1)
    mask = valid[sl, sl]
    if mask.sum() < MIN_VALID_FRACTION * mask.size:
        print(f"SKIP box={distance}: valid pixels={mask.sum()}/{mask.size} ({mask.mean():.1%})", flush=True)
        return (None,) * 5
      
    t = cut["temperature"][:, sl, sl]
    q = cut["humidity"][:, sl, sl]
  
    with np.errstate(invalid="ignore"):
        t_prof = np.nanmean(np.where(mask[None, :, :], t, np.nan), axis=(1, 2))
        q_prof = np.nanmean(np.where(mask[None, :, :], q, np.nan), axis=(1, 2))
        t2 = float(np.nanmean(np.where(mask, cut["t2"][sl, sl], np.nan)))
        q2 = float(np.nanmean(np.where(mask, cut["q2"][sl, sl], np.nan)))
        ps = float(np.nanmean(np.where(mask, cut["ps"][sl, sl], np.nan)))
    return cut["pressure"], t_prof, q_prof, ps, (t2, q2)


def cape_pair_for_sources(pressure, temperature, humidity, ps, surface_state) -> tuple[float, float]:
    if pressure is None or not np.all(np.isfinite([ps, *surface_state])):
        return np.nan, np.nan
    p, t, q = map(lambda x: np.asarray(x, dtype=float), (pressure, temperature, humidity))
    good = np.isfinite(p) & np.isfinite(t) & np.isfinite(q) & (p < ps - 0.05)
    p, t, q = p[good], t[good], q[good]
    if p.size < 5:
        return np.nan, np.nan
    order = np.argsort(p)[::-1]
    p, t, q = p[order], t[order], q[order]
    p = np.concatenate(([ps], p))
    t = np.concatenate(([surface_state[0]], t))
    q = np.concatenate(([surface_state[1]], q))
    if p.min() > PARTIAL_TOP_HPA or not np.any(np.isclose(p, SOURCE_TOP_HPA, atol=0.6)):
        return np.nan, np.nan

    p_units, t_units, q_units = p * units.hPa, t * units.K, q * units("kg/kg")
    td = calc.dewpoint_from_specific_humidity(p_units, t_units, q_units).to("K")
    td = np.minimum(td.magnitude, t) * units.K
    source_indices = np.where(p >= SOURCE_TOP_HPA - 0.01)[0]
    if source_indices.size < 2:
        return np.nan, np.nan

    source_p, cape500, capefull = [], [], []
    for index in source_indices:
        p_up, t_up, td_up = p_units[index:], t_units[index:], td[index:]
        if p_up.size < 3:
            return np.nan, np.nan
        parcel = calc.parcel_profile(p_up, t_units[index], td[index]).to("K")
        full, _ = calc.cape_cin(p_up, t_up, td_up, parcel, which_lfc="bottom", which_el="top")
        partial = p_up >= PARTIAL_TOP_HPA * units.hPa
        if int(partial.sum()) < 3:
            return np.nan, np.nan
        low, _ = calc.cape_cin(p_up[partial], t_up[partial], td_up[partial], parcel[partial],
                               which_lfc="bottom", which_el="top")
        source_p.append(float(p[index]))
        cape500.append(max(0.0, float(low.to("J/kg").magnitude)))
        capefull.append(max(0.0, float(full.to("J/kg").magnitude)))

    source_p = np.asarray(source_p)
    order = np.argsort(source_p)
    integrate = getattr(np, "trapezoid", None) or np.trapz
    icape500 = integrate(np.asarray(cape500)[order], source_p[order] * 100.0) / G
    icapefull = integrate(np.asarray(capefull)[order], source_p[order] * 100.0) / G
    return float(icape500), float(icapefull)


def process_case(row: dict, fields: dict, box_root: Path, tag: str) -> dict:
    result = {"case_id": row["case_id"], "error": ""}
    result.update({column: np.nan for column in OUTPUT_COLUMNS})
    try:
        file_2d = find_2d_file(row, box_root, tag)
        valid = noon_valid_mask(file_2d)
        cut = centred_cut(fields, row["lon"], row["lat"])
        if valid.shape != cut["temperature"].shape[1:]:
            raise ValueError(f"Mask {valid.shape} and raw cut {cut['temperature'].shape[1:]} differ")
        scale_errors = []
        for scale, distance in BOXES.items():
            try:
                profile = area_mean_profile(cut, valid, distance)
                low, full = cape_pair_for_sources(*profile)
                result[f"ICAPE500_{scale}"] = low
                result[f"ICAPEfull_{scale}"] = full
            except Exception as exc:
                scale_errors.append(f"{scale}: {type(exc).__name__}: {exc}")
        missing = [column for column in OUTPUT_COLUMNS if not np.isfinite(result[column])]
        if scale_errors or missing:
            details = "; ".join(scale_errors)
            result["error"] = f"{details}; missing={','.join(missing)}".strip("; ")
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
    return result


def process_day(task: DayTask) -> dict:
    output = Path(task.daily_file)
    if output.exists() and not task.overwrite:
        old = pd.read_csv(output)
        return {"tag": task.tag, "date": task.date, "cases": len(old), "seconds": 0.0, "status": "existing"}
    start = time.time()
    try:
        fields = load_day_fields(Path(task.raw_root), task.tag, task.date, task.surface_time_mode)
        results = [process_case(row, fields, Path(task.box_root), task.tag) for row in task.rows]
    except Exception as exc:
        message = f"DAY FAILURE {type(exc).__name__}: {exc}"
        results = [{"case_id": row["case_id"], "error": message,
                    **{column: np.nan for column in OUTPUT_COLUMNS}} for row in task.rows]
    frame = pd.DataFrame(results)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + f".tmp.{os.getpid()}")
    frame.to_csv(temporary, index=False)
    os.replace(temporary, output)
    return {"tag": task.tag, "date": task.date, "cases": len(frame),
            "seconds": time.time() - start, "status": "processed",
            "failures": int(frame["error"].astype(bool).sum())}


def save_enriched_table(original: Path, frame: pd.DataFrame, transposed: bool,
                        daily_files: list[Path], output: Path, overwrite: bool) -> None:
    if output.exists() and not overwrite:
        raise FileExistsError(f"Output exists: {output}; use --overwrite-output to replace it")
    results = pd.concat([pd.read_csv(path) for path in daily_files], ignore_index=True)
    if results["case_id"].duplicated().any() or set(results["case_id"]) != set(range(len(frame))):
        raise ValueError("Daily results do not provide exactly one row for every table case")
    results = results.set_index("case_id").sort_index()
    enriched = frame.copy()
    for column in OUTPUT_COLUMNS:
        enriched[column] = results[column].to_numpy()
    output.parent.mkdir(parents=True, exist_ok=True)
    (enriched.T if transposed else enriched).to_csv(output)

    failures = results.loc[results["error"].fillna("").astype(bool), ["error", *OUTPUT_COLUMNS]]
    failure_path = output.with_name(output.stem + "_failures.csv")
    if len(failures):
        failures.to_csv(failure_path)
    elif failure_path.exists():
        failure_path.unlink()
    print(f"Saved {output} ({len(enriched)} cases; {len(failures)} failures)")


def run_tag(args, tag: str, explicit_table: Path | None) -> None:
    table = resolve_table(args.table_dir, tag, explicit_table)
    frame, transposed = read_table(table)
    records = table_records(frame)
    grouped: dict[str, list[dict]] = {}
    for row in records:
        grouped.setdefault(row["date"], []).append(row)
    daily_dir = args.work_dir / table.stem
    daily_dir.mkdir(parents=True, exist_ok=True)
    tasks = [DayTask(tag, date, rows, str(args.raw_root), str(args.box_root),
                     str(daily_dir / f"{date}.csv"), args.surface_time_mode, args.overwrite_dates)
             for date, rows in sorted(grouped.items())]
    def needs_run(task):
        path = Path(task.daily_file)
        if task.overwrite or not path.exists():
            return True
        if args.retry_failed_dates:
            checkpoint = pd.read_csv(path, usecols=["error"])
            failed = checkpoint["error"].fillna("").astype(bool).any()
            if failed:
                task.overwrite = True
            return failed
        return False

    missing_tasks = [task for task in tasks if needs_run(task)]
    tasks_to_run = missing_tasks[:args.benchmark_dates] if args.benchmark_dates else missing_tasks

    start = time.time()
    summaries = []
    if args.processes == 1:
        for task in tasks_to_run:
            summary = process_day(task)
            summaries.append(summary)
            print(summary)
    else:
        with multiprocessing.Pool(args.processes) as pool:
            for summary in pool.imap_unordered(process_day, tasks_to_run):
                summaries.append(summary)
                print(summary)
    elapsed = time.time() - start
    processed_cases = sum(s["cases"] for s in summaries if s["status"] == "processed")
    if processed_cases:
        print(f"{tag}: processed {processed_cases} cases in {elapsed/60:.1f} min "
              f"({elapsed/processed_cases:.3f} wall-s per case with {args.processes} worker(s))")

    daily_files = [Path(task.daily_file) for task in tasks]
    missing = [path for path in daily_files if not path.exists()]
    if missing:
        print(f"{tag}: {len(missing)} date checkpoints still missing; table merge deferred.")
        return
    output = (args.hist_output if tag == "hist" else args.fut_output)
    output = output or table.with_name(table.stem + "_ICAPE.csv")
    save_enriched_table(table, frame, transposed, daily_files, output, args.overwrite_output)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tags", nargs="+", choices=("hist", "fut"), default=("hist", "fut"))
    parser.add_argument("--processes", type=int, default=4)
    parser.add_argument("--raw-root", type=Path, default=RAW_ROOT)
    parser.add_argument("--box-root", type=Path, default=BOX_ROOT)
    parser.add_argument("--table-dir", type=Path, default=TABLE_DIR)
    parser.add_argument("--work-dir", type=Path, default=TABLE_DIR / "ICAPE_daily")
    parser.add_argument("--hist-table", type=Path)
    parser.add_argument("--fut-table", type=Path)
    parser.add_argument("--hist-output", type=Path)
    parser.add_argument("--fut-output", type=Path)
    parser.add_argument("--surface-time-mode", choices=("instant", "3hmean"), default="instant",
                        help="Use instantaneous 12 UTC surface fields (recommended) or their 10-12 UTC mean.")
    parser.add_argument("--benchmark-dates", type=int, default=0,
                        help="Process only this many missing dates; checkpoints are reused by the full run.")
    parser.add_argument("--overwrite-dates", action="store_true")
    parser.add_argument("--retry-failed-dates", action="store_true",
                        help="Recalculate only date checkpoints containing failed cases.")
    parser.add_argument("--overwrite-output", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.processes < 1:
        raise ValueError("--processes must be at least 1")
    for tag in args.tags:
        run_tag(args, tag, args.hist_table if tag == "hist" else args.fut_table)


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
