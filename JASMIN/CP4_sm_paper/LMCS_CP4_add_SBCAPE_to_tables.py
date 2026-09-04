#!/usr/bin/env python3
"""Add surface-based CAPE to the final CP4 MCS mean tables.

This is the surface-parcel companion to LMCS_CP4_add_ICAPE_to_tables.py.  It
uses the same cases, raw 12-UTC fields, fixed spatial boxes, area-mean
soundings, and storm-day no-noon-rain mask.  The only physical difference is
that one parcel is lifted: the parcel defined by area-mean t2, q2 and p_srfc.

For every 0.25deg, 1deg, 2deg and fullBox sounding, two diagnostics are added:

  SBCAPE500  surface-based CAPE with the buoyancy calculation ending at 500 hPa
  SBCAPEfull surface-based CAPE through all available pressure levels

Both are conventional CAPE in J kg-1.  Existing ICAPE-enriched tables are
preferred automatically, so the default output normally contains both ICAPE
and SBCAPE.  Date checkpoints make the calculation safely resumable.
"""

from __future__ import annotations

import argparse
import multiprocessing
import os
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from metpy import calc
from metpy.units import units

import LMCS_CP4_add_ICAPE_to_tables as common


TABLE_DIR = common.TABLE_DIR
RAW_ROOT = common.RAW_ROOT
BOX_ROOT = common.BOX_ROOT
BOXES = common.BOXES
PARTIAL_TOP_HPA = common.PARTIAL_TOP_HPA
HALF_FULL_BOX = common.HALF_FULL_BOX
MIN_VALID_FRACTION = common.MIN_VALID_FRACTION
OUTPUT_COLUMNS = [f"{kind}_{scale}" for kind in ("SBCAPE500", "SBCAPEfull") for scale in BOXES]


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


def _cape_scalar(value) -> float:
    """Return a finite non-negative CAPE scalar without assuming NumPy scalar output."""
    result = float(np.asarray(value.to("J/kg").magnitude))
    return max(0.0, result) if np.isfinite(result) else np.nan


def _finite_mean(values, axis=None):
    """NaN-aware mean without RuntimeWarning for an entirely empty pressure level."""
    values = np.asarray(values, dtype=float)
    count = np.isfinite(values).sum(axis=axis)
    total = np.nansum(values, axis=axis)
    return np.divide(total, count, out=np.full(np.shape(total), np.nan, dtype=float), where=count > 0)


def area_mean_profile(cut: dict, valid: np.ndarray, distance: int | None) -> tuple:
    """Apply the same 50% rain-mask rule while retaining an explicit surface level."""
    size = valid.shape[0]
    if valid.shape != (size, size) or size != HALF_FULL_BOX * 2 + 1:
        raise ValueError(f"2-D mask/raw cut shape mismatch: {valid.shape}")
    centre = size // 2
    sl = slice(None) if distance is None else slice(centre - distance, centre + distance + 1)
    mask = valid[sl, sl]
    if mask.sum() < MIN_VALID_FRACTION * mask.size:
        return (None,) * 5

    t = np.where(mask[None, :, :], cut["temperature"][:, sl, sl], np.nan)
    q = np.where(mask[None, :, :], cut["humidity"][:, sl, sl], np.nan)
    t2 = _finite_mean(np.where(mask, cut["t2"][sl, sl], np.nan))
    q2 = _finite_mean(np.where(mask, cut["q2"][sl, sl], np.nan))
    ps = _finite_mean(np.where(mask, cut["ps"][sl, sl], np.nan))
    return cut["pressure"], _finite_mean(t, axis=(1, 2)), _finite_mean(q, axis=(1, 2)), \
        float(ps), (float(t2), float(q2))


def surface_cape_pair(pressure, temperature, humidity, ps, surface_state) -> tuple[float, float]:
    """Calculate 500-hPa-limited and full CAPE for the area-mean surface parcel."""
    if pressure is None or not np.all(np.isfinite([ps, *surface_state])):
        return np.nan, np.nan

    p, t, q = map(lambda values: np.asarray(values, dtype=float),
                  (pressure, temperature, humidity))
    # Pressure-level values at/below the mean surface are not part of the
    # atmospheric sounding.  t2/q2/p_srfc supply the always-above-ground base.
    good = np.isfinite(p) & np.isfinite(t) & np.isfinite(q) & (p < ps - 0.05)
    p, t, q = p[good], t[good], q[good]
    if p.size < 5:
        return np.nan, np.nan

    order = np.argsort(p)[::-1]
    p, t, q = p[order], t[order], q[order]
    p = np.concatenate(([ps], p))
    t = np.concatenate(([surface_state[0]], t))
    q = np.concatenate(([surface_state[1]], q))

    p_units = p * units.hPa
    t_units = t * units.K
    q_units = q * units("kg/kg")
    td = calc.dewpoint_from_specific_humidity(p_units, t_units, q_units).to("K")
    # Numerical/interpolation noise must not produce supersaturation.
    td = np.minimum(np.asarray(td.magnitude, dtype=float), t) * units.K
    parcel = calc.parcel_profile(p_units, t_units[0], td[0]).to("K")

    full, _ = calc.cape_cin(p_units, t_units, td, parcel,
                            which_lfc="bottom", which_el="top")
    full_value = _cape_scalar(full)

    partial = p_units >= PARTIAL_TOP_HPA * units.hPa
    if int(np.count_nonzero(partial)) < 3 or float(np.nanmin(p)) > PARTIAL_TOP_HPA:
        limited_value = np.nan
    else:
        limited, _ = calc.cape_cin(p_units[partial], t_units[partial], td[partial],
                                   parcel[partial], which_lfc="bottom", which_el="top")
        limited_value = _cape_scalar(limited)
    return limited_value, full_value


def process_case(row: dict, fields: dict, box_root: Path, tag: str) -> dict:
    result = {"case_id": row["case_id"], "error": ""}
    result.update({column: np.nan for column in OUTPUT_COLUMNS})
    try:
        file_2d = common.find_2d_file(row, box_root, tag)
        valid = common.noon_valid_mask(file_2d)
        cut = common.centred_cut(fields, row["lon"], row["lat"])
        if valid.shape != cut["temperature"].shape[1:]:
            raise ValueError(f"Mask {valid.shape} and raw cut {cut['temperature'].shape[1:]} differ")
        scale_errors = []
        for scale, distance in BOXES.items():
            try:
                profile = area_mean_profile(cut, valid, distance)
                limited, full = surface_cape_pair(*profile)
                result[f"SBCAPE500_{scale}"] = limited
                result[f"SBCAPEfull_{scale}"] = full
            except Exception as exc:
                scale_errors.append(f"{scale}: {type(exc).__name__}: {exc}")
        missing = [column for column in OUTPUT_COLUMNS if not np.isfinite(result[column])]
        if scale_errors or missing:
            detail = "; ".join(scale_errors)
            result["error"] = f"{detail}; missing={','.join(missing)}".strip("; ")
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
    return result


def process_day(task: DayTask) -> dict:
    output = Path(task.daily_file)
    if output.exists() and not task.overwrite:
        old = pd.read_csv(output)
        if set(OUTPUT_COLUMNS).issubset(old.columns):
            return {"tag": task.tag, "date": task.date, "cases": len(old),
                    "seconds": 0.0, "status": "existing"}
    start = time.time()
    try:
        fields = common.load_day_fields(Path(task.raw_root), task.tag, task.date,
                                        task.surface_time_mode)
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
            "failures": int(frame["error"].fillna("").astype(bool).sum())}


def resolve_input_table(table_dir: Path, tag: str, explicit: Path | None) -> Path:
    if explicit is not None:
        if not explicit.exists():
            raise FileNotFoundError(explicit)
        return explicit

    base = table_dir / f"{tag}_mean_table_JASMIN_3hmeansVersion_rainMask_fullBox.csv"
    icape = base.with_name(base.stem + "_ICAPE.csv")
    if icape.exists():
        return icape
    if base.exists():
        return base

    icape_candidates = sorted(path for path in table_dir.glob(
        f"{tag}_mean_table_JASMIN_3hmeansVersion_rainMask_fullBox*_ICAPE.csv"
    ) if "failures" not in path.stem and "SBCAPE" not in path.stem)
    if len(icape_candidates) == 1:
        return icape_candidates[0]
    base_candidates = sorted(path for path in table_dir.glob(
        f"{tag}_mean_table_JASMIN_3hmeansVersion_rainMask_fullBox*.csv"
    ) if not any(word in path.stem for word in ("ICAPE", "SBCAPE", "failures")))
    if len(base_candidates) == 1:
        return base_candidates[0]
    raise FileNotFoundError(
        f"Could not choose one {tag} mean input table; ICAPE={icape_candidates}, base={base_candidates}"
    )


def save_enriched_table(frame: pd.DataFrame, transposed: bool, daily_files: list[Path],
                        output: Path, overwrite: bool) -> None:
    if output.exists() and not overwrite:
        raise FileExistsError(f"Output exists: {output}; use --overwrite-output to replace it")
    results = pd.concat([pd.read_csv(path) for path in daily_files], ignore_index=True)
    expected_ids = set(range(len(frame)))
    if results["case_id"].duplicated().any() or set(results["case_id"]) != expected_ids:
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
    table = resolve_input_table(args.table_dir, tag, explicit_table)
    frame, transposed = common.read_table(table)
    records = common.table_records(frame)
    grouped: dict[str, list[dict]] = {}
    for row in records:
        grouped.setdefault(row["date"], []).append(row)

    daily_dir = args.work_dir / table.stem
    daily_dir.mkdir(parents=True, exist_ok=True)
    tasks = [DayTask(tag, date, rows, str(args.raw_root), str(args.box_root),
                     str(daily_dir / f"{date}.csv"), args.surface_time_mode,
                     args.overwrite_dates)
             for date, rows in sorted(grouped.items())]

    def needs_run(task: DayTask) -> bool:
        path = Path(task.daily_file)
        if task.overwrite or not path.exists():
            return True
        checkpoint = pd.read_csv(path)
        if not set(OUTPUT_COLUMNS).issubset(checkpoint.columns):
            task.overwrite = True
            return True
        if args.retry_failed_dates and checkpoint["error"].fillna("").astype(bool).any():
            task.overwrite = True
            return True
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
    processed_cases = sum(item["cases"] for item in summaries if item["status"] == "processed")
    if processed_cases:
        print(f"{tag}: processed {processed_cases} cases in {elapsed / 60:.1f} min "
              f"({elapsed / processed_cases:.3f} wall-s per case with {args.processes} worker(s))")

    daily_files = [Path(task.daily_file) for task in tasks]
    missing = [path for path in daily_files if not path.exists()]
    if missing:
        print(f"{tag}: {len(missing)} date checkpoints still missing; table merge deferred.")
        return
    output = args.hist_output if tag == "hist" else args.fut_output
    output = output or table.with_name(table.stem + "_SBCAPE.csv")
    save_enriched_table(frame, transposed, daily_files, output, args.overwrite_output)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tags", nargs="+", choices=("hist", "fut"), default=("hist", "fut"))
    parser.add_argument("--processes", type=int, default=4)
    parser.add_argument("--raw-root", type=Path, default=RAW_ROOT)
    parser.add_argument("--box-root", type=Path, default=BOX_ROOT)
    parser.add_argument("--table-dir", type=Path, default=TABLE_DIR)
    parser.add_argument("--work-dir", type=Path, default=TABLE_DIR / "SBCAPE_daily")
    parser.add_argument("--hist-table", type=Path)
    parser.add_argument("--fut-table", type=Path)
    parser.add_argument("--hist-output", type=Path)
    parser.add_argument("--fut-output", type=Path)
    parser.add_argument("--surface-time-mode", choices=("instant", "3hmean"), default="instant",
                        help="Must match the surface sampling used for the ICAPE calculation.")
    parser.add_argument("--benchmark-dates", type=int, default=0,
                        help="Process only this many missing dates; full-table merge is deferred.")
    parser.add_argument("--overwrite-dates", action="store_true")
    parser.add_argument("--retry-failed-dates", action="store_true")
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
