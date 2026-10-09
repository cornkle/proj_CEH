#!/usr/bin/env python3
"""
Analyse wetting-event return times and subsequent soil drydown around August CP4 storms.

Two independent wetting-event definitions are used:

1. SM-defined wetting:
       deltaSM >= 2.5 kg m-2

   where:
       SM_before = minimum SM during t-5 ... t
       SM_after  = maximum SM during t ... t+5
       deltaSM   = SM_after - SM_before

2. Rain-defined wetting:
       accumulated rainfall >= 5 mm

   where rainfall is summed during:
       t-5 ... t+5 h

Candidate windows belonging to the same episode are merged using a 12 h
minimum event-separation criterion.

After every wetting event, the script checks whether SM subsequently falls
by at least 2.5 kg m-2 before the next wetting event of the SAME definition.

Outputs
-------
1. One row per storm:
   - previous / t0 / next event for each definition
   - return intervals
   - t0 SM/rain characteristics
   - whether a >=2.5 kg m-2 drydown occurs after t0

2. One row per detected wetting event:
   - event amplitude
   - rainfall
   - next-event interval
   - maximum drydown before next event
   - whether >=2.5 kg m-2 drying is reached
   - time required to reach that drying threshold
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
from pathlib import Path

REPO_ROOT = Path("/home/users/cornkle/pythonWorkspace/proj_CEH")
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import numpy as np
import pandas as pd
import xarray as xr

from shared.paths import load_paths

paths = load_paths()

ROOT = Path(paths.cp4_2d) / "CP4_box_JASMIN/mean3h_v2/timeseries"


# ============================================================
# ARGUMENTS
# ============================================================

def parse_args():
    p = argparse.ArgumentParser(description=__doc__)

    p.add_argument("--hours", default="17,18,19")
    p.add_argument("--month", type=int, default=8)
    p.add_argument("--days", type=int, default=30)
    p.add_argument("--box-pixels", type=int, default=23)

    p.add_argument(
        "--sm-wetting-threshold",
        type=float,
        default=2.5,
        help="Minimum event-scale SM rise in kg m-2. Default 2.5.",
    )

    p.add_argument(
        "--rain-threshold",
        type=float,
        default=5.0,
        help="Minimum accumulated rainfall over t-5...t+5 h in mm. Default 5.",
    )

    p.add_argument(
        "--drying-threshold",
        type=float,
        default=2.5,
        help="Meaningful post-event SM decline in kg m-2. Default 2.5.",
    )

    p.add_argument(
        "--pre-hours",
        type=float,
        default=5,
    )

    p.add_argument(
        "--post-hours",
        type=float,
        default=5,
    )

    p.add_argument(
        "--event-gap-hours",
        type=float,
        default=12,
        help="Minimum dry/event-free gap separating distinct events. Default 12 h.",
    )

    p.add_argument(
        "--max-storms",
        type=int,
        default=None,
        help="For testing only.",
    )

    return p.parse_args()


def parse_hours(s):
    return [int(x) for x in s.split(",")]


# ============================================================
# FILE DISCOVERY
# ============================================================

def files_for(label, hour, month, days, box_pixels):
    ftag = "hist" if label == "hist" else "fut"

    folder = ROOT / (
        f"SM_sh_lsRain_Aug{month:02d}_{ftag}_{hour:02d}UTC_"
        f"pm{days}d_box{box_pixels:02d}px"
    )

    files = sorted(glob.glob(str(folder / "*_timeseries.nc")))

    print(f"{label} {hour}: {len(files)} files")

    return files


# ============================================================
# READ VARIABLES
# ============================================================

def get_var_lag(ds, var):
    vals = np.asarray(ds[var].values, float).squeeze()
    lag_name = f"lag_hour_{var}"

    if lag_name in ds[var].coords:
        lag = np.asarray(ds[var][lag_name].values, float).squeeze()
    elif lag_name in ds:
        lag = np.asarray(ds[lag_name].values, float).squeeze()
    else:
        raise KeyError(f"{lag_name} not found")

    if vals.ndim != 1 or lag.ndim != 1:
        raise ValueError(f"{var}: expected 1D arrays, got {vals.shape}, {lag.shape}")

    return vals, lag


def nearest_lag_index(lag, target, tolerance=0.75):
    idx = int(np.nanargmin(np.abs(lag - target)))
    return idx if abs(lag[idx] - target) <= tolerance else None


# ============================================================
# EVENT-WINDOW DIAGNOSTICS
# ============================================================

def event_metrics_at_index(
    sm,
    sm_lag,
    rain,
    rain_lag,
    idx,
    pre_hours=5,
    post_hours=5,
):
    """
    Calculate SM and rainfall metrics centred on sm_lag[idx].

    deltaSM:
        min SM in t-5...t
        -> max SM in t...t+5

    rain:
        accumulation in t-5...t+5
    """

    t = sm_lag[idx]

    pre = (
        (sm_lag >= t - pre_hours) &
        (sm_lag <= t)
    )

    post = (
        (sm_lag >= t) &
        (sm_lag <= t + post_hours)
    )

    rr = (
        (rain_lag >= t - pre_hours) &
        (rain_lag <= t + post_hours)
    )

    if pre.sum() == 0 or post.sum() == 0:
        return None

    pre_indices = np.where(pre)[0]
    post_indices = np.where(post)[0]

    if (
        np.all(~np.isfinite(sm[pre_indices])) or
        np.all(~np.isfinite(sm[post_indices]))
    ):
        return None

    pre_idx = int(
        pre_indices[np.nanargmin(sm[pre_indices])]
    )

    post_idx = int(
        post_indices[np.nanargmax(sm[post_indices])]
    )

    sm_before = float(sm[pre_idx])
    sm_after = float(sm[post_idx])

    rain_accum = (
        float(np.nansum(rain[rr]))
        if rr.sum() else np.nan
    )

    return {
        "centre_idx": int(idx),
        "centre_lag": float(t),

        "pre_idx": pre_idx,
        "post_idx": post_idx,

        "pre_lag": float(sm_lag[pre_idx]),
        "post_lag": float(sm_lag[post_idx]),

        "SM_before": sm_before,
        "SM_after": sm_after,

        "dSM": sm_after - sm_before,
        "rain_accum": rain_accum,
    }


def calculate_candidate_metrics(
    sm,
    sm_lag,
    rain,
    rain_lag,
    pre_hours,
    post_hours,
):
    """
    Calculate t-5...t+5 metrics for every possible candidate centre.
    """

    candidates = []

    for i in range(len(sm)):
        metric = event_metrics_at_index(
            sm,
            sm_lag,
            rain,
            rain_lag,
            i,
            pre_hours,
            post_hours,
        )

        if metric is not None:
            candidates.append(metric)

    return candidates


# ============================================================
# DISTINCT EVENT CATALOGUE
# ============================================================

def make_event_catalogue(
    candidates,
    definition,
    threshold,
    event_gap_hours=12,
):
    """
    Select qualifying candidate windows and merge candidates belonging
    to the same episode.

    definition:
        "sm"   -> dSM >= threshold
        "rain" -> rainfall accumulation >= threshold

    Representative event:
        largest dSM for SM-defined events
        largest rain accumulation for rain-defined events
    """

    if definition == "sm":
        qualifying = [
            c for c in candidates
            if np.isfinite(c["dSM"]) and
            c["dSM"] >= threshold
        ]
        ranking_key = "dSM"

    elif definition == "rain":
        qualifying = [
            c for c in candidates
            if np.isfinite(c["rain_accum"]) and
            c["rain_accum"] >= threshold
        ]
        ranking_key = "rain_accum"

    else:
        raise ValueError(definition)

    if len(qualifying) == 0:
        return []

    qualifying = sorted(
        qualifying,
        key=lambda x: x["centre_lag"]
    )

    groups = []
    current = [qualifying[0]]

    for c in qualifying[1:]:

        gap = (
            c["centre_lag"] -
            current[-1]["centre_lag"]
        )

        if gap <= event_gap_hours:
            current.append(c)

        else:
            groups.append(current)
            current = [c]

    groups.append(current)

    events = []

    for group in groups:

        representative = max(
            group,
            key=lambda x: x[ranking_key]
        ).copy()

        representative["definition"] = definition

        representative["group_start_lag"] = float(
            group[0]["centre_lag"]
        )

        representative["group_end_lag"] = float(
            group[-1]["centre_lag"]
        )

        events.append(representative)

    return events


# ============================================================
# DRYDOWN AFTER EACH EVENT
# ============================================================

def add_drydown_metrics(
    events,
    sm,
    sm_lag,
    drying_threshold=2.5,
):
    """
    For each event:

    Start from its post-event SM maximum.

    Search forward until the next event of the SAME definition.

    Calculate:
        maximum SM decline achieved before next event
        whether decline >= drying_threshold
        first time that drying threshold is reached
        interval to next wetting event
    """

    for i, event in enumerate(events):

        start_idx = event["post_idx"]
        start_sm = event["SM_after"]

        if i < len(events) - 1:

            next_event = events[i + 1]

            # Stop before the pre-event minimum/window belonging
            # to the next distinct wetting episode.
            end_idx = next_event["pre_idx"]

            event["next_event_lag"] = next_event["centre_lag"]

            event["hours_to_next_event"] = (
                next_event["centre_lag"] -
                event["centre_lag"]
            )

        else:

            end_idx = len(sm) - 1

            event["next_event_lag"] = np.nan
            event["hours_to_next_event"] = np.nan

        if end_idx <= start_idx:

            event["max_drydown"] = np.nan
            event["drydown_ge_threshold"] = False
            event["hours_to_dry_threshold"] = np.nan

            continue

        ii = np.arange(start_idx, end_idx + 1)

        valid = np.isfinite(sm[ii])

        if valid.sum() == 0:

            event["max_drydown"] = np.nan
            event["drydown_ge_threshold"] = False
            event["hours_to_dry_threshold"] = np.nan

            continue

        ii = ii[valid]

        minimum_sm = np.nanmin(sm[ii])

        max_drydown = (
            start_sm - minimum_sm
        )

        event["max_drydown"] = float(max_drydown)

        reached = (
            start_sm - sm[ii]
        ) >= drying_threshold

        if np.any(reached):

            dry_idx = int(ii[np.where(reached)[0][0]])

            event["drydown_ge_threshold"] = True

            event["hours_to_dry_threshold"] = float(
                sm_lag[dry_idx] -
                sm_lag[start_idx]
            )

            event["days_to_dry_threshold"] = (
                event["hours_to_dry_threshold"] / 24
            )

            event["dry_threshold_lag"] = float(
                sm_lag[dry_idx]
            )

        else:

            event["drydown_ge_threshold"] = False
            event["hours_to_dry_threshold"] = np.nan
            event["days_to_dry_threshold"] = np.nan
            event["dry_threshold_lag"] = np.nan

    return events


# ============================================================
# T0 / PREVIOUS / NEXT EVENT
# ============================================================

def classify_relative_to_t0(events, event_gap_hours=12):
    """
    Identify the event belonging to the t0 convective cycle,
    plus previous and next distinct event.

    First preference:
        event group actually overlaps lag=0

    Fallback:
        representative event within +/- event_gap_hours of t0.
    """

    if not events:
        return None, None, None

    t0_candidates = [
        (i, e)
        for i, e in enumerate(events)
        if (
            e["group_start_lag"] <= 0 <= e["group_end_lag"]
        )
    ]

    if not t0_candidates:

        t0_candidates = [
            (i, e)
            for i, e in enumerate(events)
            if abs(e["centre_lag"]) <= event_gap_hours
        ]

    if t0_candidates:

        i0, t0_event = min(
            t0_candidates,
            key=lambda x: abs(x[1]["centre_lag"])
        )

        prev_event = (
            events[i0 - 1]
            if i0 > 0 else None
        )

        next_event = (
            events[i0 + 1]
            if i0 < len(events) - 1 else None
        )

    else:

        t0_event = None

        previous = [
            e for e in events
            if e["centre_lag"] < 0
        ]

        following = [
            e for e in events
            if e["centre_lag"] > 0
        ]

        prev_event = (
            previous[-1]
            if previous else None
        )

        next_event = (
            following[0]
            if following else None
        )

    return prev_event, t0_event, next_event


# ============================================================
# EVENT -> OUTPUT ROW
# ============================================================

def event_to_row(
    event,
    label,
    hour,
    filename,
    definition,
):
    return {
        "label": label,
        "hour": hour,
        "file": filename,
        "definition": definition,

        "event_lag_h": event["centre_lag"],

        "event_group_start_lag_h":
            event["group_start_lag"],

        "event_group_end_lag_h":
            event["group_end_lag"],

        "SM_before": event["SM_before"],
        "SM_after": event["SM_after"],
        "dSM": event["dSM"],

        "rain_accum_m5_p5":
            event["rain_accum"],

        "max_drydown":
            event.get("max_drydown", np.nan),

        "drydown_ge_2p5":
            event.get(
                "drydown_ge_threshold",
                False
            ),

        "hours_to_dry_2p5":
            event.get(
                "hours_to_dry_threshold",
                np.nan
            ),

        "days_to_dry_2p5":
            event.get(
                "days_to_dry_threshold",
                np.nan
            ),

        "next_event_lag_h":
            event.get(
                "next_event_lag",
                np.nan
            ),

        "hours_to_next_event":
            event.get(
                "hours_to_next_event",
                np.nan
            ),
    }


# ============================================================
# STORM-CENTRED SUMMARY
# ============================================================

def add_relative_event_columns(
    out,
    prefix,
    prev_event,
    t0_event,
    next_event,
):
    """
    Add previous / t0 / next-event diagnostics using a prefix:
        sm_...
        rain_...
    """

    out[f"{prefix}_t0_qualifies"] = (
        t0_event is not None
    )

    if prev_event is not None:

        out[f"{prefix}_prev_event_lag_h"] = (
            prev_event["centre_lag"]
        )

        out[f"{prefix}_prev_event_dSM"] = (
            prev_event["dSM"]
        )

        out[f"{prefix}_prev_event_rain_mm"] = (
            prev_event["rain_accum"]
        )

    else:

        out[f"{prefix}_prev_event_lag_h"] = np.nan
        out[f"{prefix}_prev_event_dSM"] = np.nan
        out[f"{prefix}_prev_event_rain_mm"] = np.nan

    if t0_event is not None:

        out[f"{prefix}_t0_event_lag_h"] = (
            t0_event["centre_lag"]
        )

        out[f"{prefix}_t0_event_dSM"] = (
            t0_event["dSM"]
        )

        out[f"{prefix}_t0_event_rain_mm"] = (
            t0_event["rain_accum"]
        )

        out[f"{prefix}_t0_max_drydown"] = (
            t0_event.get(
                "max_drydown",
                np.nan
            )
        )

        out[f"{prefix}_t0_drydown_ge_2p5"] = (
            t0_event.get(
                "drydown_ge_threshold",
                False
            )
        )

        out[f"{prefix}_t0_hours_to_dry_2p5"] = (
            t0_event.get(
                "hours_to_dry_threshold",
                np.nan
            )
        )

    else:

        out[f"{prefix}_t0_event_lag_h"] = np.nan
        out[f"{prefix}_t0_event_dSM"] = np.nan
        out[f"{prefix}_t0_event_rain_mm"] = np.nan
        out[f"{prefix}_t0_max_drydown"] = np.nan
        out[f"{prefix}_t0_drydown_ge_2p5"] = False
        out[f"{prefix}_t0_hours_to_dry_2p5"] = np.nan

    if next_event is not None:

        out[f"{prefix}_next_event_lag_h"] = (
            next_event["centre_lag"]
        )

        out[f"{prefix}_next_event_dSM"] = (
            next_event["dSM"]
        )

        out[f"{prefix}_next_event_rain_mm"] = (
            next_event["rain_accum"]
        )

    else:

        out[f"{prefix}_next_event_lag_h"] = np.nan
        out[f"{prefix}_next_event_dSM"] = np.nan
        out[f"{prefix}_next_event_rain_mm"] = np.nan

    # Event-to-event return intervals.

    if prev_event is not None and t0_event is not None:

        out[f"{prefix}_days_prev_to_t0"] = (
            t0_event["centre_lag"] -
            prev_event["centre_lag"]
        ) / 24

    else:

        out[f"{prefix}_days_prev_to_t0"] = np.nan

    if t0_event is not None and next_event is not None:

        out[f"{prefix}_days_t0_to_next"] = (
            next_event["centre_lag"] -
            t0_event["centre_lag"]
        ) / 24

    else:

        out[f"{prefix}_days_t0_to_next"] = np.nan


# ============================================================
# ONE STORM FILE
# ============================================================

def analyse_file(
    path,
    label,
    hour,
    sm_threshold,
    rain_threshold,
    drying_threshold,
    pre_hours,
    post_hours,
    event_gap_hours,
):

    with xr.open_dataset(path) as ds:

        sm, sm_lag = get_var_lag(ds, "SM")
        rain, rain_lag = get_var_lag(ds, "lsRain")

        # ----------------------------------------------------
        # Metrics at actual storm t0, independent of threshold.
        # ----------------------------------------------------

        t0_idx = nearest_lag_index(sm_lag, 0)

        if t0_idx is None:
            raise ValueError("Could not find t0")

        actual_t0 = event_metrics_at_index(
            sm,
            sm_lag,
            rain,
            rain_lag,
            t0_idx,
            pre_hours,
            post_hours,
        )

        # ----------------------------------------------------
        # Calculate all moving-window candidates once.
        # ----------------------------------------------------

        candidates = calculate_candidate_metrics(
            sm,
            sm_lag,
            rain,
            rain_lag,
            pre_hours,
            post_hours,
        )

        # ----------------------------------------------------
        # SM-defined event catalogue.
        # ----------------------------------------------------

        sm_events = make_event_catalogue(
            candidates,
            definition="sm",
            threshold=sm_threshold,
            event_gap_hours=event_gap_hours,
        )

        sm_events = add_drydown_metrics(
            sm_events,
            sm,
            sm_lag,
            drying_threshold=drying_threshold,
        )

        # ----------------------------------------------------
        # Rain-defined event catalogue.
        # ----------------------------------------------------

        rain_events = make_event_catalogue(
            candidates,
            definition="rain",
            threshold=rain_threshold,
            event_gap_hours=event_gap_hours,
        )

        rain_events = add_drydown_metrics(
            rain_events,
            sm,
            sm_lag,
            drying_threshold=drying_threshold,
        )

        # ----------------------------------------------------
        # Identify previous / t0 / next events.
        # ----------------------------------------------------

        sm_prev, sm_t0, sm_next = (
            classify_relative_to_t0(
                sm_events,
                event_gap_hours,
            )
        )

        rain_prev, rain_t0, rain_next = (
            classify_relative_to_t0(
                rain_events,
                event_gap_hours,
            )
        )

        # ----------------------------------------------------
        # Storm-centred output.
        # ----------------------------------------------------

        out = {
            "label": label,
            "hour": hour,
            "file": os.path.basename(path),

            "storm_time":
                ds.attrs.get("storm_time", ""),

            "storm_lon":
                float(
                    ds.attrs.get(
                        "storm_lon",
                        np.nan
                    )
                ),

            "storm_lat":
                float(
                    ds.attrs.get(
                        "storm_lat",
                        np.nan
                    )
                ),

            # Actual t0 properties, regardless of whether
            # either threshold is met.

            "t0_SM_before":
                actual_t0["SM_before"],

            "t0_SM_after":
                actual_t0["SM_after"],

            "t0_dSM":
                actual_t0["dSM"],

            "t0_rain_m5_p5":
                actual_t0["rain_accum"],

            "n_SM_events":
                len(sm_events),

            "n_rain_events":
                len(rain_events),
        }

        add_relative_event_columns(
            out,
            "sm",
            sm_prev,
            sm_t0,
            sm_next,
        )

        add_relative_event_columns(
            out,
            "rain",
            rain_prev,
            rain_t0,
            rain_next,
        )

        # ----------------------------------------------------
        # Long-form event rows.
        # ----------------------------------------------------

        event_rows = []

        for event in sm_events:

            event_rows.append(
                event_to_row(
                    event,
                    label,
                    hour,
                    os.path.basename(path),
                    "SM>=2.5",
                )
            )

        for event in rain_events:

            event_rows.append(
                event_to_row(
                    event,
                    label,
                    hour,
                    os.path.basename(path),
                    "Rain>=5mm",
                )
            )

        return out, event_rows


# ============================================================
# MAIN
# ============================================================

def main():

    args = parse_args()
    hours = parse_hours(args.hours)

    storm_rows = []
    all_event_rows = []

    for label in ["hist", "future"]:

        for hour in hours:

            files = files_for(
                label,
                hour,
                args.month,
                args.days,
                args.box_pixels,
            )

            if args.max_storms is not None:
                files = files[:args.max_storms]

            for n, path in enumerate(files):

                if n % 500 == 0:
                    print(
                        f"{label} {hour}: "
                        f"{n}/{len(files)}"
                    )

                try:

                    storm_row, event_rows = (
                        analyse_file(
                            path,
                            label,
                            hour,
                            args.sm_wetting_threshold,
                            args.rain_threshold,
                            args.drying_threshold,
                            args.pre_hours,
                            args.post_hours,
                            args.event_gap_hours,
                        )
                    )

                    storm_rows.append(storm_row)
                    all_event_rows.extend(event_rows)

                except Exception as exc:

                    print(
                        f"FAILED {Path(path).name}: "
                        f"{exc}"
                    )

    storm_df = pd.DataFrame(storm_rows)
    event_df = pd.DataFrame(all_event_rows)

    hours_tag = "".join(
        str(h) for h in hours
    )

    sm_tag = str(
        args.sm_wetting_threshold
    ).replace(".", "p")

    rain_tag = str(
        args.rain_threshold
    ).replace(".", "p")

    dry_tag = str(
        args.drying_threshold
    ).replace(".", "p")

    prefix = (
        f"Aug{args.month:02d}_"
        f"hours{hours_tag}_"
        f"SMwet{sm_tag}_"
        f"Rain{rain_tag}mm_"
        f"Dry{dry_tag}_"
        f"gap{int(args.event_gap_hours)}h_"
        f"pm{args.days}d_"
        f"box{args.box_pixels:02d}px"
    )

    storm_file = ROOT / (
        f"storm_event_summary_{prefix}.csv"
    )

    event_file = ROOT / (
        f"all_wetting_drydown_events_{prefix}.csv"
    )

    storm_df.to_csv(
        storm_file,
        index=False
    )

    event_df.to_csv(
        event_file,
        index=False
    )

    print("\nWrote:")
    print(storm_file)
    print(event_file)

    # --------------------------------------------------------
    # Useful summary output
    # --------------------------------------------------------

    if len(event_df):

        print("\nNumber of events:")

        print(
            event_df.groupby(
                ["label", "definition"]
            ).size()
        )

        print(
            "\nFraction of wetting events followed "
            "by >=2.5 kg m-2 drying before next event:"
        )

        print(
            event_df.groupby(
                ["label", "definition"]
            )["drydown_ge_2p5"].mean()
        )

        print(
            "\nMedian days to >=2.5 kg m-2 drying:"
        )

        reached = event_df[
            event_df["drydown_ge_2p5"] == True
        ]

        print(
            reached.groupby(
                ["label", "definition"]
            )["days_to_dry_2p5"].median()
        )

        print(
            "\nMedian interval to next wetting event (days):"
        )

        tmp = event_df.copy()

        tmp["days_to_next_event"] = (
            tmp["hours_to_next_event"] / 24
        )

        print(
            tmp.groupby(
                ["label", "definition"]
            )["days_to_next_event"].median()
        )


if __name__ == "__main__":
    main()