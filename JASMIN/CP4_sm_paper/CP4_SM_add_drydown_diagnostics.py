#!/usr/bin/env python3

"""
Add detailed drydown diagnostics to all_wetting_drydown_events_*.csv.

For every SM-defined and rain-defined wetting event:

  * locate the post-event SM maximum within event centre ... +5 h
  * determine the start of the next qualifying event episode
  * calculate the available uninterrupted drying interval
  * calculate fixed-time drying after:
        12, 24, 36, 48, 72, 96, 120, 168, 240 h
  * for each horizon save:
        dSM_Xh         = SM_peak - SM exactly X h later
        max_dSM_Xh     = SM_peak - minimum SM within first X h
        drying_slope_Xh = linear SM decline rate over first X h
        valid_Xh       = whether next event occurs after X h
  * calculate:
        minimum SM before next event
        total maximum drydown
        time peak -> minimum
        average peak -> minimum drying rate

The resulting table is intended to separate:

    available drying time
        versus
    drying response at fixed duration
        versus
    total realised drydown.
"""

from pathlib import Path
import argparse
import glob
import sys

import numpy as np
import pandas as pd
import xarray as xr

REPO_ROOT = Path("/home/users/cornkle/pythonWorkspace/proj_CEH")
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from shared.paths import load_paths

paths = load_paths()

ROOT = Path(paths.cp4_2d) / "CP4_box_JASMIN/mean3h_v2/timeseries"

MONTH = 8
DAYS = 30
BOX_PIXELS = 23
POST_HOURS = 5

# Includes 5 and 10 days explicitly.
FIXED_HOURS = [12, 24, 36, 48, 72, 96, 120, 168, 240]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--event-file", default=None)
    p.add_argument("--month", type=int, default=MONTH)
    p.add_argument("--days", type=int, default=DAYS)
    p.add_argument("--box-pixels", type=int, default=BOX_PIXELS)
    p.add_argument("--post-hours", type=float, default=POST_HOURS)
    return p.parse_args()


def find_latest_event_file():
    files = sorted(glob.glob(str(ROOT / "all_wetting_drydown_events_*.csv")))

    # Do not accidentally use an already augmented table.
    files = [
        f for f in files
        if "_detailed_diagnostics" not in Path(f).stem
    ]

    if not files:
        raise FileNotFoundError(
            f"No all_wetting_drydown_events_*.csv found in {ROOT}"
        )

    return Path(files[-1])


def timeseries_path(label, hour, filename, month, days, box_pixels):
    ftag = "hist" if label == "hist" else "fut"

    folder = ROOT / (
        f"SM_sh_lsRain_Aug{month:02d}_{ftag}_{int(hour):02d}UTC_"
        f"pm{days}d_box{box_pixels:02d}px"
    )

    return folder / filename


def get_sm(ds):
    sm = np.asarray(ds["SM"].values, float).squeeze()

    if "lag_hour_SM" in ds["SM"].coords:
        lag = np.asarray(
            ds["SM"]["lag_hour_SM"].values, float
        ).squeeze()
    else:
        lag = np.asarray(
            ds["lag_hour_SM"].values, float
        ).squeeze()

    return sm, lag


def nearest_idx(lag, target, tolerance=0.75):
    i = int(np.nanargmin(np.abs(lag - target)))

    if abs(lag[i] - target) > tolerance:
        return None

    return i


def post_event_peak(sm, lag, event_lag, post_hours=5):
    """
    Post-wetting SM maximum between event centre and +post_hours.
    """

    mask = (
        (lag >= event_lag) &
        (lag <= event_lag + post_hours) &
        np.isfinite(sm)
    )

    ii = np.where(mask)[0]

    if len(ii) == 0:
        return None, np.nan, np.nan

    idx = int(ii[np.nanargmax(sm[ii])])

    return idx, float(lag[idx]), float(sm[idx])


def fixed_horizon_metrics(
    sm,
    lag,
    peak_lag,
    peak_sm,
    horizon_h,
    next_event_start_lag,
):
    """
    Drying over an exactly fixed amount of available time.

    A horizon is valid only if the next qualifying event starts AFTER
    the target time.
    """

    target_lag = peak_lag + horizon_h

    if (
        np.isfinite(next_event_start_lag) and
        next_event_start_lag <= target_lag
    ):
        return False, np.nan, np.nan, np.nan

    target_idx = nearest_idx(lag, target_lag)

    if target_idx is None or not np.isfinite(sm[target_idx]):
        return False, np.nan, np.nan, np.nan

    # SM loss at exactly X hours.
    dsm = peak_sm - sm[target_idx]

    # Largest SM loss achieved at any point during first X hours.
    mask = (
        (lag >= peak_lag) &
        (lag <= target_lag) &
        np.isfinite(sm)
    )

    ii = np.where(mask)[0]

    if len(ii) == 0:
        return False, np.nan, np.nan, np.nan

    max_dsm = peak_sm - np.nanmin(sm[ii])

    # Fitted rate: positive values mean drying.
    if len(ii) >= 3:
        x = (lag[ii] - peak_lag) / 24.0
        y = sm[ii]

        slope_sm = np.polyfit(x, y, 1)[0]
        drying_slope = -slope_sm
    else:
        drying_slope = np.nan

    return (
        True,
        float(dsm),
        float(max_dsm),
        float(drying_slope),
    )


def drying_to_minimum(
    sm,
    lag,
    peak_lag,
    peak_sm,
    next_event_start_lag,
):
    """
    Lowest SM reached after the post-wetting peak and before the next
    qualifying event starts.
    """

    if np.isfinite(next_event_start_lag):
        mask = (
            (lag >= peak_lag) &
            (lag < next_event_start_lag) &
            np.isfinite(sm)
        )
    else:
        mask = (
            (lag >= peak_lag) &
            np.isfinite(sm)
        )

    ii = np.where(mask)[0]

    if len(ii) == 0:
        return (np.nan,) * 6

    min_idx = int(ii[np.nanargmin(sm[ii])])

    min_sm = float(sm[min_idx])
    min_lag = float(lag[min_idx])

    dsm = peak_sm - min_sm
    hours = min_lag - peak_lag
    days = hours / 24.0

    rate = dsm / days if days > 0 else np.nan

    return (
        min_sm,
        min_lag,
        float(dsm),
        float(hours),
        float(days),
        float(rate),
    )


def main():

    args = parse_args()

    event_file = (
        Path(args.event_file)
        if args.event_file is not None
        else find_latest_event_file()
    )

    print("Reading:")
    print(event_file)

    events = pd.read_csv(event_file)

    # --------------------------------------------------------
    # IMPORTANT:
    # Determine the actual start of the next event episode.
    #
    # Do this independently for SM-defined and rain-defined
    # catalogues.
    # --------------------------------------------------------

    events = events.sort_values(
        ["label", "hour", "file", "definition", "event_lag_h"]
    ).reset_index(drop=True)

    events["next_event_start_lag_h"] = (
        events.groupby(
            ["label", "hour", "file", "definition"],
            sort=False
        )["event_group_start_lag_h"]
        .shift(-1)
    )

    # Actual uninterrupted drying time will be calculated from
    # POST-EVENT SM PEAK to next-event start.
    events["available_drying_hours"] = np.nan
    events["available_drying_days"] = np.nan

    # Peak / ultimate minimum diagnostics.
    for col in [
        "SM_post_peak",
        "post_peak_lag_h",
        "SM_min_before_next",
        "min_SM_lag_h",
        "dSM_to_min",
        "hours_peak_to_min",
        "days_peak_to_min",
        "drying_rate_to_min",
    ]:
        events[col] = np.nan

    # Fixed-duration diagnostics.
    for h in FIXED_HOURS:
        events[f"valid_{h}h"] = False
        events[f"dSM_{h}h"] = np.nan
        events[f"max_dSM_{h}h"] = np.nan
        events[f"drying_slope_{h}h"] = np.nan

    # --------------------------------------------------------
    # Open every source NetCDF once.
    # --------------------------------------------------------

    groups = events.groupby(
        ["label", "hour", "file"],
        sort=False
    )

    print(f"Unique source files: {groups.ngroups}")
    print(f"Event rows:          {len(events)}")

    for ng, ((label, hour, filename), group) in enumerate(groups):

        if ng % 500 == 0:
            print(f"{ng}/{groups.ngroups}")

        path = timeseries_path(
            label,
            hour,
            filename,
            args.month,
            args.days,
            args.box_pixels,
        )

        if not path.exists():
            print("Missing:", path)
            continue

        try:
            with xr.open_dataset(path) as ds:
                sm, lag = get_sm(ds)

        except Exception as exc:
            print(f"FAILED {path.name}: {exc}")
            continue

        for row_idx, row in group.iterrows():

            event_lag = row["event_lag_h"]
            next_start = row["next_event_start_lag_h"]

            if not np.isfinite(event_lag):
                continue

            _, peak_lag, peak_sm = post_event_peak(
                sm,
                lag,
                event_lag,
                post_hours=args.post_hours,
            )

            if not np.isfinite(peak_lag):
                continue

            events.loc[row_idx, "SM_post_peak"] = peak_sm
            events.loc[row_idx, "post_peak_lag_h"] = peak_lag

            # ------------------------------------------------
            # Actual time available for uninterrupted drying.
            # ------------------------------------------------

            if np.isfinite(next_start):

                available_h = next_start - peak_lag

                if available_h > 0:
                    events.loc[
                        row_idx,
                        "available_drying_hours"
                    ] = available_h

                    events.loc[
                        row_idx,
                        "available_drying_days"
                    ] = available_h / 24.0

            # ------------------------------------------------
            # Fixed-time diagnostics.
            # ------------------------------------------------

            for h in FIXED_HOURS:

                valid, dsm, max_dsm, slope = (
                    fixed_horizon_metrics(
                        sm,
                        lag,
                        peak_lag,
                        peak_sm,
                        h,
                        next_start,
                    )
                )

                events.loc[row_idx, f"valid_{h}h"] = valid
                events.loc[row_idx, f"dSM_{h}h"] = dsm
                events.loc[row_idx, f"max_dSM_{h}h"] = max_dsm
                events.loc[row_idx, f"drying_slope_{h}h"] = slope

            # ------------------------------------------------
            # Total realised drying before next event.
            #
            # Do NOT use cases with no following event for the
            # duration decomposition because they are censored.
            # ------------------------------------------------

            if np.isfinite(next_start):

                (
                    min_sm,
                    min_lag,
                    dsm_min,
                    hours_min,
                    days_min,
                    rate_min,
                ) = drying_to_minimum(
                    sm,
                    lag,
                    peak_lag,
                    peak_sm,
                    next_start,
                )

                events.loc[
                    row_idx,
                    "SM_min_before_next"
                ] = min_sm

                events.loc[
                    row_idx,
                    "min_SM_lag_h"
                ] = min_lag

                events.loc[
                    row_idx,
                    "dSM_to_min"
                ] = dsm_min

                events.loc[
                    row_idx,
                    "hours_peak_to_min"
                ] = hours_min

                events.loc[
                    row_idx,
                    "days_peak_to_min"
                ] = days_min

                events.loc[
                    row_idx,
                    "drying_rate_to_min"
                ] = rate_min

    # --------------------------------------------------------
    # Save.
    # --------------------------------------------------------

    output = event_file.with_name(
        event_file.stem +
        "_detailed_diagnostics.csv"
    )

    events.to_csv(output, index=False)

    print("\nWrote:")
    print(output)

    # --------------------------------------------------------
    # Quick diagnostics.
    # --------------------------------------------------------

    print("\nMedian available drying days:")
    print(
        events.groupby(
            ["label", "definition"]
        )["available_drying_days"].median()
    )

    for h in [48, 120, 240]:

        print(f"\n===== {h} h =====")

        print("Valid samples:")
        print(
            events.groupby(
                ["label", "definition"]
            )[f"valid_{h}h"].sum()
        )

        print("\nMedian fixed-time decline:")
        print(
            events.groupby(
                ["label", "definition"]
            )[f"dSM_{h}h"].median()
        )

        print("\nMean fitted drying slope:")
        print(
            events.groupby(
                ["label", "definition"]
            )[f"drying_slope_{h}h"].mean()
        )

    print("\nMedian total drydown:")
    print(
        events.groupby(
            ["label", "definition"]
        )["dSM_to_min"].median()
    )


if __name__ == "__main__":
    main()