#!/usr/bin/env python
"""Merged CP4 Figure 1: six spatial composites plus the change decomposition.

Notebook use
------------
Run ``%matplotlib widget`` in a cell before this code if you want interactive
zooming/panning.  The code itself is ordinary Python, so it can also be run as
a script.
"""

from pathlib import Path

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd
import xarray as xr
from metpy import calc
from metpy.units import units

from shared.paths import load_paths


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
paths = load_paths()
ROOT = Path(paths.cp4_2d) / "CP4_box_JASMIN" / "mean3h_v2"
TABLE_DIR = ROOT / "tables"
OUTPUT_FIG = Path("cp4_fig1_merged.png")

HOURS = (17, 18)
GRID_KM = 4.4
SAMPLE_BOX_KM = 200
DTAG = "2deg"
PMAX_MIN = 1.0

# Add a variable here only if you want its physically favourable direction to
# plot upward, e.g. SIGN_FLIP = {"div_2deg"} for stronger convergence.
SIGN_FLIP = set()

BAR_VARS = [
    "SM_2deg", "sh_2deg", "lh_2deg", "ef_2deg", "sw_net_2deg",
    "t_srfc_2deg", "q_srfc_2deg", "tcwv_2deg", "vpd_srfc_2deg", "div_2deg",
]
BAR_LABELS = {
    "SM_2deg": "SM", "sh_2deg": "SH", "lh_2deg": "LH",
    "ef_2deg": "EF", "sw_net_2deg": r"SW$_{net}$",
    "t_srfc_2deg": r"T$_{sfc}$", "q_srfc_2deg": r"q$_{sfc}$",
    "tcwv_2deg": "TCWV", "vpd_srfc_2deg": "VPD", "div_2deg": "DIV",
}

COLORS = {
    "blue": "#31688E", "orange": "#D97757", "ink": "#27313A",
    "muted": "#68737D", "sm": "#2468A2", "paper": "#FAFAF8",
}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9.5,
    "axes.titlesize": 11,
    "axes.labelsize": 9.5,
    "axes.titleweight": "semibold",
    "axes.edgecolor": "#6B737A",
    "axes.linewidth": 0.8,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8.5,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})


# -----------------------------------------------------------------------------
# Derived variables
# -----------------------------------------------------------------------------
def _as_float(values):
    return np.asarray(values, dtype=float)


def _normalise_t_q(t, q):
    """Return temperature in K and specific humidity in kg kg-1."""
    t, q = _as_float(t), _as_float(q)
    t_k = np.where(t < 150, t + 273.15, t)
    q_kgkg = np.where(q > 0.1, q / 1000.0, q)
    q_kgkg = np.where((q_kgkg > 0) & (q_kgkg < 1), q_kgkg, np.nan)
    return t_k, q_kgkg


def calc_vpd_from_t_q(t, q, p_hpa=925.0):
    """Vapour-pressure deficit (hPa) from T and specific humidity."""
    t_k, q = _normalise_t_q(t, q)
    t_c = t_k - 273.15
    es = 6.112 * np.exp(17.67 * t_c / (t_c + 243.5))
    e = q * p_hpa / (0.622 + 0.378 * q)
    return es - e


def calc_ef(lh, sh):
    """Evaporative fraction in the range implied by the input fluxes."""
    lh, sh = _as_float(lh), _as_float(sh)
    denom = lh + sh
    return np.divide(lh, denom, out=np.full_like(lh, np.nan),
                     where=np.isfinite(denom) & ~np.isclose(denom, 0))


def saturation_vapour_pressure_hpa(t_k):
    t_c = t_k - 273.15
    return 6.112 * np.exp(17.67 * t_c / (t_c + 243.5))


def calc_thetae_from_t_q(t, q, p_hpa=925.0):
    """Bolton-style equivalent potential temperature (K)."""
    t_k, q = _normalise_t_q(t, q)
    r = q / (1 - q)
    e = q * p_hpa / (0.622 + 0.378 * q)
    e = np.where(e > 0, e, np.nan)
    ln_ratio = np.log(e / 6.112)
    td_k = 243.5 * ln_ratio / (17.67 - ln_ratio) + 273.15
    tlcl = 1 / (1 / (td_k - 56) + np.log(t_k / td_k) / 800) + 56
    theta = t_k * (1000 / p_hpa) ** 0.2854
    return theta * np.exp((3376 / tlcl - 2.54) * r * (1 + 0.81 * r))


def calc_thetaes_from_t(t, p_hpa=650.0):
    """Saturated equivalent potential temperature (K)."""
    t = _as_float(t)
    t_k = np.where(t < 150, t + 273.15, t)
    es = saturation_vapour_pressure_hpa(t_k)
    es = np.where(es < 0.99 * p_hpa, es, np.nan)
    rs = 0.622 * es / (p_hpa - es)
    theta = t_k * (1000 / p_hpa) ** 0.2854
    return theta * np.exp((3376 / t_k - 2.54) * rs * (1 + 0.81 * rs))


def add_derived_variables(abs_df, anom_df, suffix=DTAG):
    """Add VPD, EF, theta-e and the two-level CAPE proxy to one climate."""
    t_s, q_s = f"t_srfc_{suffix}", f"q_srfc_{suffix}"
    lh, sh, t_m = f"lh_{suffix}", f"sh_{suffix}", f"t_mid_{suffix}"

    t_clim = abs_df[t_s] - anom_df[t_s]
    q_clim = abs_df[q_s] - anom_df[q_s]
    lh_clim = abs_df[lh] - anom_df[lh]
    sh_clim = abs_df[sh] - anom_df[sh]
    tm_clim = abs_df[t_m] - anom_df[t_m]

    derived = {
        f"vpd_srfc_{suffix}": (
            calc_vpd_from_t_q(abs_df[t_s], abs_df[q_s]),
            calc_vpd_from_t_q(t_clim, q_clim),
            1.0,
        ),
        f"ef_{suffix}": (
            calc_ef(abs_df[lh], abs_df[sh]), calc_ef(lh_clim, sh_clim), 100.0,
        ),
        f"thetae_{suffix}": (
            calc_thetae_from_t_q(abs_df[t_s], abs_df[q_s]),
            calc_thetae_from_t_q(t_clim, q_clim),
            1.0,
        ),
    }

    cape_abs = (calc_thetae_from_t_q(abs_df[t_s], abs_df[q_s])
                - calc_thetaes_from_t(abs_df[t_m]))
    cape_clim = calc_thetae_from_t_q(t_clim, q_clim) - calc_thetaes_from_t(tm_clim)
    derived[f"cape_proxy_calc_{suffix}"] = (cape_abs, cape_clim, 1.0)

    for name, (absolute, climatology, scale) in derived.items():
        abs_df[name] = absolute * scale
        anom_df[name] = (absolute - climatology) * scale
    return abs_df, anom_df


def calc_divergence(ds, grid_m=GRID_KM * 1000):
    """Horizontal divergence on the composite grid (s-1)."""
    u = units.Quantity(ds["u_srfc"].squeeze().values, "m/s")
    v = units.Quantity(ds["v_srfc"].squeeze().values, "m/s")
    return calc.divergence(u, v, dx=grid_m * units.m, dy=grid_m * units.m).magnitude


# -----------------------------------------------------------------------------
# Data loading
# -----------------------------------------------------------------------------
def mean_hours(period, kind):
    """Load and average the requested composite hours into memory."""
    datasets = []
    try:
        for hour in HOURS:
            filename = ROOT / f"{period}_{kind}_Jul-Sep_{hour}h_2d.nc"
            datasets.append(xr.open_dataset(filename))
        return xr.concat(datasets, dim="composite_hour").mean("composite_hour").load()
    finally:
        for ds in datasets:
            ds.close()


def read_numeric_table(filename):
    df = pd.read_csv(filename, index_col=0, low_memory=False).T
    for col in df.columns:
        try:
            df[col] = pd.to_numeric(df[col])
        except (TypeError, ValueError):
            pass
    return df


def align_absolute_and_anomaly(abs_df, anom_df):
    """Retain storm cases present in both tables, following the original logic."""
    keys = ["pmax", "area", "area_1mm"]
    common = pd.merge(anom_df[keys].drop_duplicates(),
                      abs_df[keys].drop_duplicates(), on=keys, how="inner")
    return (abs_df.merge(common, on=keys, how="inner"),
            anom_df.merge(common, on=keys, how="inner"))


def load_tables():
    mean_name = "{}_mean_table_JASMIN_3hmeansVersion_rainMask_fullBox_ICAPE.csv"
    anom_name = "{}_anom_table_JASMIN_3hmeansVersion_rainMask_fullBox.csv"
    output = {}
    for period in ("hist", "fut"):
        absolute = read_numeric_table(TABLE_DIR / mean_name.format(period))
        anomaly = read_numeric_table(TABLE_DIR / anom_name.format(period))
        absolute, anomaly = align_absolute_and_anomaly(absolute, anomaly)
        absolute, anomaly = add_derived_variables(absolute, anomaly)
        keep = absolute["pmax"] >= PMAX_MIN
        output[period] = (absolute.loc[keep].copy(), anomaly.loc[keep].copy())
        print(f"{period}: retained {keep.sum():,} matched cases")
    return output


def calculate_bar_contributions(lhist, lahist, lfut, lafut):
    """Decompose total storm-day change into climate and anomaly components."""
    rows = []
    for var in BAR_VARS:
        storm_hist, storm_fut = np.nanmean(lhist[var]), np.nanmean(lfut[var])
        anom_hist, anom_fut = np.nanmean(lahist[var]), np.nanmean(lafut[var])
        clim_hist, clim_fut = storm_hist - anom_hist, storm_fut - anom_fut
        sign = -1 if var in SIGN_FLIP else 1
        dclim = sign * (clim_fut - clim_hist)
        danom = sign * (anom_fut - anom_hist)
        dstorm = sign * (storm_fut - storm_hist)
        denom = abs(dstorm)
        climate_pct = np.nan if np.isclose(denom, 0) else 100 * dclim / denom
        anomaly_pct = np.nan if np.isclose(denom, 0) else 100 * danom / denom
        rows.append((var, climate_pct, anomaly_pct, dstorm))
        print(f"{var:18s} dclim={dclim: .4g}  danom={danom: .4g}  dstorm={dstorm: .4g}")
    return pd.DataFrame(rows, columns=["variable", "climate", "anomaly", "total_change"])


# -----------------------------------------------------------------------------
# Plot helpers
# -----------------------------------------------------------------------------
def field_values(data):
    return data.squeeze().values if isinstance(data, xr.DataArray) else _as_float(data).squeeze()


def symmetric_levels(data, percentile=98, n=13):
    values = field_values(data)
    limit = np.nanpercentile(np.abs(values), percentile)
    if not np.isfinite(limit) or np.isclose(limit, 0):
        limit = 1.0
    return np.linspace(-limit, limit, n)


def spatial_xy(data):
    ny, nx = field_values(data).shape
    return ((np.arange(nx) - (nx - 1) / 2) * GRID_KM,
            (np.arange(ny) - (ny - 1) / 2) * GRID_KM)


def add_panel_label(ax, label, dark_background=False):
    """Place a black panel letter just outside the upper-left plot corner."""
    ax.text(-0.045, 1.035, f"({label})", transform=ax.transAxes,
            va="bottom", ha="right", color="black", fontsize=11,
            fontweight="bold", clip_on=False, zorder=20)


def add_contour_labels(ax, contour, color):
    labels = ax.clabel(contour, inline=True, inline_spacing=3, fontsize=6.8,
                       fmt="%g", colors=color)
    for label in labels:
        label.set_path_effects([pe.withStroke(linewidth=1.8,
                                              foreground="white" if color != "white" else "#333333")])


def add_sample_box(ax):
    half = SAMPLE_BOX_KM / 2
    ax.add_patch(Rectangle((-half, -half), SAMPLE_BOX_KM, SAMPLE_BOX_KM,
                           facecolor="none", edgecolor="white", linewidth=1.2,
                           linestyle=(0, (4, 2)), zorder=12,
                           path_effects=[pe.withStroke(linewidth=2.2, foreground="#333333")]))


def style_map_axis(ax, row, col):
    ticks = [-250, -125, 0, 125, 250]
    ax.set_aspect("equal")
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xlim(-251, 251)
    ax.set_ylim(-251, 251)
    ax.tick_params(length=3, width=0.7, color="#66717A")
    if row == 0:
        ax.tick_params(labelbottom=False)
    else:
        ax.set_xlabel("Distance from storm centre (km)")
    if col == 0:
        ax.set_ylabel("Distance from storm centre (km)")
    else:
        ax.tick_params(labelleft=False)
    ax.plot(0, 0, marker="o", ms=4.2, mfc="white", mec=COLORS["ink"], mew=1,
            zorder=15)


def add_colorbar(fig, ax, artist, label):
    # An inset axis keeps the colour bar attached to the map.  Letting the
    # colour bar participate in constrained_layout leaves a deceptively large
    # gap even when fig.colorbar's ``pad`` is very small.
    cax = ax.inset_axes([1.018, 0.04, 0.035, 0.92])
    cbar = fig.colorbar(artist, cax=cax, orientation="vertical")
    cbar.set_label(label, labelpad=6)
    cbar.ax.tick_params(length=2, pad=2, labelsize=7.5)
    cbar.outline.set_linewidth(0.6)


def add_overlay_note(ax, text, dark=False):
    ax.text(0.98, 0.025, text, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=6.8, color="white" if dark else COLORS["ink"], zorder=18,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "#20272D" if dark else "white",
                  "edgecolor": "none", "alpha": 0.76})


def plot_composites(fig, axes, hist, fut, hist_anom, fut_anom):
    t_hist = hist_anom["t_srfc"]
    t_change = fut_anom["t_srfc"] - hist_anom["t_srfc"]
    div_hist = calc_divergence(hist) * 1e5
    div_change = (calc_divergence(fut) - calc_divergence(hist)) * 1e5
    q_hist = hist_anom["q_srfc"] * 1000
    q_change = (fut_anom["q_srfc"] - hist_anom["q_srfc"]) * 1000
    x, y = spatial_xy(t_hist)

    # Filled fields: historic composite (top) and future - historic (bottom).
    fields = [t_hist, div_hist, q_hist, t_change, div_change, q_change]
    levels = [np.arange(-0.9, 0.91, 0.1), symmetric_levels(div_hist, 95),
              symmetric_levels(q_hist, 95), np.arange(-0.9, 0.91, 0.1),
              symmetric_levels(div_change, 95), symmetric_levels(q_change, 95)]
    cmaps = ["RdBu_r", "BrBG", "RdBu", "RdBu_r", "BrBG", "RdBu"]
    cbar_labels = [
        r"925-hPa T anomaly (K)", r"Divergence ($10^{-5}$ s$^{-1}$)",
        r"925-hPa q anomaly (g kg$^{-1}$)", r"Change in T anomaly (K)",
        r"Change in divergence ($10^{-5}$ s$^{-1}$)",
        r"Change in q anomaly (g kg$^{-1}$)",
    ]

    for i, (ax, data, lev, cmap, label) in enumerate(
            zip(axes.flat, fields, levels, cmaps, cbar_labels)):
        artist = ax.contourf(x, y, field_values(data), levels=lev, cmap=cmap,
                             extend="both")
        style_map_axis(ax, i // 3, i % 3)
        add_colorbar(fig, ax, artist, label)
        add_panel_label(ax, "abcdef"[i], dark_background=True)

    axes[0, 0].set_title("Surface temperature and land fluxes", pad=8)
    axes[0, 1].set_title("Low-level divergence and flow", pad=8)
    axes[0, 2].set_title("Near-surface humidity and flow", pad=8)

    # (a): SH and soil-moisture anomalies.
    cs = axes[0, 0].contour(x, y, field_values(hist_anom["sh"]),
                            levels=[5, 10, 15, 20], colors=COLORS["ink"], linewidths=1.0)
    add_contour_labels(axes[0, 0], cs, COLORS["ink"])
    cs = axes[0, 0].contour(x, y, field_values(hist_anom["SM"]),
                            levels=[-1.5, -1, 0, 1, 1.5], colors=COLORS["sm"], linewidths=1.0)
    add_contour_labels(axes[0, 0], cs, COLORS["sm"])
    add_sample_box(axes[0, 0])
    add_overlay_note(axes[0, 0], "Contours: SH (black), SM (blue)")

    # (b): T anomaly and anomalous surface winds.
    t_levels = [-0.9, -0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.7, 0.9]
    cs = axes[0, 1].contour(x, y, field_values(t_hist), levels=t_levels,
                            colors=COLORS["ink"], linewidths=0.9)
    add_contour_labels(axes[0, 1], cs, COLORS["ink"])
    skip = (slice(None, None, 10), slice(None, None, 10))
    q = axes[0, 1].quiver(x[::10], y[::10], field_values(hist_anom["u_srfc"])[skip],
                          field_values(hist_anom["v_srfc"])[skip], color="#343A40",
                          scale=14, width=0.003, headwidth=3.4, zorder=11)
    axes[0, 1].quiverkey(q, 0.87, 1.025, 1, r"1 m s$^{-1}$", labelpos="E",
                         coordinates="axes", fontproperties={"size": 7})
    add_overlay_note(axes[0, 1], "Contours: T anomaly; arrows: wind anomaly")

    # (c): T anomaly and historical mean meridional wind.
    cs = axes[0, 2].contour(x, y, field_values(t_hist), levels=t_levels,
                            colors=COLORS["ink"], linewidths=0.9)
    add_contour_labels(axes[0, 2], cs, COLORS["ink"])
    cs = axes[0, 2].contour(x, y, field_values(hist["v_srfc"]),
                            levels=np.arange(-1.5, 1.51, 0.5), colors="white", linewidths=1.0)
    add_contour_labels(axes[0, 2], cs, "white")
    add_overlay_note(axes[0, 2], "Contours: T anomaly (black), historic v (white)", dark=True)

    # (d): changes in SH and soil-moisture anomalies.
    cs = axes[1, 0].contour(x, y, field_values(fut_anom["sh"] - hist_anom["sh"]),
                            levels=[5, 10, 15], colors=COLORS["ink"], linewidths=1.0)
    add_contour_labels(axes[1, 0], cs, COLORS["ink"])
    cs = axes[1, 0].contour(x, y, field_values(fut_anom["SM"] - hist_anom["SM"]),
                            levels=[-1.5, -1, -0.5], colors=COLORS["sm"], linewidths=1.0)
    add_contour_labels(axes[1, 0], cs, COLORS["sm"])
    add_sample_box(axes[1, 0])
    add_overlay_note(axes[1, 0], "Contours: change in SH (black), SM (blue)")

    # (e): change in T anomaly and anomalous winds.
    cs = axes[1, 1].contour(x, y, field_values(t_change), levels=t_levels,
                            colors=COLORS["ink"], linewidths=0.9)
    add_contour_labels(axes[1, 1], cs, COLORS["ink"])
    du = field_values(fut_anom["u_srfc"] - hist_anom["u_srfc"])
    dv = field_values(fut_anom["v_srfc"] - hist_anom["v_srfc"])
    q = axes[1, 1].quiver(x[::10], y[::10], du[skip], dv[skip], color="#343A40",
                          scale=14, width=0.003, headwidth=3.4, zorder=11)
    axes[1, 1].quiverkey(q, 0.87, 1.025, 1, r"1 m s$^{-1}$", labelpos="E",
                         coordinates="axes", fontproperties={"size": 7})
    add_overlay_note(axes[1, 1], "Contours/arrows: future - historical anomaly")

    # (f): change in T anomaly and future mean meridional wind.  This preserves
    # the intended absolute-wind context of the original panel.
    cs = axes[1, 2].contour(x, y, field_values(t_change), levels=t_levels,
                            colors=COLORS["ink"], linewidths=0.9)
    add_contour_labels(axes[1, 2], cs, COLORS["ink"])
    cs = axes[1, 2].contour(x, y, field_values(fut["v_srfc"]),
                            levels=np.arange(-1.5, 1.51, 0.5), colors="white", linewidths=1.0)
    add_contour_labels(axes[1, 2], cs, "white")
    add_overlay_note(axes[1, 2], "Contours: change in T (black), future v (white)", dark=True)

    axes[0, 0].text(-0.35, 0.5, "HISTORICAL\nCOMPOSITE", transform=axes[0, 0].transAxes,
                    rotation=90, ha="center", va="center", color=COLORS["muted"],
                    fontsize=8, fontweight="semibold")
    axes[1, 0].text(-0.35, 0.5, "FUTURE - HISTORICAL", transform=axes[1, 0].transAxes,
                    rotation=90, ha="center", va="center", color=COLORS["muted"],
                    fontsize=8, fontweight="semibold")


def plot_bar_panel(ax, contributions):
    x = np.arange(len(contributions))
    climate = contributions["climate"].to_numpy()
    anomaly = contributions["anomaly"].to_numpy()
    climate_pos, climate_neg = np.maximum(climate, 0), np.minimum(climate, 0)
    anomaly_pos, anomaly_neg = np.maximum(anomaly, 0), np.minimum(anomaly, 0)

    width = 0.64
    bar_kw = {"width": width, "edgecolor": "white", "linewidth": 0.75, "zorder": 3}
    ax.bar(x, climate_pos, color=COLORS["blue"], **bar_kw)
    ax.bar(x, climate_neg, color=COLORS["blue"], **bar_kw)
    ax.bar(x, anomaly_pos, bottom=climate_pos, color=COLORS["orange"], hatch="////", **bar_kw)
    ax.bar(x, anomaly_neg, bottom=climate_neg, color=COLORS["orange"], hatch="////", **bar_kw)

    # The diamond is the net of components; it sits at +100 or -100 unless the
    # total storm-day change is numerically zero.
    net = climate + anomaly
    ax.scatter(x, net, marker="D", s=24, color=COLORS["ink"], edgecolor="white",
               linewidth=0.6, zorder=5)

    ax.axhline(0, color=COLORS["ink"], lw=0.9, zorder=2)
    ax.axhline(100, color="#AAB1B7", lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.axhline(-100, color="#AAB1B7", lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.yaxis.grid(True, color="#D9DEE2", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)

    extent = np.nanmax(np.abs(np.r_[climate_pos + anomaly_pos,
                                    climate_neg + anomaly_neg, net]))
    ylim = max(125, 25 * np.ceil((extent + 12) / 25))
    ax.set_ylim(-ylim, ylim)
    ax.set_xlim(-0.65, len(x) - 0.35)
    ax.set_xticks(x, [BAR_LABELS[v] for v in contributions["variable"]])
    ax.set_ylabel("Contribution relative to |total storm-day change| (%)")
    ax.set_title("Storm-day change: background climate shift and anomaly enhancement",
                 loc="left", pad=10)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", length=0, pad=6)
    add_panel_label(ax, "g")

    # Soft separators make the variable families easier to scan without adding
    # another set of axes or a heavy table-like frame.
    for xpos in (4.5, 8.5):
        ax.axvline(xpos, color="#C9CED2", lw=0.7, ls=(0, (2, 3)), zorder=1)
    transform = ax.get_xaxis_transform()
    ax.text(2.0, -0.12, "LAND / SURFACE FLUXES", transform=transform, ha="center",
            va="top", fontsize=7, color=COLORS["muted"], fontweight="semibold")
    ax.text(6.5, -0.12, "THERMODYNAMIC ENVIRONMENT", transform=transform, ha="center",
            va="top", fontsize=7, color=COLORS["muted"], fontweight="semibold")
    ax.text(9.0, -0.12, "DYNAMICS", transform=transform, ha="center", va="top",
            fontsize=7, color=COLORS["muted"], fontweight="semibold")

    handles = [
        Patch(facecolor=COLORS["blue"], edgecolor="none", label="Background climate shift"),
        Patch(facecolor=COLORS["orange"], edgecolor="white", hatch="////",
              label="Storm-day anomaly enhancement"),
        Line2D([], [], marker="D", linestyle="none", ms=5, color=COLORS["ink"],
               label="Net storm-day change"),
    ]
    ax.legend(handles=handles, frameon=False, ncol=3, loc="upper right",
              bbox_to_anchor=(1, 1.12), handlelength=1.5, columnspacing=1.6)


def make_figure(hist, fut, hist_anom, fut_anom, contributions):
    fig = plt.figure(figsize=(15.2, 11.5), dpi=160, layout="constrained")
    grid = fig.add_gridspec(3, 3, height_ratios=(1, 1, 0.80), hspace=0.02, wspace=0.04)
    map_axes = np.array([[fig.add_subplot(grid[r, c]) for c in range(3)] for r in range(2)])
    bar_ax = fig.add_subplot(grid[2, :])

    plot_composites(fig, map_axes, hist, fut, hist_anom, fut_anom)
    plot_bar_panel(bar_ax, contributions)
    return fig


# -----------------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------------
def main():
    da_hist_anom = mean_hours("hist", "anom")
    da_fut_anom = mean_hours("fut", "anom")
    da_hist = mean_hours("hist", "mean")
    da_fut = mean_hours("fut", "mean")  # both 17 and 18 UTC are future files

    tables = load_tables()
    lhist, lahist = tables["hist"]
    lfut, lafut = tables["fut"]
    contributions = calculate_bar_contributions(lhist, lahist, lfut, lafut)

    fig = make_figure(da_hist, da_fut, da_hist_anom, da_fut_anom, contributions)
    fig.savefig(OUTPUT_FIG, dpi=300, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
