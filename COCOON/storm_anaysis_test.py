"""
Synthetic data demo — convective storm shear–moisture interaction analysis
Run directly: python demo_plots.py
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import warnings
warnings.filterwarnings("ignore")

# ── Region definitions ───────────────────────────────────────────────────────
MREGIONS = {
    "GPlains":   [[-105, -85,  25,  50], "Great Plains"],
    "china":     [[ 100, 125,  20,  45], "China"],
    "india":     [[  65,  90,   5,  30], "India"],
    "WAf":       [[ -15,  25,   0,  20], "West Africa"],
    "australia": [[ 115, 155, -35, -10], "Australia"],
    "SAf":       [[  15,  40, -35,  -5], "South Africa"],
    "sub_SA":    [[ -75, -40, -35,   5], "Sub-trop S.Am."],
    "trop_SA":   [[ -75, -40,   5,  15], "Trop. S.Am."],
    "CAf":       [[  10,  35,  -5,  10], "Central Africa"],
    "EAf":       [[  30,  50, -10,  15], "East Africa"],
    "Atl":       [[ -60, -20,   5,  25], "Atlantic"],
    "Pcf":       [[ 140, 180, -20,   5], "Pacific"],
    "InO":       [[  60,  90, -15,  10], "Indian Ocean"],
}

REGION_COLORS = {k: plt.cm.tab20(i / len(MREGIONS))
                 for i, k in enumerate(MREGIONS)}

# ── Synthetic data generation ────────────────────────────────────────────────
def make_synthetic(n_per_region=300, seed=42):
    """
    Generate physically plausible synthetic storm data.

    Each region gets its own TCW climatology, shear climatology, and
    a true β₃ drawn from a distribution that is mostly positive
    (mimicking the Alfaro mechanism) but allows two negative outliers
    (China, SAf) as described in the scientific context.
    """
    rng = np.random.default_rng(seed)

    # Regional climatological parameters (mean_tcw, mean_shear, true_beta3)
    region_params = {
        "GPlains":   dict(mu_tcw=28, mu_shear=12, beta3= 0.018),
        "china":     dict(mu_tcw=32, mu_shear=10, beta3=-0.012),  # dry mid-levels
        "india":     dict(mu_tcw=42, mu_shear= 8, beta3= 0.022),
        "WAf":       dict(mu_tcw=45, mu_shear= 6, beta3= 0.030),
        "australia": dict(mu_tcw=25, mu_shear=11, beta3= 0.010),
        "SAf":       dict(mu_tcw=22, mu_shear=14, beta3=-0.008),  # dry mid-levels
        "sub_SA":    dict(mu_tcw=35, mu_shear= 9, beta3= 0.020),
        "trop_SA":   dict(mu_tcw=48, mu_shear= 5, beta3= 0.028),
        "CAf":       dict(mu_tcw=50, mu_shear= 4, beta3= 0.035),
        "EAf":       dict(mu_tcw=38, mu_shear= 7, beta3= 0.016),
        "Atl":       dict(mu_tcw=44, mu_shear= 6, beta3= 0.025),
        "Pcf":       dict(mu_tcw=46, mu_shear= 5, beta3= 0.027),
        "InO":       dict(mu_tcw=43, mu_shear= 6, beta3= 0.023),
    }

    frames = []
    for region, p in region_params.items():
        bbox = MREGIONS[region][0]
        lon0, lon1, lat0, lat1 = bbox

        n = n_per_region
        tcwv  = rng.normal(p["mu_tcw"],  6.0, n).clip(5, 80)
        shear = rng.normal(p["mu_shear"], 4.0, n).clip(0, 35)

        # Log-area: linear in tcwv + shear + interaction + noise
        log_area = (
            8.5
            + 0.04  * (tcwv  - p["mu_tcw"])
            + 0.03  * (shear - p["mu_shear"])
            + p["beta3"] * (tcwv - p["mu_tcw"]) * (shear - p["mu_shear"])
            + rng.normal(0, 0.6, n)
        )

        # Log-precip: similar structure, slightly different coefficients
        log_prcp = (
            3.2
            + 0.025 * (tcwv  - p["mu_tcw"])
            + 0.015 * (shear - p["mu_shear"])
            + p["beta3"] * 0.8 * (tcwv - p["mu_tcw"]) * (shear - p["mu_shear"])
            + rng.normal(0, 0.4, n)
        )

        lons = rng.uniform(lon0, lon1, n)
        lats = rng.uniform(lat0, lat1, n)

        frames.append(pd.DataFrame({
            "region":   region,
            "tcwv":     tcwv,
            "shear":    shear,
            "log_area": log_area,
            "log_prcp": log_prcp,
            "lon":      lons,
            "lat":      lats,
        }))

    df = pd.concat(frames, ignore_index=True)

    # Centre TCW and z-score shear globally
    df["tcw_c"]   = df["tcwv"] - df["tcwv"].mean()
    df["shear_z"] = (df["shear"] - df["shear"].mean()) / df["shear"].std()
    df["interact"] = df["tcw_c"] * df["shear_z"]

    return df


# ── Regional OLS β₃ table ────────────────────────────────────────────────────
def regional_beta3(df, outcome):
    records = []
    for region, grp in df.groupby("region"):
        X = sm.add_constant(grp[["tcw_c", "shear_z", "interact"]])
        y = grp[outcome]
        ols   = sm.OLS(y, X).fit()
        ci    = ols.conf_int()
        records.append(dict(
            region = region,
            beta3  = ols.params["interact"],
            ci_lo  = ci.loc["interact", 0],
            ci_hi  = ci.loc["interact", 1],
            pval   = ols.pvalues["interact"],
            n      = len(grp),
        ))
    return pd.DataFrame(records).sort_values("beta3", ascending=True).reset_index(drop=True)


# ── Plot 1: Forest plot ──────────────────────────────────────────────────────
def plot_forest(df):
    fig, axes = plt.subplots(1, 2, figsize=(11, 6))
    fig.patch.set_facecolor("#fafafa")

    configs = [
        ("log_area", "β₃  for log(storm area)",  "#2166ac"),
        ("log_prcp", "β₃  for log(max precip)",  "#d6604d"),
    ]

    for ax, (outcome, xlabel, spine_color) in zip(axes, configs):
        res = regional_beta3(df, outcome)
        y   = np.arange(len(res))

        ax.set_facecolor("#fafafa")
        ax.axvline(0, color="#555", lw=1.0, ls="--", alpha=0.7, zorder=1)

        for i, row in res.iterrows():
            sig   = row["pval"] < 0.05
            color = REGION_COLORS.get(row["region"], "grey")
            # CI bar
            ax.plot([row["ci_lo"], row["ci_hi"]], [i, i],
                    color=color, lw=2.2, alpha=0.9, solid_capstyle="round", zorder=2)
            # Point estimate
            ax.scatter(row["beta3"], i,
                       color=color,
                       marker="D" if sig else "o",
                       s=70 if sig else 40,
                       edgecolors="white", linewidths=0.6,
                       zorder=3)
            # n label
            ax.text(row["ci_hi"] + 0.0005, i, f"n={row['n']}",
                    va="center", fontsize=7, color="#666")

        label_map = {k: v[1] for k, v in MREGIONS.items()}
        ax.set_yticks(y)
        ax.set_yticklabels([label_map.get(r, r) for r in res["region"]], fontsize=9)
        ax.set_xlabel(xlabel, fontsize=10)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines["left"].set_color(spine_color)
        ax.spines["bottom"].set_color("#aaa")
        ax.tick_params(axis="x", labelsize=8)
        ax.grid(axis="x", lw=0.4, color="#ccc", alpha=0.7)

    # Legend
    sig_patch = plt.Line2D([0], [0], marker="D", color="w",
                           markerfacecolor="#555", markersize=8, label="p < 0.05")
    ns_patch  = plt.Line2D([0], [0], marker="o", color="w",
                           markerfacecolor="#555", markersize=7, label="p ≥ 0.05")
    fig.legend(handles=[sig_patch, ns_patch],
               loc="lower center", ncol=2, fontsize=9,
               frameon=False, bbox_to_anchor=(0.5, -0.02))

    fig.suptitle("Regional β₃  (TCW × shear interaction)\n"
                 "bars = 95 % CI,  diamonds = significant",
                 fontsize=12, y=1.01)
    plt.tight_layout()
    fig.savefig("forest_beta3.png", dpi=180, bbox_inches="tight")
    print("  Saved: forest_beta3.png")
    return fig


# ── Plot 2 & 3: Scatter panels ───────────────────────────────────────────────
def plot_scatter_panels(df, outcome="log_area", ncols=4):
    label_map = {k: v[1] for k, v in MREGIONS.items()}
    regions   = sorted(df["region"].unique())
    nrows     = int(np.ceil(len(regions) / ncols))
    ylabel    = "log(storm area)" if outcome == "log_area" else "log(max precip)"

    fig, axes = plt.subplots(nrows, ncols, figsize=(13, 8),
                             constrained_layout=True)
    fig.patch.set_facecolor("#fafafa")
    axes_flat = axes.flatten()

    for i, region in enumerate(regions):
        ax    = axes_flat[i]
        grp   = df[df["region"] == region]
        color = REGION_COLORS.get(region, "steelblue")

        x = grp["tcw_c"].values
        y = grp[outcome].values

        ax.scatter(x, y, alpha=0.20, s=6, color=color, rasterized=True)

        # OLS + CI band
        X_fit  = sm.add_constant(x)
        ols    = sm.OLS(y, X_fit).fit()
        x_grid = np.linspace(x.min(), x.max(), 200)
        X_grid = sm.add_constant(x_grid)
        pred   = ols.get_prediction(X_grid).summary_frame(alpha=0.05)

        ax.plot(x_grid, pred["mean"], color=color, lw=1.8)
        ax.fill_between(x_grid,
                        pred["mean_ci_lower"],
                        pred["mean_ci_upper"],
                        color=color, alpha=0.25)
        ax.axvline(0, color="#aaa", lw=0.7, ls=":")

        # β₃ annotation
        X2  = sm.add_constant(grp[["tcw_c", "shear_z", "interact"]])
        o2  = sm.OLS(y, X2).fit()
        b3  = o2.params["interact"]
        pv  = o2.pvalues["interact"]
        star = "**" if pv < 0.01 else ("*" if pv < 0.05 else "")
        ax.text(0.97, 0.05, f"β₃={b3:.3f}{star}",
                transform=ax.transAxes, ha="right", fontsize=7,
                color=color, fontweight="bold")

        ax.set_facecolor("#fafafa")
        ax.set_title(label_map.get(region, region), fontsize=9, pad=3)
        ax.tick_params(labelsize=7)
        ax.spines[["top", "right"]].set_visible(False)
        if i % ncols == 0:
            ax.set_ylabel(ylabel, fontsize=7)

    for j in range(i + 1, len(axes_flat)):
        axes_flat[j].set_visible(False)

    fig.supxlabel("Centred TCW  (kg m⁻²)", fontsize=10)
    fig.suptitle(f"Within-region: centred TCW vs {ylabel}\n"
                 f"(OLS line + 95 % CI;  β₃ from full interaction model;  * p<0.05  ** p<0.01)",
                 fontsize=11)

    fname = f"scatter_{outcome.replace('log_','')}.png"
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    print(f"  Saved: {fname}")
    return fig


# ── Plot 4: Global frequency map ─────────────────────────────────────────────
def plot_frequency_map(df, bin_deg=2.0):
    fig = plt.figure(figsize=(14, 6.5))
    ax  = fig.add_subplot(1, 1, 1,
                          projection=ccrs.Robinson(central_longitude=10))
    ax.set_global()

    ax.add_feature(cfeature.LAND,      facecolor="#f0ede4", zorder=0)
    ax.add_feature(cfeature.OCEAN,     facecolor="#d4eaf5", zorder=0)
    ax.add_feature(cfeature.COASTLINE, lw=0.5, edgecolor="#888", zorder=2)
    ax.add_feature(cfeature.BORDERS,   lw=0.3, edgecolor="#aaa", zorder=2)
    ax.gridlines(draw_labels=False, lw=0.3, color="grey", alpha=0.4)

    # 2-D histogram
    lons = df["lon"].values
    lats = df["lat"].values
    lon_bins = np.arange(-180, 181, bin_deg)
    lat_bins = np.arange( -90,  91, bin_deg)
    H, _, _ = np.histogram2d(lons, lats, bins=[lon_bins, lat_bins])

    lon_c = 0.5 * (lon_bins[:-1] + lon_bins[1:])
    lat_c = 0.5 * (lat_bins[:-1] + lat_bins[1:])
    LON, LAT = np.meshgrid(lon_c, lat_c)

    pcm = ax.pcolormesh(
        LON, LAT, H.T,
        cmap="YlOrRd",
        vmin=0, vmax=np.percentile(H[H > 0], 95),
        transform=ccrs.PlateCarree(),
        zorder=1, alpha=0.88,
    )
    cbar = plt.colorbar(pcm, ax=ax, shrink=0.5, pad=0.02, aspect=20)
    cbar.set_label(f"Storm count per {bin_deg}° cell", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    # Region bounding boxes
    for key, (bbox, label) in {k: (v[0], v[1]) for k, v in MREGIONS.items()}.items():
        lon0, lon1, lat0, lat1 = bbox
        color = REGION_COLORS[key]
        rect_lons = [lon0, lon1, lon1, lon0, lon0]
        rect_lats = [lat0, lat0, lat1, lat1, lat0]
        ax.plot(rect_lons, rect_lats,
                transform=ccrs.PlateCarree(),
                color=color, lw=2.0, zorder=4)
        ax.text((lon0 + lon1) / 2, lat1 + 1.2, label,
                transform=ccrs.PlateCarree(),
                fontsize=6.5, ha="center", va="bottom",
                color=color, fontweight="bold", zorder=5,
                bbox=dict(boxstyle="round,pad=0.15",
                          fc="white", ec="none", alpha=0.75))

    ax.set_title("Global convective storm frequency  |  analysis regions",
                 fontsize=13, pad=10)

    fig.savefig("storm_frequency_map.png", dpi=180, bbox_inches="tight")
    print("  Saved: storm_frequency_map.png")
    return fig


# ── Main ─────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("Generating synthetic storm dataset …")
    df = make_synthetic(n_per_region=300)
    print(f"  {len(df):,} storms  |  {df['region'].nunique()} regions\n")

    print("Plot 1: Forest plot …")
    plot_forest(df)

    print("Plot 2: Area scatter panels …")
    plot_scatter_panels(df, outcome="log_area")

    print("Plot 3: Precip scatter panels …")
    plot_scatter_panels(df, outcome="log_prcp")

    print("Plot 4: Frequency map …")
    plot_frequency_map(df)

    print("\nDone. Showing figures …")
    plt.show()
