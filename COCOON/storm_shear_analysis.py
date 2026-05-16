"""
Convective storm shear–moisture interaction analysis
=====================================================
Fits a mixed-effects interaction model:
    Pmax / Area = β₁·TCW + β₂·shear + β₃·(TCW × shear)
with region as a random effect, then produces:
  1. Forest plot of regional β₃ (area and precip)
  2. Within-region scatter plots with regression lines + CIs
  3. Global storm frequency map with region bounding boxes
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy import stats

# statsmodels mixed-effects
import statsmodels.formula.api as smf
import statsmodels.api as sm

# ── Region bounding boxes ────────────────────────────────────────────────────
# Each entry: [lon_min, lon_max, lat_min, lat_max], label
MREGIONS = {
    "GPlains":  [[-105, -85,  25,  50], "Great Plains"],
    "china":    [[ 100, 125,  20,  45], "China"],
    "india":    [[  65,  90,   5,  30], "India"],
    "WAf":      [[ -15,  25,   0,  20], "West Africa"],
    "australia":[[ 115, 155, -35, -10], "Australia"],
    "SAf":      [[  15,  40, -35,  -5], "South Africa"],
    "sub_SA":   [[ -75, -40, -35,   5], "Sub-trop S.Am."],
    "trop_SA":  [[ -75, -40,   5,  15], "Trop. S.Am."],
    "CAf":      [[  10,  35,  -5,  10], "Central Africa"],
    "EAf":      [[  30,  50, -10,  15], "East Africa"],
    "Atl":      [[ -60, -20,   5,  25], "Atlantic"],
    "Pcf":      [[ 140, 180, -20,   5], "Pacific"],
    "InO":      [[  60,  90, -15,  10], "Indian Ocean"],
}
REGION_COLORS = {k: plt.cm.tab20(i / len(MREGIONS))
                 for i, k in enumerate(MREGIONS)}

# ── Load data ────────────────────────────────────────────────────────────────

def load_data(path):
    """Load storm tracking CSV. Adjust path / sep as needed."""
    df = pd.read_csv(path)
    # Normalise column names to lowercase
    df.columns = df.columns.str.lower().str.strip()
    return df


def prepare(df: pd.DataFrame,
            moisture_col: str = "tcwv",
            area_col: str = "area",
            precip_col: str = "prcp",
            shear_col: str = "shear") -> pd.DataFrame:
    """
    Filter, engineer features, centre moisture at global mean.
    Returns a clean DataFrame ready for modelling.
    """
    required = [moisture_col, area_col, precip_col, shear_col, "region"]
    df = df.dropna(subset=required).copy()

    # Keep physically plausible values
    df = df[(df[area_col] > 0) & (df[precip_col] > 0) & (df[shear_col] >= 0)]

    # Centre moisture (so β₂ is interpretable at mean TCW)
    global_mean_tcw = df[moisture_col].mean()
    df["tcw_c"] = df[moisture_col] - global_mean_tcw  # centred

    # Standardise shear (z-score) for comparable β magnitudes
    df["shear_z"] = (df[shear_col] - df[shear_col].mean()) / df[shear_col].std()

    # Interaction term
    df["interact"] = df["tcw_c"] * df["shear_z"]

    # Log-transform outcomes (right-skewed)
    df["log_area"]  = np.log(df[area_col])
    df["log_prcp"]  = np.log(df[precip_col])

    return df, global_mean_tcw


# ── Mixed-effects model ──────────────────────────────────────────────────────
def fit_me_model(df: pd.DataFrame, outcome: str) -> smf.mixedlm:
    """
    Fit:
        outcome ~ tcw_c + shear_z + tcw_c:shear_z
        random intercept + random slope on interact per region.
    Returns fitted MixedLMResults object.
    """
    formula = f"{outcome} ~ tcw_c + shear_z + tcw_c:shear_z"
    md = smf.mixedlm(
        formula, df,
        groups=df["region"],
        exog_re=df[["interact"]]   # random slope on interaction
    )
    try:
        result = md.fit(reml=True, method="lbfgs", maxiter=500)
    except Exception:
        # Fall back to simpler random-intercept-only if convergence fails
        md = smf.mixedlm(formula, df, groups=df["region"])
        result = md.fit(reml=True, method="lbfgs", maxiter=500)
    return result


def regional_beta3(df: pd.DataFrame,
                   outcome: str) -> pd.DataFrame:
    """
    Fit OLS per region to obtain region-specific β₃ estimates + 95 % CI.
    These are the values shown in the forest plot.
    """
    records = []
    for region, grp in df.groupby("region"):
        if len(grp) < 20:
            continue
        X = sm.add_constant(grp[["tcw_c", "shear_z", "interact"]])
        y = grp[outcome]
        try:
            ols = sm.OLS(y, X).fit()
            coef  = ols.params["interact"]
            ci_lo = ols.conf_int().loc["interact", 0]
            ci_hi = ols.conf_int().loc["interact", 1]
            p     = ols.pvalues["interact"]
            records.append(dict(region=region, beta3=coef,
                                ci_lo=ci_lo, ci_hi=ci_hi, pval=p,
                                n=len(grp)))
        except Exception:
            pass
    return pd.DataFrame(records).sort_values("beta3", ascending=True)


# ── 1. Forest plot ───────────────────────────────────────────────────────────
def forest_plot(df: pd.DataFrame, figsize=(10, 7)):
    """Forest plot of regional β₃ for both area and precipitation."""
    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=False)

    outcomes = [("log_area", "β₃  (storm area)", "#2166ac"),
                ("log_prcp", "β₃  (max precip)", "#d6604d")]

    for ax, (outcome, xlabel, color) in zip(axes, outcomes):
        res = regional_beta3(df, outcome)
        y   = np.arange(len(res))
        ax.axvline(0, color="k", lw=0.8, ls="--", alpha=0.6)

        for i, row in res.reset_index(drop=True).iterrows():
            sig  = row["pval"] < 0.05
            ec   = REGION_COLORS.get(row["region"], "grey")
            ax.plot([row["ci_lo"], row["ci_hi"]], [i, i],
                    color=ec, lw=2, alpha=0.8)
            ax.scatter(row["beta3"], i,
                       color=ec,
                       marker="D" if sig else "o",
                       s=60 if sig else 35,
                       zorder=5, edgecolors="k", linewidths=0.4)

        ax.set_yticks(y)
        ax.set_yticklabels(
            [MREGIONS.get(r, {1: r})[1]
             if isinstance(MREGIONS.get(r, r), list)
             else r
             for r in res["region"]])
        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_title(xlabel, fontsize=12, pad=8)
        ax.tick_params(axis="y", labelsize=9)
        ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle("Regional β₃ (TCW × shear interaction)\n"
                 "Diamonds = p < 0.05; bars = 95 % CI",
                 fontsize=13, y=1.01)
    plt.tight_layout()
    return fig


# ── 2. Within-region scatter plots ──────────────────────────────────────────
def _plot_one_region(ax, grp, outcome, region_name, color):
    """Scatter + OLS regression line + 95 % CI band on centred TCW."""
    x = grp["tcw_c"].values
    y = grp[outcome].values

    ax.scatter(x, y, alpha=0.25, s=8, color=color, rasterized=True)

    # OLS fit for overlay
    X  = sm.add_constant(x)
    ols = sm.OLS(y, X).fit()
    x_grid = np.linspace(x.min(), x.max(), 200)
    X_grid = sm.add_constant(x_grid)
    pred   = ols.get_prediction(X_grid).summary_frame(alpha=0.05)

    ax.plot(x_grid, pred["mean"], color=color, lw=1.8)
    ax.fill_between(x_grid,
                    pred["mean_ci_lower"],
                    pred["mean_ci_upper"],
                    color=color, alpha=0.2)

    ax.set_title(region_name, fontsize=9, pad=3)
    ax.tick_params(labelsize=7)
    ax.spines[["top", "right"]].set_visible(False)


def scatter_panels(df: pd.DataFrame,
                   outcome: str = "log_area",
                   ncols: int = 4,
                   figsize=(14, 10)):
    """Grid of within-region scatter plots (TCW_c vs log outcome)."""
    regions   = sorted(df["region"].unique())
    nrows     = int(np.ceil(len(regions) / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=figsize,
                             constrained_layout=True)
    axes_flat = axes.flatten()

    ylabel = "log(storm area)" if outcome == "log_area" else "log(max precip)"

    for i, region in enumerate(regions):
        grp  = df[df["region"] == region]
        name = (MREGIONS[region][1]
                if region in MREGIONS and isinstance(MREGIONS[region], list)
                else region)
        color = REGION_COLORS.get(region, "steelblue")
        _plot_one_region(axes_flat[i], grp, outcome, name, color)
        if i % ncols == 0:
            axes_flat[i].set_ylabel(ylabel, fontsize=7)

    # Hide unused subplots
    for j in range(i + 1, len(axes_flat)):
        axes_flat[j].set_visible(False)

    fig.supxlabel("Centred TCW  (kg m⁻²)", fontsize=10)
    fig.suptitle(f"Within-region: centred TCW vs {ylabel}\n"
                 "(line = OLS; shading = 95 % CI)",
                 fontsize=12)
    return fig


# ── 3. Global storm frequency map ───────────────────────────────────────────
def frequency_map(df: pd.DataFrame,
                  lon_col: str = "lon",
                  lat_col: str = "lat",
                  bin_deg: float = 2.0,
                  figsize=(14, 7)):
    """
    2-D histogram of storm frequency + region bounding boxes.
    Requires lon/lat columns in df; falls back to centroid columns if named
    differently (lon_centroid, lat_centroid, clon, clat).
    """
    # Flexible column detection
    for lc in [lon_col, "lon_centroid", "clon", "longitude"]:
        if lc in df.columns:
            lon_col = lc; break
    for lac in [lat_col, "lat_centroid", "clat", "latitude"]:
        if lac in df.columns:
            lat_col = lac; break

    fig = plt.figure(figsize=figsize)
    ax  = fig.add_subplot(1, 1, 1,
                          projection=ccrs.Robinson(central_longitude=0))
    ax.set_global()
    ax.add_feature(cfeature.LAND,   facecolor="#f0efe4", zorder=0)
    ax.add_feature(cfeature.OCEAN,  facecolor="#d0e8f0", zorder=0)
    ax.add_feature(cfeature.COASTLINE, lw=0.4, zorder=1)
    ax.add_feature(cfeature.BORDERS,   lw=0.2, zorder=1)
    ax.gridlines(draw_labels=False, lw=0.3, color="grey", alpha=0.5)

    if lon_col in df.columns and lat_col in df.columns:
        lons = df[lon_col].dropna().values
        lats = df[lat_col].dropna().values
        lon_bins = np.arange(-180, 181, bin_deg)
        lat_bins = np.arange( -90,  91, bin_deg)
        H, _, _ = np.histogram2d(lons, lats,
                                  bins=[lon_bins, lat_bins])
        lon_centres = 0.5 * (lon_bins[:-1] + lon_bins[1:])
        lat_centres = 0.5 * (lat_bins[:-1] + lat_bins[1:])
        LON, LAT = np.meshgrid(lon_centres, lat_centres)
        pcm = ax.pcolormesh(LON, LAT, H.T,
                            cmap="YlOrRd",
                            norm=plt.Normalize(vmin=0,
                                               vmax=np.percentile(H[H > 0], 95)),
                            transform=ccrs.PlateCarree(),
                            zorder=2, alpha=0.85)
        cbar = plt.colorbar(pcm, ax=ax, shrink=0.55, pad=0.02)
        cbar.set_label("Storm count per 2° cell", fontsize=9)
    else:
        ax.text(0, 0, "No lon/lat columns found\n"
                "(add lon/lat centroid columns to df)",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=12, color="red")

    # Draw region bounding boxes
    for key, (bbox, label) in {k: (v[0], v[1])
                                for k, v in MREGIONS.items()}.items():
        lon0, lon1, lat0, lat1 = bbox
        color = REGION_COLORS[key]
        rect_lons = [lon0, lon1, lon1, lon0, lon0]
        rect_lats = [lat0, lat0, lat1, lat1, lat0]
        ax.plot(rect_lons, rect_lats,
                transform=ccrs.PlateCarree(),
                color=color, lw=1.8, zorder=3)
        ax.text((lon0 + lon1) / 2, lat1 + 1, label,
                transform=ccrs.PlateCarree(),
                fontsize=6, ha="center", va="bottom",
                color=color, fontweight="bold", zorder=4,
                bbox=dict(boxstyle="round,pad=0.1",
                          fc="white", ec="none", alpha=0.6))

    ax.set_title("Global convective storm frequency  +  analysis regions",
                 fontsize=13, pad=10)
    return fig


# ── Main ─────────────────────────────────────────────────────────────────────
def main(data_path: str,
         moisture_col: str = "tcwv",   # swap to "q925" to test Alfaro signal
         save_figs: bool = True):

    print("Loading data …")
    raw = load_data(data_path)
    df, global_mean_tcw = prepare(raw, moisture_col=moisture_col)
    print(f"  {len(df):,} storms across {df['region'].nunique()} regions")
    print(f"  Global mean {moisture_col}: {global_mean_tcw:.2f}")

    # ── Mixed-effects models ────────────────────────────────────────────────
    print("\nFitting mixed-effects model: log_area …")
    me_area = fit_me_model(df, "log_area")
    print(me_area.summary())

    print("\nFitting mixed-effects model: log_prcp …")
    me_prcp = fit_me_model(df, "log_prcp")
    print(me_prcp.summary())

    # Extract global β₃ + SE from mixed-effects models
    for label, result in [("AREA", me_area), ("PRECIP", me_prcp)]:
        b3   = result.params.get("tcw_c:shear_z",
               result.params.get("interact", np.nan))
        se3  = result.bse.get("tcw_c:shear_z",
               result.bse.get("interact", np.nan))
        pv3  = result.pvalues.get("tcw_c:shear_z",
               result.pvalues.get("interact", np.nan))
        print(f"\n  [{label}] global β₃ = {b3:.4f}  "
              f"SE = {se3:.4f}  p = {pv3:.4g}")

    # ── Figures ─────────────────────────────────────────────────────────────
    print("\nGenerating forest plot …")
    fig1 = forest_plot(df)
    if save_figs:
        fig1.savefig("forest_beta3.png", dpi=180, bbox_inches="tight")

    print("Generating area scatter panels …")
    fig2 = scatter_panels(df, outcome="log_area")
    if save_figs:
        fig2.savefig("scatter_area.png", dpi=150, bbox_inches="tight")

    print("Generating precip scatter panels …")
    fig3 = scatter_panels(df, outcome="log_prcp")
    if save_figs:
        fig3.savefig("scatter_prcp.png", dpi=150, bbox_inches="tight")

    print("Generating frequency map …")
    fig4 = frequency_map(df)
    if save_figs:
        fig4.savefig("storm_frequency_map.png", dpi=180, bbox_inches="tight")

    plt.show()
    return me_area, me_prcp, df


if __name__ == "__main__":
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else "storms.csv"
    main(path)
