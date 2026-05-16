"""
Shear scaling plots for low / mid / high TCW decile bins
across convective regions.
 
Expected input structure:
    regions = {
        'west africa': pd.DataFrame,   # columns defined in CONFIG below
        'india':       pd.DataFrame,
        ...
    }
 
Call:
    from shear_scaling_plots import make_plots
    make_plots(regions, outdir='/your/path')
 
Outputs:
    shear_scaling_<AREA_COL>.png
    shear_scaling_<PRECIP_COL>.png
"""
 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
 
# ── CONFIG — only edit here ───────────────────────────────────────────────────
SHEAR_COL  = 'shear'
TCW_COL    = 'tcwv'
PRECIP_COL = 'prcp'
AREA_COL   = 'area'
 
TCW_BINS = [
    ('Low\n(0–10th)',    0,  10),
    ('Mid\n(45–55th)',  45,  55),
    ('High\n(90–100th)',90, 100),
]
 
COLORS = {
    'Low\n(0–10th)':    '#4393c3',
    'Mid\n(45–55th)':   '#f4a582',
    'High\n(90–100th)': '#d6604d',
}
# ─────────────────────────────────────────────────────────────────────────────
 
 
def shear_slope(df, response_col, tcw_lo, tcw_hi):
    """
    OLS slope of response_col ~ SHEAR_COL within a TCW percentile band.
    Returns (slope, 95%-CI half-width), or (nan, nan) if too few points.
    """

    lo = np.percentile(df[TCW_COL], tcw_lo)
    hi = np.percentile(df[TCW_COL], tcw_hi)
    sub = df[(df[TCW_COL] >= lo) & (df[TCW_COL] <= hi)].dropna(
        subset=[SHEAR_COL, response_col])
 
    if len(sub) < 5:
        return np.nan, np.nan
 
    slope, _, _, _, se = stats.linregress(sub[SHEAR_COL], sub[response_col])
    ci = stats.t.ppf(0.975, df=len(sub) - 2) * se
    return slope, ci
 
 
def plot_shear_scaling(data, response_col, ylabel, title, outpath):
    regions   = list(data.keys())
    n_regions = len(regions)
 
    bar_width = 0.22
    x = np.arange(n_regions) * (len(TCW_BINS) * bar_width + 0.25)
 
    fig, ax = plt.subplots(figsize=(max(10, n_regions * 1.8), 5.5))
 
    for i, (label, lo, hi) in enumerate(TCW_BINS):
        slopes, cis = [], []
        for region in regions:
           
            s, c = shear_slope(data[region], response_col, lo, hi)
            slopes.append(s)
            cis.append(c)
 
        offsets = x + (i - 1) * bar_width
        ax.bar(offsets, slopes, width=bar_width,
               label=label, color=COLORS[label],
               edgecolor='white', linewidth=0.6, zorder=3)
 
        for j, (s, c) in enumerate(zip(slopes, cis)):
            if not np.isnan(s) and not np.isnan(c):
                ax.errorbar(offsets[j], s, yerr=c, fmt='none',
                            ecolor='#333333', elinewidth=1.2,
                            capsize=3, capthick=1.2, zorder=4)
 
    ax.axhline(0, color='#555555', linewidth=0.8, linestyle='--', zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(regions, fontsize=9.5, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(title, fontsize=12, fontweight='bold', pad=10)
    ax.yaxis.grid(True, linestyle=':', alpha=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(title=f'{TCW_COL} bin', fontsize=8.5, title_fontsize=9,
              framealpha=0.9, edgecolor='#cccccc', loc='upper right')
 
    fig.tight_layout()
    fig.savefig(outpath, dpi=180, bbox_inches='tight')
    print(f"Saved -> {outpath}")
    plt.close(fig)
 
 
def make_plots(data, outdir='.'):
    """Pass your regions dict and output directory."""
    plot_shear_scaling(
        data,
        response_col=AREA_COL,
        ylabel=f'Shear scaling slope  ({AREA_COL} / {SHEAR_COL} unit)',
        title=f'Shear scaling of storm {AREA_COL} by {TCW_COL} decile bin',
        outpath=f'{outdir}/shear_scaling_{AREA_COL}.png',
    )
    plot_shear_scaling(
        data,
        response_col=PRECIP_COL,
        ylabel=f'Shear scaling slope  ({PRECIP_COL} / {SHEAR_COL} unit)',
        title=f'Shear scaling of storm {PRECIP_COL} by {TCW_COL} decile bin',
        outpath=f'{outdir}/shear_scaling_{PRECIP_COL}.png',
    )
 
 
# ── DEMO WITH SYNTHETIC DATA (remove when using real data) ───────────────────
if __name__ == '__main__':
 
    DEMO_REGIONS = ['west africa', 'india', 'south america',
                    'south africa', 'australia', 'us plains', 'east asia']
 
    def synthetic_region(seed, beta1, beta2, beta3, n=2000):
        rng = np.random.default_rng(seed)
        tcwv  = rng.uniform(20, 70, n)
        shear = rng.uniform(0, 20, n)
        area  = (beta1 * tcwv + beta2 * shear
                 + beta3 * tcwv * shear
                 + rng.normal(0, 50, n))
        prcp  = (0.8 * beta1 * tcwv + 0.6 * beta2 * shear
                 + 0.5 * beta3 * tcwv * shear
                 + rng.normal(0, 30, n))
        # column names taken directly from CONFIG — never hardcoded
        return pd.DataFrame({
            SHEAR_COL:  shear,
            TCW_COL:    tcwv,
            AREA_COL:   area,
            PRECIP_COL: prcp,
        })
 
    params = [
        (0, 2.0, 1.5, 0.04),
        (1, 1.8, 1.2, 0.03),
        (2, 2.2, 1.8, 0.05),
        (3, 1.5, 1.0, 0.02),
        (4, 2.5, 2.0, 0.06),
        (5, 1.7, 1.3, 0.035),
        (6, 2.1, 1.6, 0.045),
    ]
 
    demo_data = {
        region: synthetic_region(seed, b1, b2, b3)
        for region, (seed, b1, b2, b3) in zip(DEMO_REGIONS, params)
    }
 
    make_plots(demo_data, outdir='/mnt/user-data/outputs')
 
