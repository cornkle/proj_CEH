import numpy as np
import matplotlib.pyplot as plt
from COCOON import shear_scaling_plots_core_v2 as core

TERCILE_BINS = [
    ('Low\n(0–20th)',   0,  20),
    ('Mid\n(40–60th)', 40,  60),
    ('High\n(80–100th)', 80, 100),
]

COLORS = {
    'Low\n(0–20th)':    '#4393c3',
    'Mid\n(40–60th)':   '#f4a582',
    'High\n(80–100th)': '#d6604d',
}


def plot_fractional_bars(regions, var1='shear', var2='tcwv', target_var='area'):
    """
    Grouped bar chart: x = regions, grouped bars = low/mid/high moisture bin.
    Y axis = fractional shear sensitivity above lowest bin (lowest = 0).
    Error bars = 95% CI.
    """
    df = core.compute_regional(regions, var2, var1, target_var, TERCILE_BINS)
    df = df.dropna(subset=['frac_slope'])

    region_list = list(regions.keys())
    bin_labels  = [b[0] for b in TERCILE_BINS]
    n_regions   = len(region_list)
    bar_width   = 0.22
    x           = np.arange(n_regions) * (len(TERCILE_BINS) * bar_width + 0.25)

    fig, ax = plt.subplots(figsize=(max(9, n_regions * 1.6), 5.5))

    for i, label in enumerate(bin_labels):
        sub     = df[df['bin_label'] == label].set_index('region')
        slopes  = [sub.loc[r, 'frac_slope'] if r in sub.index else np.nan
                   for r in region_list]
        cis     = [sub.loc[r, 'frac_ci']    if r in sub.index else np.nan
                   for r in region_list]
        offsets = x + (i - 1) * bar_width

        ax.bar(offsets, slopes, width=bar_width,
               color=COLORS[label], label=label,
               edgecolor='white', linewidth=0.6, zorder=3)

        for j, (s, c) in enumerate(zip(slopes, cis)):
            if np.isfinite(s) and np.isfinite(c):
                ax.errorbar(offsets[j], s, yerr=c, fmt='none',
                            ecolor='#333333', elinewidth=1.2,
                            capsize=3, zorder=4)

    ax.axhline(0, color='#555555', linewidth=0.8, linestyle='--', zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(region_list, fontsize=9, fontweight='bold')
    ax.set_ylabel(f'Fractional shear sensitivity  (lowest bin = 0)', fontsize=10)
    ax.set_title(f'Shear amplification by {var2} tercile  —  {target_var}\n'
                 f'(error bars = 95% CI)',
                 fontsize=10, fontweight='bold')
    ax.legend(title=f'{var2} bin', fontsize=8.5, title_fontsize=9,
              framealpha=0.9, edgecolor='#cccccc')
    ax.spines[['top', 'right']].set_visible(False)
    ax.yaxis.grid(True, linestyle=':', alpha=0.45)
    ax.set_axisbelow(True)
    plt.tight_layout()
    fname = f'fractional_bars_{var2}_{target_var}.png'
    plt.savefig(fname, dpi=180, bbox_inches='tight')
    plt.show()
    print(f"Saved -> {fname}")
