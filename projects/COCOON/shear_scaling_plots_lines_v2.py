import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from projects.COCOON import shear_scaling_plots_core_v2 as core

DECILE_BINS = [
    ('0–10th',   0,  10),
    ('10–20th',  10, 20),
    ('20–30th',  20, 30),
    ('30–40th',  30, 40),
    ('40–50th',  40, 50),
    ('50–60th',  50, 60),
    ('60–70th',  60, 70),
    ('70–80th',  70, 80),
    ('80–90th',  80, 90),
    ('90–100th', 90, 100),
]


def _region_colors(region_list):
    palette = cm.get_cmap('tab10', len(region_list))
    return {r: palette(i) for i, r in enumerate(region_list)}


def plot_absolute(regions, var1='shear', var2='tcwv', target_var='area'):
    """
    Scatter: mean moisture of bin (x) vs absolute shear slope (y).
    Each point = one region x one bin. Error bars = 95% CI.
    """
    df = core.compute_regional(regions, var2, var1, target_var, DECILE_BINS)
    df = df.dropna(subset=['slope'])

    region_list  = list(regions.keys())
    colors       = _region_colors(region_list)

    fig, ax = plt.subplots(figsize=(8, 5))

    for region in region_list:
        sub = df[df['region'] == region]
        ax.errorbar(sub['mean_moisture'], sub['slope'],
                    yerr=sub['ci'], fmt='o',
                    color=colors[region], ecolor=colors[region],
                    elinewidth=1.0, capsize=3, markersize=6,
                    markeredgecolor='white', label=region, zorder=3)
        sub_s = sub.sort_values('bin_order')
        ax.plot(sub_s['mean_moisture'], sub_s['slope'],
                color=colors[region], linewidth=0.8, alpha=0.4, zorder=2)

    ax.axhline(0, color='#aaaaaa', linewidth=0.8, linestyle=':')
    ax.set_xlabel(f'{var2} bin mean', fontsize=10)
    ax.set_ylabel(f'{target_var} per unit {var1}', fontsize=10)
    ax.set_title(f'Shear sensitivity of {target_var} vs {var2}\n'
                 f'(each point = region × bin, error bars = 95% CI)',
                 fontsize=10, fontweight='bold')
    ax.legend(fontsize=7.5, framealpha=0.9, edgecolor='#cccccc')
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(True, linestyle=':', alpha=0.4)
    plt.tight_layout()
    fname = f'absolute_{var2}_{target_var}.png'
    plt.savefig(fname, dpi=180, bbox_inches='tight')
    plt.show()
    print(f"Saved -> {fname}")


def plot_fractional_lines(regions, var1='shear', var2='tcwv', target_var='area'):
    """
    Line plot: x = decile bins, y = fractional shear sensitivity (lowest bin = 0).
    One coloured line per region. Solid = p<0.05 interaction, dashed = not significant.
    """
    df = core.compute_regional(regions, var2, var1, target_var, DECILE_BINS)
    df = df.dropna(subset=['frac_slope'])

    region_list = list(regions.keys())
    colors      = _region_colors(region_list)
    bin_x       = list(range(len(DECILE_BINS)))

    fig, ax = plt.subplots(figsize=(10, 5))

    for region in region_list:
        sub   = df[df['region'] == region].sort_values('bin_order')
        color = colors[region]
        frac  = sub['frac_slope'].values
        ci    = sub['frac_ci'].values
        sig   = (sub['p'] < 0.05).any()

        ax.fill_between(bin_x, frac - ci, frac + ci,
                        color=color, alpha=0.12, zorder=1)
        ax.plot(bin_x, frac, color=color, linewidth=1.8,
                linestyle='-' if sig else '--', alpha=0.85, zorder=2)
        ax.scatter(bin_x, frac, color=color, s=50,
                   edgecolors='white', linewidths=0.6, zorder=3)
        ax.annotate(region, xy=(bin_x[-1], frac[-1]),
                    xytext=(5, 0), textcoords='offset points',
                    fontsize=7.5, color=color, va='center')

    ax.axhline(0, color='#555555', linewidth=0.9, linestyle=':', zorder=1)
    ax.set_xticks(bin_x)
    ax.set_xticklabels([b[0] for b in DECILE_BINS], fontsize=8,
                       rotation=45, ha='right')
    ax.set_ylabel(f'Fractional shear sensitivity  (lowest bin = 0)', fontsize=10)
    ax.set_title(f'Fractional shear amplification  —  {target_var} vs {var2}\n'
                 f'(solid = p<0.05, shading = 95% CI)',
                 fontsize=10, fontweight='bold')
    handles = [plt.Line2D([0], [0], color=colors[r], linewidth=2, label=r)
               for r in region_list]
    ax.legend(handles=handles, fontsize=7.5, framealpha=0.9,
              edgecolor='#cccccc', loc='upper left')
    ax.spines[['top', 'right']].set_visible(False)
    ax.yaxis.grid(True, linestyle=':', alpha=0.4)
    ax.set_axisbelow(True)
    ax.set_xlim(-0.5, len(DECILE_BINS) - 0.2)
    plt.tight_layout()
    fname = f'fractional_lines_{var2}_{target_var}.png'
    plt.savefig(fname, dpi=180, bbox_inches='tight')
    plt.show()
    print(f"Saved -> {fname}")


def plot_allregion(regions, var1='shear', var2='tcwv', target_var='area',
                   n_bins=8, exclude_regions=None):
    """
    All-region pooled plot with leave-one-out sensitivity lines.
    Absolute moisture bins derived from pooled data range.
    """
    full, loo_dict, edges = core.compute_allregion(
        regions, var2, var1, target_var,
        n_bins=n_bins, exclude_regions=exclude_regions
    )

    region_list = list({k: None for k in regions
                        if k not in (exclude_regions or [])})
    palette     = cm.get_cmap('tab10', len(region_list))
    x           = np.arange(len(full))

    fig, ax = plt.subplots(figsize=(9, 5))

    # full combined
    ax.fill_between(x,
                    full['frac_slope'] - full['frac_ci'],
                    full['frac_slope'] + full['frac_ci'],
                    color='#333333', alpha=0.15, zorder=1)
    ax.errorbar(x, full['frac_slope'], yerr=full['frac_ci'],
                fmt='o', color='#333333', ecolor='#333333',
                elinewidth=1.2, capsize=4, markersize=7,
                markeredgecolor='white', linewidth=2.5,
                label='All regions', zorder=4)

    # leave-one-out
    for j, (region, loo_df) in enumerate(loo_dict.items()):
        ax.plot(x, loo_df['frac_slope'],
                color=palette(j), linewidth=1.2,
                linestyle='--', alpha=0.85,
                label=f'excl. {region}', zorder=2)

    ax.axhline(0, color='#aaaaaa', linewidth=0.8, linestyle=':', zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(full['bin_label'], fontsize=8, rotation=45, ha='right')
    ax.set_xlabel(f'{var2}  (absolute bins, derived from pooled data)', fontsize=10)
    ax.set_ylabel(f'Fractional shear sensitivity\n(lowest bin = 0)', fontsize=10)
    ax.set_title(f'All-region shear amplification  —  {target_var}\n'
                 f'moisture: {var2}  |  dashed = leave-one-out',
                 fontsize=10, fontweight='bold')
    ax.legend(fontsize=7.5, framealpha=0.9, edgecolor='#cccccc', loc='upper left')
    ax.spines[['top', 'right']].set_visible(False)
    ax.yaxis.grid(True, linestyle=':', alpha=0.4)
    ax.set_axisbelow(True)
    plt.tight_layout()
    excl_str = '_excl_' + '_'.join(exclude_regions) if exclude_regions else ''
    fname = f'allregion_{var2}_{target_var}{excl_str}.png'
    plt.savefig(fname, dpi=180, bbox_inches='tight')
    plt.show()
    print(f"Saved -> {fname}")
