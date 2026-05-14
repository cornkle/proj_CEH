import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from scipy.stats import zscore

# ── CONFIG ────────────────────────────────────────────────────────────────────
SHEAR_COL  = 'shear'
TCW_COL    = 'tcwv'
PRECIP_COL = 'prcp'
AREA_COL   = 'area'
# ─────────────────────────────────────────────────────────────────────────────

# ── FIT INTERACTION MODELS ACROSS REGIONS ────────────────────────────────────
records = []

for region, df in regions.items():
    # z-score normalise so beta3 is comparable across regions
    d = df[[TCW_COL, SHEAR_COL, AREA_COL, PRECIP_COL]].dropna().copy()
    d['tcwv_z'] = zscore(d[TCW_COL])
    d['shear_z'] = zscore(d[SHEAR_COL])

    for response, response_col in [('area', AREA_COL), ('precip', PRECIP_COL)]:
        m = smf.ols(f'{response_col} ~ tcwv_z * shear_z', data=d).fit()
        records.append({
            'region':        region,
            'response':      response,
            'beta_tcwv':     m.params['tcwv_z'],
            'beta_shear':    m.params['shear_z'],
            'beta3':         m.params['tcwv_z:shear_z'],
            'p_beta3':       m.pvalues['tcwv_z:shear_z'],
            'r2':            m.rsquared,
            'mean_tcwv':     d[TCW_COL].mean(),
        })

summary = pd.DataFrame(records)
print(summary.to_string(index=False))

# ── PLOT: beta3 per region, side by side for area and precip ─────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=False)

for ax, response in zip(axes, ['area', 'precip']):
    sub = summary[summary['response'] == response].copy()
    colors = ['#d6604d' if p < 0.05 else '#aaaaaa' for p in sub['p_beta3']]

    bars = ax.barh(sub['region'], sub['beta3'], color=colors,
                   edgecolor='white', linewidth=0.6)

    # mark significance
    for i, (val, p) in enumerate(zip(sub['beta3'], sub['p_beta3'])):
        if p < 0.05:
            ax.text(val + 0.001, i, '  *', va='center', fontsize=11,
                    color='#d6604d', fontweight='bold')

    ax.axvline(0, color='#555555', linewidth=0.8, linestyle='--')
    ax.set_title(f'Interaction β₃  (tcwv × shear → {response})',
                 fontsize=11, fontweight='bold')
    ax.set_xlabel('β₃ (normalised units)', fontsize=9)
    ax.spines[['top', 'right']].set_visible(False)
    ax.xaxis.grid(True, linestyle=':', alpha=0.5)
    ax.set_axisbelow(True)

fig.suptitle('Shear sensitivity increases with TCWV (β₃ > 0 = interaction present)',
             fontsize=12, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('interaction_beta3.png', dpi=180, bbox_inches='tight')
plt.show()

# ── BONUS: headline scatter — mean TCWV vs beta3 across regions ───────────────
fig2, axes2 = plt.subplots(1, 2, figsize=(11, 4.5))

for ax, response in zip(axes2, ['area', 'precip']):
    sub = summary[summary['response'] == response]
    ax.scatter(sub['mean_tcwv'], sub['beta3'], s=80,
               color='#4393c3', edgecolors='#1a5276', linewidths=0.8, zorder=3)

    # label each point
    for _, row in sub.iterrows():
        ax.annotate(row['region'], (row['mean_tcwv'], row['beta3']),
                    fontsize=7.5, textcoords='offset points', xytext=(4, 3))

    # regression line through the 7 region points
    if len(sub) > 2:
        m_scatter = smf.ols('beta3 ~ mean_tcwv', data=sub).fit()
        x_range = np.linspace(sub['mean_tcwv'].min(), sub['mean_tcwv'].max(), 100)
        ax.plot(x_range, m_scatter.predict(pd.DataFrame({'mean_tcwv': x_range})),
                color='#d6604d', linewidth=1.5, linestyle='--', alpha=0.8)
        r2 = m_scatter.rsquared
        ax.text(0.05, 0.92, f'R²={r2:.2f}', transform=ax.transAxes,
                fontsize=9, color='#d6604d')

    ax.axhline(0, color='#555555', linewidth=0.8, linestyle='--')
    ax.set_xlabel(f'Mean {TCW_COL} (kg/m²)', fontsize=9)
    ax.set_ylabel('β₃ interaction term', fontsize=9)
    ax.set_title(f'Moister regions → stronger shear coupling ({response})',
                 fontsize=10, fontweight='bold')
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(True, linestyle=':', alpha=0.4)

plt.tight_layout()
plt.savefig('interaction_scatter.png', dpi=180, bbox_inches='tight')
plt.show()
