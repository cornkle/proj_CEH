import numpy as np
import pandas as pd
from scipy import stats

CLIP_PERCENTILES = (1, 99)
MIN_BIN_N = 30

def preprocess(df, moisture_col, shear_col, target_col):
    """Clip moisture outliers, return clean copy with required columns."""
    cols = [moisture_col, shear_col, target_col]
    d = df[cols].dropna().copy()
    lo = np.percentile(d[moisture_col], CLIP_PERCENTILES[0])
    hi = np.percentile(d[moisture_col], CLIP_PERCENTILES[1])
    return d[d[moisture_col].between(lo, hi)].copy()


def ols_slope(x, y):
    """OLS slope with 95% CI and p-value. Returns (nan,nan,nan) if too few points."""
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < MIN_BIN_N:
        return np.nan, np.nan, np.nan
    slope, _, _, p, se = stats.linregress(x[mask], y[mask])
    ci = stats.t.ppf(0.975, df=mask.sum() - 2) * se
    return slope, ci, p


def bin_slopes(df, moisture_col, shear_col, target_col, bins):
    """
    Compute OLS shear slope within each moisture bin.
    bins: list of (label, lo_pct, hi_pct)
    Returns DataFrame with slope, ci, p, mean_moisture, n per bin.
    """
    records = []
    for i, (label, lo, hi) in enumerate(bins):
        lo_val = np.percentile(df[moisture_col], lo)
        hi_val = np.percentile(df[moisture_col], hi)
        sub    = df[df[moisture_col].between(lo_val, hi_val)]
        slope, ci, p = ols_slope(sub[shear_col].values, sub[target_col].values)
        records.append({
            'bin_label':    label,
            'bin_order':    i,
            'mean_moisture': sub[moisture_col].mean(),
            'slope':        slope,
            'ci':           ci,
            'p':            p,
            'n':            len(sub),
        })
    return pd.DataFrame(records)


def to_fractional(df):
    """
    Convert absolute slopes to fractional increase above lowest bin (= 0).
    Operates on a DataFrame returned by bin_slopes.
    """
    df  = df.sort_values('bin_order').copy()
    ref = df.iloc[0]['slope']
    if ref == 0 or np.isnan(ref):
        df['frac_slope'] = np.nan
        df['frac_ci']    = np.nan
    else:
        df['frac_slope'] = (df['slope'] - ref) / abs(ref)
        df['frac_ci']    = df['ci'] / abs(ref)
    return df


def compute_regional(regions, moisture_col, shear_col, target_col, bins):
    """
    Run bin_slopes + to_fractional for every region.
    Returns long-form DataFrame with region column added.
    """
    records = []
    for region, df in regions.items():
        d    = preprocess(df, moisture_col, shear_col, target_col)
        sl   = bin_slopes(d, moisture_col, shear_col, target_col, bins)
        sl   = to_fractional(sl)
        sl['region'] = region
        records.append(sl)
    return pd.concat(records, ignore_index=True)


def compute_allregion(regions, moisture_col, shear_col, target_col,
                      n_bins=8, exclude_regions=None):
    """
    Pool all regions, derive absolute moisture bins from pooled data range.
    Returns (full_df, loo_dict, edges) where loo_dict maps region->DataFrame.
    """
    if exclude_regions is None:
        exclude_regions = []

    active = {k: v for k, v in regions.items() if k not in exclude_regions}

    pooled = pd.concat(
        [preprocess(df, moisture_col, shear_col, target_col)
         for df in active.values()],
        ignore_index=True
    )

    if moisture_col not in pooled.columns:
        raise ValueError(f"'{moisture_col}' not found. "
                         f"Columns: {pooled.columns.tolist()}")

    m_lo  = np.percentile(pooled[moisture_col], 1)
    m_hi  = np.percentile(pooled[moisture_col], 99)
    edges = np.linspace(m_lo, m_hi, n_bins + 1)

    def _slope_from_pool(data):
        records = []
        for i in range(len(edges) - 1):
            sub = data[data[moisture_col].between(edges[i], edges[i+1])]
            slope, ci, p = ols_slope(sub[shear_col].values,
                                     sub[target_col].values)
            records.append({
                'bin_label':    f'{edges[i]:.1f}–{edges[i+1]:.1f}',
                'bin_order':    i,
                'mean_moisture': sub[moisture_col].mean() if len(sub) > 0 else np.nan,
                'slope':        slope,
                'ci':           ci,
                'p':            p,
                'n':            len(sub),
            })
        return to_fractional(pd.DataFrame(records))

    full_df  = _slope_from_pool(pooled)
    loo_dict = {}
    for region in active:
        loo_data = pd.concat(
            [preprocess(df, moisture_col, shear_col, target_col)
             for k, df in active.items() if k != region],
            ignore_index=True
        )
        loo_dict[region] = _slope_from_pool(loo_data)

    return full_df, loo_dict, edges
