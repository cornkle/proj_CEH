import numpy as np
from numpy import ma
from matplotlib.colors import Normalize
import scipy.stats as stats
import statsmodels.api as sm
import pandas as pd
from sklearn import linear_model

class MidPointNorm(Normalize):
    def __init__(self, midpoint=0, vmin=None, vmax=None, clip=False):
        Normalize.__init__(self,vmin, vmax, clip)
        self.midpoint = midpoint

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if not (vmin < midpoint < vmax):
            raise ValueError("midpoint must be between maxvalue and minvalue.")
        elif vmin == vmax:
            result.fill(0) # Or should it be all masked? Or 0.5?
        elif vmin > vmax:
            raise ValueError("maxvalue must be bigger than minvalue")
        else:
            vmin = float(vmin)
            vmax = float(vmax)
            if clip:
                mask = ma.getmask(result)
                result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)

            # ma division is very slow; we can take a shortcut
            resdat = result.data

            #First scale to -1 to 1 range, than to from 0 to 1.
            resdat -= midpoint
            resdat[resdat>0] /= abs(vmax - midpoint)
            resdat[resdat<0] /= abs(vmin - midpoint)

            resdat /= 2.
            resdat += 0.5
            result = ma.array(resdat, mask=result.mask, copy=False)

        if is_scalar:
            result = result[0]
        return result



def histo_frequency(data, **kwargs):
    weights = np.ones_like(data) / float(len(data))
    hist, h = np.histogram(data, weights=weights, **kwargs)
    count, h = np.histogram(data, **kwargs)
    return hist, count, h


def linear_trend(x):
    pf = np.polyfit(np.arange(len(x)), x, 1)
    return pf


def fdr_threshold(pvalues, alpha=0.05):
    """Computes the FDR threshod after Wilks (2016)."""
    p = np.sort(np.asarray(pvalues).flatten())
    n = len(p)

    return np.max(np.where(p <= (np.arange(1, n+1) / n * alpha), p, 0))



def cor(x, y):
    """It is annoying that np.corrcoef returns a matrix, this returns a float."""

    return np.corrcoef(x, y)[0, 1]


def pcor(x, y, c):
    """Partial correlation (r2) of x and y when the effect of C is removed.

    Couldn't find a routine to do exactly this in statsmodels, so I
    rolled my own. (F. Maussion)

    y and x are the variables from which we want to compute the correlation,
    when the effect of the controlling variables in C is removed.

    The residuals after regressing X/Y on Ci are the parts of X/Y that
    cannot be predicted by Ci. The partial correlation coefficient between
    Y and X adjusted for Ci is the correlation between these two sets of
    residuals.

    Returns
    -------
    tuple (r2, pvalue)

    """

    # Degrees of freedom
    if len(c.shape) == 1:
        df = len(x) - 2 - 1
    else:
        df = len(x) - 2 - c.shape[1]

    # Dont forget the constant
    _c = sm.add_constant(c)
    fity = sm.OLS(y, _c).fit()
    fitx = sm.OLS(x, _c).fit()

    r = cor(fitx.resid, fity.resid)
    t = r / np.sqrt((1. - r ** 2) / df)
    p_e = stats.t.sf(np.abs(t), df) * 2  # error probability (two tailed)

    return r, p_e


def multi_partial_correlation(input_df):
    """
    Returns the sample linear partial correlation coefficients between pairs of variables,
    controlling for all other remaining variables

    Parameters
    ----------
    input_df : array-like, shape (n, p)
        Array with the different variables. Each column is taken as a variable.

    Returns
    -------
    P : array-like, shape (p, p)
        P[i, j] contains the partial correlation of input_df[:, i] and input_df[:, j]
        controlling for all other remaining variables.
    """
    partial_corr_matrix = np.zeros((input_df.shape[1], input_df.shape[1]));
    for i, column1 in enumerate(input_df):
        for j, column2 in enumerate(input_df):
            control_variables = np.delete(np.arange(input_df.shape[1]), [i, j]);
            if i==j:
                partial_corr_matrix[i, j] = 1;
                continue
            data_control_variable = input_df.iloc[:, control_variables]
            data_column1 = input_df[column1].values
            data_column2 = input_df[column2].values
            fit1 = linear_model.LinearRegression(fit_intercept=True)
            fit2 = linear_model.LinearRegression(fit_intercept=True)
            fit1.fit(data_control_variable, data_column1)
            fit2.fit(data_control_variable, data_column2)
            residual1 = data_column1 - (np.dot(data_control_variable, fit1.coef_) + fit1.intercept_)
            residual2 = data_column2 - (np.dot(data_control_variable, fit2.coef_) + fit2.intercept_)
            partial_corr_matrix[i,j] = stats.pearsonr(residual1, residual2)[0]
    return pd.DataFrame(partial_corr_matrix, columns = input_df.columns, index = input_df.columns)

# # Generating data in our minion world
# test_sample = 10000;
# Math_score = np.random.randint(100,600, size=test_sample) + 20 * np.random.random(size=test_sample)
# Eng_score = np.random.randint(100,600, size=test_sample) - 10 * Math_score + 20 * np.random.random(size=test_sample)
# Phys_score = Math_score * 5 - Eng_score + np.random.randint(100,600, size=test_sample) + 20 * np.random.random(size=test_sample)
# Econ_score = np.random.randint(100,200, size=test_sample) + 20 * np.random.random(size=test_sample)
# Hist_score = Econ_score + 100 * np.random.random(size=test_sample)
#
# minions_df = pd.DataFrame(np.vstack((Math_score, Eng_score, Phys_score, Econ_score, Hist_score)).T,
#                           columns=['Math', 'Eng', 'Phys', 'Econ', 'Hist'])
#
# calculate_partial_correlation(minions_df)
