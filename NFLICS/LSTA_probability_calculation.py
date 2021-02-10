import numpy as np
import pandas as pd
import os
from utils import constants as cnst

pfile_path = cnst.elements_drive + '/Africa/WestAfrica/NFLICS/tables/prob_dictionary/'

def fit_curve_poly(x, y, h):
    # lower order polyfit for nighttime/morning curves to reduce noise
    if (h <= 10) | (h >= 23):
        order = 2
    else:
        order = 6
    z = np.polyfit(x, y, order)
    f = np.poly1d(z)

    return f


def run(h, LSTA_array, table_path):
    """

    :param h: Hour to calculate probability for
    :param LSTA_array: Any LSTA array (in degC!), numpy array or data array (2D)
    :param table_path: directory path to pre-calculated hourly LSTA probability tables
    :return: Array with probability factors of LSTA-array shape
    """

    tab = pd.read_table(table_path + os.sep + 'BigDomain_LSTAprobability_17-0W_9-20N_' + str(h).zfill(2) + '.csv', delimiter=',')

    f = fit_curve_poly(tab['bins'], tab['factor'], h)
    probability_values = f(LSTA_array)

    # Probability adjustment only valid within reference LSTA ranges. Extremes are set to last valid value.
    p10, p90 = np.percentile(tab['bins'],[10,90])
    probability_values[(LSTA_array>=p90)] = f(p90)
    probability_values[(LSTA_array<=p10)]= f(p10)

    return probability_values

