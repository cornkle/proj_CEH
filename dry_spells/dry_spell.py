import pandas as pd
import numpy as np


class RainySeason(pd.Series):
    """ Class doc"""

    def __init__(self, data):
        """ Class initialiser """
        assert isinstance(data, pd.Series), "Exp series got {0.__class__}".format(data)
        if data.index.dtype == int:
            idx = data.index
        else:
            idx = np.arange(data.index.size) + 1

        assert idx.min() == 1 and idx.max() in (365, 366) #testing for right format, else AssertionError
        ts = data.fillna(0)
        mrr = ts.mean()
        diff = ts - mrr
        cumdif = diff.cumsum()

        pd.Series.__init__(self, cumdif.values, index=idx)