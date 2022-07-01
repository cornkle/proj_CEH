import multiprocessing
import numpy as np
import ipdb

def mp():
    pool = multiprocessing.get_context('spawn').Pool(processes=2)
    res = pool.map(func, np.arange(10))
    pool.close()
    print(res)


def func(x):
    return x + 1

