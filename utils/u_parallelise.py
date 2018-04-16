import multiprocessing
import numpy as np
import pdb


def run_arrays(nb_processors, func, data, dic_names):

    pool = multiprocessing.Pool(processes=nb_processors)

    res = pool.map(func, data)
    pool.close()

    res = [x for x in res if x is not None]
    dic = {}

    res = np.array(res)

    for id, l in enumerate(dic_names):

            dic[l] = np.squeeze(res[:,id,...])

    return dic

def run_mixed(nb_processors, func, data, dic_names):

    pool = multiprocessing.Pool(processes=nb_processors)

    res = pool.map(func, data)
    pool.close()
    pdb.set_trace()
    if len(res[0]) != len(dic_names):
        print('Error with creating output dic. Need same number of function output as dictionary names')

    res = [x for x in res if x is not None]
    dic = {}
    for k in dic_names:
        dic[k] = []

    for r in res:
        for id, l in enumerate(dic_names):
            dic[l].append(np.squeeze(r[id]))

    return dic

