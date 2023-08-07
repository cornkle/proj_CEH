import multiprocessing
import numpy as np
import ipdb


def multi_test():

    def f(x):
        return x + 1

    pool = multiprocessing.get_context('spawn').Pool(processes=1)
    res = pool.map(f, np.arange(10))
    pool.close()
    print(res)


def run_arrays(nb_processors, func, data, dic_names):

    pool = multiprocessing.get_context('fork').Pool(processes=nb_processors)
    print('Pooled')
    res = pool.map(func, data)
    pool.close()

    res = [x for x in res if x is not None]


    # res = []
    # for d in data:
    #     r = func(d)
    #     res.append(r)

    dic = {}

    res = np.array(res)

    for id, l in enumerate(dic_names):

            dic[l] = np.squeeze(res[:,id,...])

    return dic


def era_run_arrays(nb_processors, func, data):

    pool = multiprocessing.Pool(processes=nb_processors)

    res = pool.map(func, data)
    pool.close()

    print('Returned from parallel')
    #ipdb.set_trace()
    res = [x for x in res if x is not None]
    dic = {}

    rres = []
    #ipdb.set_trace()
    dic_names = (res[0])[1]
    for r in res:
        rres.append(np.array(r[0]))


    vars = np.array(rres)
    for id, l in enumerate(dic_names):
        try:
            dic[l] = dic[l] + np.nansum(np.squeeze(vars[:,id,...]), axis=0)
        except KeyError:
            dic[l] = np.nansum(np.squeeze(vars[:,id,...]), axis=0)

    return dic


def run_mixed(nb_processors, func, data, dic_names):

    pool = multiprocessing.Pool(processes=nb_processors)

    res = pool.map(func, data)
    pool.close()

    if len(res[0]) != len(dic_names):
        print('Error with creating output dic. Need same number of function output as dictionary names')

    res = [x for x in res if x is not None]
    dic = {}
    for k in dic_names:
        dic[k] = []

    for r in res:
        for id, l in enumerate(dic_names):
            dic[l].append(r[id])

    return dic


def run_flat(nb_processors, func, data, dic_names):

    pool = multiprocessing.Pool(processes=nb_processors)

    res = pool.map(func, data)
    pool.close()

    res = [x for x in res if x is not None]
    dic = {}

    res = np.array(res)

    for id, l in enumerate(dic_names):

            dic[l] = np.squeeze(res[:,id,...])

    return dic
