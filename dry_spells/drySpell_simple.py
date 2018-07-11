import numpy as np



def spell(timeseries, days, type='dry', rain_thresh=1):
    """
    Calculates number of spells of a certain length, type can be 'wet' or 'dry' spell
    :param timeseries: the timeseries to check
    :param days: the number of days of the spell
    :keyword type: string stating 'dry' or 'wet' for dry or wet spell, default is dry
    :keyword type: threshold to be considered as "rainfall", default is 1mm
    :return: the number of spells over the time period

    """
    if type == 'dry':
        pos = np.where(timeseries>=rain_thresh)  # used for calculating the differences between indices
    if type == 'wet':
        pos = np.where(timeseries<rain_thresh)
    y = np.append(np.append(0, pos[0]), timeseries.size-1)

    gaps = y[1::]-y[0:-1]

    is_spell = np.sum(gaps>=days)

    return is_spell