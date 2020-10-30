from CCores import powerUtils as utils

def _create_dic(tc, tw, dx, dist, start, nb):

    dic = {}
    dic['Tcut'] = tc
    dic['Twav'] = tw
    dic['dx'] = dx
    dic['dist'] = dist
    dic['start'] = start
    dic['nb'] = nb

    return dic


############ Frequently used datasets
# cloud edge (determines cloud area), wavelet max T edge (determines cloud area for wavelet application & buffer zone)
# , resolution, distance between scales,start scale, number of scales

NAMES = {
    'METEOSAT5K': _create_dic(-40, -50, 5, 1 / 12., 15, 45),
    'METEOSAT5K_vera': _create_dic(-40, -50,5, 0.5, 25, 2),
    'METEOSAT5K_veraLS': _create_dic(-40, -50,5, 0.5, 25, 4),
    'METEOSAT3K_veraLS': _create_dic(-40, -50,3, 0.48, 9, 9),  #(-40, -50,3, 0.48, 9, 9)#(-40, -50,3, 0.60, 12, 6)
    'METEOSAT8K': _create_dic(-40, -50,8, 1 / 12., 24, 40),
    'METEOSAT10K': _create_dic(-40, -50, 10, 1 / 12., 30, 40),
    'GRIDSAT': _create_dic(-40, -50,8, 1 / 12., 24, 40),
}



UTILS = {
    'sum' : utils.find_power_sum,
    'ind' : utils.find_power_individual,
    'nflics' : utils.find_power_nflics,
    'nflicsv2' : utils.find_power_nflicsv2,
    'dominant' : utils.find_power_dominant
}