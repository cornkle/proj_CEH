import os


def _create_dic(dx, dist, start, nb):

    dic = {}
    dic['dx'] = dx
    dic['dist'] = dist
    dic['start'] = start
    dic['nb'] = nb

    return dic


############ Frequently used DATASET SETUPS.
# Defines surface thresholds, data resolution and (number of) scales
# _create_dic entries as follows:
#
# (1) resolution, (2) distance between scales, (3) start scale, (4) number of scales

NAMES = {
    'SM5k': _create_dic(5, 0.25, 5, 35),  # to avoid bias in full wavelet reconstruction, (2) = 0.25 and a minimum of 35 scales is needed.
    'SM5k_sensiTest': _create_dic(5, 0.25, 5, 20),
    'CP4_VARS': _create_dic(4.4, 0.5, 12, 10),  # nb =14 also tested
    'CP4_OBS_VARS': _create_dic(3, 0.5, 12, 10),  # nb =14 also tested
    'P25_VARS': _create_dic(25, 0.5, 25, 8),  # nb =14 also tested
    'organisation' : _create_dic(1, 0.1, 5, 43),  # nb =14 also tested
}

########## Test case data

TESTDATA = os.path.abspath(os.path.dirname(__file__)) + os.sep + 'testdata' + os.sep + 'sm_testfile.nc'
