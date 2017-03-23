from wavelet import wavelet_composite_no_overlap as wcno
from wavelet import wavelet_composite_Tonly as wct
import numpy as np
from utils import u_arrays as ua
import pdb





files = ua.locate(".nc", '/users/global/cornkle/MCSfiles/WA15_big_-40_15W-20E_size/')

for f in files:
    no = wcno.file_loop(f)
    t = wct.file_loop(f)

    pdb.set_trace()

    if (no == []) & (t == []):
        continue

    print('Area T, no', t[0][11], no[0][11])
    print('Max T, no', t[0][12], no[0][12])
    print('Bulkg30 T, no', t[0][17], no[0][17])

    c = 0
    for s in np.arange(len(t)): c+=t[s][24]
    print('All g30 T:', c)
    k = 0
    for s in np.arange(len(no)): k+=no[s][25]
    print('All g30 no:', k)


    if c != k:
        pdb.set_trace()

#yp, xp = np.where(no > 30)
#f = plt.figure()
        # plt.imshow(labels_blob)
        # plt.plot(xp, yp, 'yo', markersize=3)
        # plt.title(str(sc-1))
        #
        # f = plt.figure()
        # plt.imshow(outint)
        # plt.plot(xp, yp, 'yo', markersize=3)
        # plt.title(str(sc-1))
        #
        # f = plt.figure()
        # plt.imshow(figure)
        # plt.plot(xp, yp, 'yo', markersize=3)
        # plt.title(str(sc-1))