from wavelet import wavelet_composite_no_overlap as wcno
from wavelet import wavelet_composite_Tonly as wct
import numpy as np
from utils import u_arrays as ua
import pdb
import matplotlib.pyplot as plt


files = ua.locate(".nc", '/users/global/cornkle/MCSfiles/WA15_big_-40_15W-20E_size/')

for id, f in enumerate(files[1043::]):

    no = wcno.file_loop(f)
    t = wct.file_loop(f)
    print(id+1043)

    if ((t == None) & (no == None)):
        continue

    if ((t == None) & (no != None)):
        print('T isnt')
        pdb.set_trace()

    if ((t != None) & (no == None)):
        print('No isnt')
        pdb.set_trace()
        continue

    if ((t == []) & (no == [])):
        continue

    if ((t == []) & (no != [])):
        print('T isnt')
        pdb.set_trace()

    if ((t != []) & (no == [])):
        print('No isnt')
        pdb.set_trace()
        continue


    print('got through')


    print('Area T, no', t[0][11], no[0][11])
    print('Max T, no', t[0][12], no[0][12])
    print('Bulkg30 T, no', t[0][17], no[0][17])

    c = 0
    for s in np.arange(len(t)):
        c+=t[s][24]
    print('All g30 T:', c)
    k = 0
    for s in np.arange(len(no)):
        k+=no[s][24]
    print('All g30 no:', k)

    if t[0][17] != no[0][17]:
        print('bulkg30 not same!')
        pdb.set_trace()


    if c < k:
        print('Circles not same!')
        pdb.set_trace()

    if t[0][17] != c:
        print('T bulk and circle not same!')
        pdb.set_trace()

    if np.min(t[0][19]) < 0:
        pdb.set_trace()
#yp, xp = np.where(no > 30)
#