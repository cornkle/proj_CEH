
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl







dic = pkl.load(open('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000_blobs.p', 'rb'))

dic2 = pkl.load(open('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000_T.p', 'rb'))

ids = np.array(dic['id'])
scales = np.array(dic['scale'])

uids, uinds = np.unique(dic['id'], return_index=True)

udscale = np.unique(scales)

print(np.percentile(scales, np.arange(0, 101, 20)))

psum = np.array(dic['circle_p'])
tmin = np.array(dic['circle_t'])
pbulk_g30 = np.array(dic['bulk_g30'])
val = np.array(dic['circle_val'])

pp = np.concatenate(psum)
tt = np.concatenate(tmin)
pall_g30 = np.sum(pp > 30)

pp15 = np.concatenate(psum[scales <= 30])

pt15 = (pp[(tt <= -65)])

pp = np.concatenate(psum[scales <= 30])
tt = np.concatenate(tmin[scales <= 30])
pts15 = (pp[(tt <= -65)])

print('Nb 30mm identified', pall_g30)
print('Nb 30mm bulk', np.sum(pbulk_g30[uinds]))
print('Nb 30mm identified to bulk', pall_g30 / np.sum(pbulk_g30[uinds]))
print('Nb 30mm identified lt 40km', np.sum(pp15 >= 30))
print('Nb 30mm identified lt 40km to identified', np.sum(pp15 >= 30) / pall_g30)
print('Nb 30mm pixel identified lt 40km to bulk', np.sum(pp15 >= 30) / np.sum(pbulk_g30[uinds]))

print('Nb 30mm pixel identified lt T-65 to identified', np.sum(pt15 >= 30) / pall_g30)

print ('Fraction to valid 40scale', np.sum(pp15 >= 30) / np.sum(pp15 >= 0))

print ('Fraction to valid 40scaleT', np.sum(pts15 >= 30) / np.sum(pts15 >= 0))



tsum = np.array(dic2['circle_p'])
tmin = np.array(dic2['circle_t'])

pp = np.concatenate(tsum)
tt = np.concatenate(tmin)
ppall_g30 = np.sum(pp > 30)

pt15 = (pp[tt <= -65])

print('T 30mm identified', ppall_g30)
print('T 30mm pixel identified lt -65 to to identified', np.sum(pt15 >= 30) / ppall_g30)
print ('Fraction to valid T65', np.sum(pt15 >= 30) / np.sum(pt15 >= 0))