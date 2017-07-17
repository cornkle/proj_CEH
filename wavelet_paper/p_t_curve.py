import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import pdb



dic = pkl.load( open ('/users/global/cornkle/C_paper/wavelet/saves/bulk_40big_size_zR.p', 'rb'))
scf =  pkl.load(open('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000_noR.p', 'rb'))

hour = np.array(dic['hour'])
p = np.array(dic['p'])
t = np.array(dic['t'])

p = np.concatenate(p[hour==18])

t = np.concatenate(t[hour==18])

hours = np.array(scf['hour'])
pscf = np.array(scf['circle_p'])
tscf = np.array(scf['circle_t'])
pscf = np.concatenate(pscf[hours==18])
tscf = np.concatenate(tscf[hours==18])


trange = np.arange(-90,-40,1)

xtick = trange[0:-1]+0.5

tlist = []
tscflist = []

for id, ti in enumerate(trange):

    if id == 0:
        continue

    pm = np.nanmean(p[(t>trange[id-1]) & (t <= ti) & (p>=0.1)])
    pscfm = np.nanmean(pscf[(tscf > trange[id - 1]) & (tscf <= ti) & (pscf>=0.1)])

    tlist.append(pm)
    tscflist.append(pscfm)


f = plt.figure()
plt.plot(xtick, tlist)

f = plt.figure()
plt.plot(xtick, tscflist)

