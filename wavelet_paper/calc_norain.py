import pickle as pkl
import numpy as np



dic = pkl.load( open ('/users/global/cornkle/C_paper/wavelet/saves/bulk_40big_zR.p', 'rb'))
dic2 = pkl.load( open ('/users/global/cornkle/C_paper/wavelet/saves/bulk_40big_size_zR.p', 'rb'))


p = np.array(dic['pmax'])
p2 = np.array(dic2['pmax'])
a = np.array(dic['area'])*25
a2 = np.array(dic2['area'])*25



valid =np.sum(np.isfinite(p))
rain = np.sum(p>0.1)

valid2 = np.sum(np.isfinite(p2))
rain2 = np.sum(p2>0.1)

ext = np.sum(p2>=30)

print('Small', rain/valid)
print('Big ones', rain2/valid2)
print('Ext', ext/valid2)
print('Small area', np.min(a), np.max(a))
print('Big area', np.min(a2), np.max(a2))
print('Small nb', a.size)
print('Big nb', a2.size)

