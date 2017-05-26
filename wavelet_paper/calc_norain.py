import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import pdb


dic = pkl.load( open ('/users/global/cornkle/C_paper/wavelet/saves/bulk_40big_zR.p', 'rb'))
dic2 = pkl.load( open ('/users/global/cornkle/C_paper/wavelet/saves/bulk_40big_size_zR.p', 'rb'))
scf =  pkl.load(open('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000_noR.p', 'rb'))


p = np.array(dic['pmax'])
p2 = np.array(dic2['pmax'])

p30 = np.array(dic['po30'])
p302 = np.array(dic2['po30'])
pconv = np.array(dic2['pc'])

pconv = np.concatenate(pconv)

f = plt.figure()
plt.hist(pconv, bins=np.arange(5,50,3))

print('Convective rain', np.mean(pconv[pconv>8]))
print('Convective rain', np.min(pconv[pconv>8]))

pp = np.concatenate(np.array(dic['p']))
pp2 = np.concatenate(np.array(dic2['p']))



a = np.array(dic['area'])*25
a2 = np.array(dic2['area'])*25

isfin = np.array(dic['isfin'])
isfin2 = np.array(dic2['isfin'])
isnz = np.array(dic['isnz'])
isnz2 = np.array(dic2['isnz'])



hours = np.array(dic2['hour'])

ids = np.array(scf['id'])
area_scf = np.array(scf['area'])
tmin = np.array(scf['circle_Tcentre'])
isfin_scf = np.array(scf['circle_val'])
isnz_scf = np.array(scf['circle_nz'])
p_scf = np.array(scf['circle_g30'])
pmax_scf = np.array(scf['bulk_pmax'])
pcmax_scf = np.array(scf['circle_max'])
pp_scf = np.concatenate(np.array(scf['circle_p']))
ppc_scf = np.concatenate(np.array(scf['circle_pc']))

print('Convective rain prob', np.sum(pconv>8)/np.sum(pp2>=0))
print('Convective SCF rain prob', np.sum(ppc_scf>8)/np.sum(pp_scf>=0))

uni,idd = np.unique(ids,return_index=True)

aarea_scf = np.sum(area_scf[idd])*25

valid =np.sum(np.isfinite(p))
rain = np.sum(p>0.1)

valid2 = np.sum(np.isfinite(p2))
rain2 = np.sum(p2>0.1)

ext = np.sum(p2>=30)

print('Small', rain/valid)
print('Big ones', rain2/valid2)
print ('SCF', np.sum(pmax_scf[idd]>0.1)/ np.unique(ids).size)

print('Small area', np.min(a), np.max(a))
print('Big area', np.min(a2), np.max(a2))

print('Rainy pix small', np.nansum(isnz)/np.nansum(isfin))
print('Rainy pix big', np.nansum(isnz2)/np.nansum(isfin2))
print('Rainy pix SCF', np.nansum(isnz_scf)/np.nansum(isfin_scf))

print('Mean rain small', np.nanmean(pp[pp>0.1]))
print('Mean rain big', np.nanmean(pp2[pp2>0.1]))
print('Mean rain SCF', np.nanmean(pp_scf[pp_scf>0.1]))

print('Max rain small', np.nanmean(p))
print('Max rain big', np.nanmean(p2))
print('Max rain SCF', np.nanmean(pcmax_scf))

print('30/val small', np.nansum(p30)/np.nansum(isfin))
print('30/val big', np.sum(pp2>=30)/np.sum(pp2>=0))
print('30/val SCF', np.sum(pp_scf>=30)/np.sum(pp_scf>=0))

print('Cloud pixels with TRMM overlap, "valid pixels"')
print('Valid pix small', np.nansum(isfin)*25)
print('Valid pix big', np.nansum(isfin2)*25)
print('Valid pix SCF', np.nansum(isfin_scf)*25)

print('Cloud are no TRMM overlap, MSG only (actual cloud size)')
print('cloud area small', np.sum(a))
print('cloud area big', np.sum(a2))
print('cloud area SCF', aarea_scf)

print('Small nb', a.size)
print('Big nb', a2.size)
print('SCF nb', np.unique(ids).size)


print('Lostpixel big', 100- (np.nansum(p302)/np.nansum(p30))*100)
print('Lostpixel SCF', 100- (np.nansum(p_scf)/np.nansum(p30))*100)
# f=plt.figure()
# plt.hist(hours,bins=np.arange(-0.5,24, 1)+0.5, edgecolor='black')
#
#
# f=plt.figure()
# plt.scatter(hours, a2)
#
# f=plt.figure()
# for h in range(0,24,1):
#     plt.scatter(h,np.percentile(a2[hours==h], 95))
# plt.show()

