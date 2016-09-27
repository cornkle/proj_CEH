# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:20:36 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
from utils import u_arrays as ua
import matplotlib.pyplot as plt
from scipy import ndimage
import pickle as pkl
from scipy.stats import pearsonr

def readMCS_label():
             
    files = ua.locate(".nc", '/users/global/cornkle/MCSfiles/')
    
    wave={}
    
    strarr=['0', '1', '2', '3', 'x']
       
    wave['porig']=[]

    
    for st in strarr:    
        wave['torig'+st]=[]                  
        wave['p'+st+'_max']=[]    # max p  in radius       
        wave['p'+st+'_mean']=[]    # mean p in radius
        wave['t'+st+'_mean']=[]     # t mean in radius   
        wave['t'+st+'_min']=[]     # t min in radius   
        wave['scales'+st]=[] 
        wave['pnb'+st]=[] 
        wave['xt'+st]=[] 
        wave['yt'+st]=[] 
        wave['zt'+st]=[] 
        wave['xp'+st]=[] 
        wave['yp'+st]=[] 
        wave['zp'+st]=[]        
                                     
    cntmax=0 
    cntin=0      
       
    for f in files[270:274]:
        print('Doing file: '+f)
        dic = xr.open_dataset(f)
        
     #   if (dic['time.hour'].values<15) or (dic['time.hour'].values>21):
          # print('smaller') 
     #      continue       
                
        outp=np.array(dic['p'].values.copy())   
        for nb in range(5):
                  boole=np.isnan(outp)
                  outp[boole]=-1000
                  gg=np.gradient(outp)
                  outp[boole]=np.nan
                  outp[abs(gg[1])>300]=np.nan
                  outp[abs(gg[0])>300]=np.nan               
        wave['porig'].append(outp)
        
                
        for strr in strarr[1:2]:
            
            outt=np.array(dic['t_lag'+strr].values.copy())
            outp[np.isnan(outp)]=-10**-5
            wave['torig'+strr].append(outt)    
                        
          #  outt[np.isnan(outt)]=-40
         #   outt[outt>-40]=-40
            o2 = outt.copy()                            
         #   o2[np.where(dic['pmask'].values==0)]=np.nan   
            gradi=np.gradient(o2)  
            grad=gradi[1]
                        
         #   if strr=='0':
         #     outp[np.where(dic['pmask'].values==0)]=np.nan  
            
            xt = []
            yt = []
            zt = []
            
            xp = []
            yp = []
            zp = []
            
            xg = []
            yg = []
            zg = []
            
            for ss in [8]:
                maxoutt = (outt == ndimage.minimum_filter(outt,ss, mode='constant',cval=np.amin(outt)+1))
                maxoutt = maxoutt.astype(int)
                ypks,xpks=np.where((maxoutt==1) & (outt < -50))
                zpks=[ss]*len(ypks)
                yt.extend(ypks)
                xt.extend(xpks)
                zt.extend(zpks)
            for ss in [8]:    
                maxoutp = (outp == ndimage.maximum_filter(outp,ss, mode='constant',cval=np.amax(outp)+1))
                maxoutp = maxoutp.astype(int)
                ypks,xpks=np.where((maxoutp==1) & (outp > 10))
                zpks=[ss]*len(ypks)
                yp.extend(ypks)
                xp.extend(xpks)
                zp.extend(zpks)
                
            for ss in [8]:    
                maxoutg = (grad == ndimage.minimum_filter(grad,ss, mode='constant',cval=np.amax(grad)+1))
                maxoutg = maxoutg.astype(int)
                ypks,xpks=np.where((maxoutg==1) & (grad < -10))
                zpks=[ss]*len(ypks)
                yg.extend(ypks)
                xg.extend(xpks)
                zg.extend(zpks)    

            f = plt.figure()   
            siz=5
            ax = f.add_subplot(321)
            plt.pcolormesh(outp, vmin=0, vmax=30)   
            plt.plot(xp,yp, 'yo', markersize=siz)
            plt.plot(xt,yt, 'ro', markersize=siz)        
            ax.set_title('Precip',fontsize=12)
            cbar=plt.colorbar()
            cbar.set_label('mm h-1')
    
            ax1 = f.add_subplot(322)
            plt.pcolormesh(outt, vmin=-90, vmax=-20)   
            plt.plot(xp,yp, 'yo', markersize=siz)
            plt.plot(xt,yt, 'ro', markersize=siz)        
            ax1.set_title('T',fontsize=12)
            plt.colorbar()
            
            ax1 = f.add_subplot(323)
            plt.pcolormesh(outp, vmin=0, vmax=30)   
            plt.plot(xp,yp, 'yo', markersize=siz)
      #      plt.plot(xt,yt, 'ro')        
            ax1.set_title('Precip',fontsize=12)
            plt.colorbar()
            
            ax1 = f.add_subplot(324)
            plt.pcolormesh(outt, vmin=-90, vmax=-20)   
        #    plt.plot(xp,yp, 'yo')
            plt.plot(xg,yg, 'wo', markersize=siz)     
            plt.plot(xt,yt, 'ro', markersize=siz)   
            plt.plot(xp,yp, 'yo')        
            ax1.set_title('T',fontsize=12)
            plt.colorbar()
            
            ax1 = f.add_subplot(325)
            plt.pcolormesh(outp, vmin=0, vmax=30)   
            plt.plot(xp,yp, 'yo', markersize=siz)
            plt.plot(xg,yg, 'wo')        
            ax1.set_title('Precip',fontsize=12)
            plt.colorbar()
            
            ax1 = f.add_subplot(326)
            plt.pcolormesh(grad, vmin=-20, vmax=20)   
            plt.plot(xp,yp, 'yo')
            plt.plot(xg,yg, 'wo', markersize=siz)        
            ax1.set_title('T gradient',fontsize=12)
            plt.colorbar()
          #  
         #   plt.savefig('/users/global/cornkle/figs/lmax_saves/lmax'+str(cntmax)+'.pdf')
        #    plt.close('all')
            cntmax=cntmax+1
            

#            
                
#                                       
#            for i in range(len(z)):  
#                                                
#                zz = z[i]
#                xx = x[i]
#                yy = y[i]                
#
#                                                
#                if dic['pmask'][yy,xx]==0:    # if maximum falls in region where no TRMM exists, continue                   
#                    continue
#                
#                if strr=='0':                       
#                    cntin = cntin+1
#                                                  
#                iscale = zz   
#                            
#                tmin=outt[yy,xx]     
#               
#                #Find all indices within the local circle of radius iscale...
#                # ... Then average over those indices
#                xloc1 = np.arange(xx-iscale,xx+iscale+1)
#                yloc1 = np.arange(yy-iscale,yy+iscale+1)
#                xloc,yloc = np.meshgrid(xloc1,yloc1)
#                distloc = ( (xloc-xx)**2 + (yloc-yy)**2 ) ** .5
#
#                indloc = (distloc <= iscale).nonzero()
#                ycirc = indloc[0] - iscale + yy
#                xcirc = indloc[1] - iscale + xx   
#                
#                noky=np.where(ycirc>=outt.shape[0])   # if the circle is off the edge                               
#                if noky[0].size>0:
#                    ycirc=np.delete(ycirc,noky)
#                    xcirc=np.delete(xcirc,noky)
#                    
#                nokx=np.where(xcirc>=outt.shape[1])                                  
#                if nokx[0].size>0:
#                    ycirc=np.delete(ycirc,nokx)
#                    xcirc=np.delete(xcirc,nokx)    
#                               
#                tmean=np.nanmean(dic['tc_lag'+strr].values[ycirc, xcirc])
#                pmean=np.nanmean(dic['p'].values[ycirc, xcirc])               
#                pmax=np.nanmax(outp[ycirc, xcirc])         
#                tmin=np.nanmin(outt[ycirc, xcirc])
#                pnb=ycirc.size
#              
#                wave['p'+strr+'_max'].append(pmax)
#                wave['p'+strr+'_mean'].append(pmean)
#                wave['t'+strr+'_mean'].append(tmean)
#                wave['t'+strr+'_min'].append(tmin)
#                wave['scales'+strr].append(zz*5*2)
#                wave['pnb'+strr].append(pnb)
#                               
#                #just append all the variables into a dictionary now!                                                 
#    
#    for k in wave:
#        if isinstance(wave[k][0], np.ndarray):
#            continue
#        print(k)
#        wave[k]=np.array(wave[k])
#    
#    pkl.dump(wave, open('/users/global/cornkle/MCSfiles/save/MCS_allyears_label.p', 'wb'))
#    
#    print('Saved!') 
#    print('Found '+str(cntmax)+' maxima in '+str(len(files))+' systems.')
#    print(str(cntin)+' maxima coincided with TRMM')    
    
#readMCS_getWavelet_label()
    
def plotShit_label_noW():
    
    wavelet = pkl.load( open ('/users/global/cornkle/MCSfiles/save/MCS_allyears_label.p', 'rb'))
    #plt.figure(1)
    
    strg=['0', '1', '2', '3', 'x']
    f = plt.figure()
    k=361
    cnt=0
    for strr in strg:    
    
        t=np.array(wavelet['t'+strr+'_min'])#*(-1)   
        p=np.array(wavelet['p'+strr+'_max'])
        plt100=np.where(t<-40)
        #  print(plt100[0:10])

        # m,b = np.polyfit(p,t1,1)
        r=pearsonr(p[plt100],t[plt100])
            
        title = "Tmin and Pmax (t"+strr+")"
        ax = f.add_subplot(361+cnt)
        ax.set_title(title,fontsize=12)
    
        ax.text(0.7, 0.8, 'r: '+str(np.round(r[0], decimals=2)), transform=ax.transAxes)
        ax.scatter(t[plt100], p[plt100])
        #plt.plot(p, m*p+b, '-')    
        ax.set_xlabel('T')
        ax.set_ylabel('P')
        cnt=cnt+1
    
    cnt=0
    for strr in strg:          
                  
        t=np.array(wavelet['t'+strr+'_mean'])#*(-1)   
        p=np.array(wavelet['p'+strr+'_max'])
        plt100=np.where(t<-40)
        #  print(plt100[0:10])

        # m,b = np.polyfit(p,t1,1)
        r=pearsonr(p[plt100],t[plt100])
            
        title = "Tmean and Pmax (t"+strr+")"
        ax = f.add_subplot(3,6,7+cnt)
        ax.set_title(title,fontsize=12)
    
        ax.text(0.7, 0.8, 'r: '+str(np.round(r[0], decimals=2)), transform=ax.transAxes)
        ax.scatter(t[plt100], p[plt100])
        #plt.plot(p, m*p+b, '-')    
        ax.set_xlabel('T')
        ax.set_ylabel('P')
        cnt=cnt+1
         
        
            
        
 #   plt.tight_layout()
    plt.show()          