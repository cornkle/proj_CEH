# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 10:04:08 2016

@author: cornkle
"""

from datetime import datetime, timedelta
import numpy as np

# Transforms seconds of day to hours, minutes, seconds
def sec_to_time(isec):
    sec = timedelta(seconds=int(isec))
    d = datetime(1,1,1) + sec

    #print("DAYS:HOURS:MIN:SEC")
    #print("%d:%d:%d:%d" % (d.day-1, d.hour, d.minute, d.second))
    return d
    


class date_list(object):
    
    def __init__(self):
        self.y = []
        self.m = []
        self.d = []
        self.h = []
        self.mi = []
        self.s = []
        
    def add(self,y,mon,d,h,m,s):
        self.y.append(y)
        self.m.append(mon)
        self.d.append(d)
        self.h.append(h)
        self.mi.append(m)
        self.s.append(s) 
        
    def getObj(self, ind):
        
        x = date_list()
        
        for i in ind:
           print(ind)
           x.add(self.y[i] , self.m[i] , self.d[i] , self.h[i] , self.mi[i] , self.s[i])
           
        return x 
    
       
    def getStr(self, ind= None):
        
        if ind is None: ind=(list(range(len(self.y))))
        x = []
      
        for i in ind:
           
           x.append(str(self.y[i])+str(self.m[i]).zfill(2)+str(self.d[i]).zfill(2)+' '+ str(self.h[i]).zfill(2)+':'+str(self.mi[i]).zfill(2)+':'+str(self.s[i]).zfill(2))
           
        return x   
        
    def getInd(self, yr, mo, dy, hr, mi):
        
       yr=np.atleast_1d(yr) 
       mo=np.atleast_1d(mo) 
       dy=np.atleast_1d(dy) 
       hr=np.atleast_1d(hr)
       mi=np.atleast_1d(mi) 
        
       sstr=self.getStr()
       
       b=[]
        
       for _yr, _mo, _dy, _hr, _mi in zip(yr, mo, dy, hr, mi):

         dstr=str(_yr)+str(_mo).zfill(2)+str(_dy).zfill(2)+' '+ str(_hr).zfill(2)+':'+str(_mi).zfill(2)+':'+'00'          
         
         try:
            bb = sstr.index(dstr)      
         except ValueError:
            print('Date was not found') 
            bb = False 
         b.append(bb)      
       return b   
               
               
        
       
        
        
      
