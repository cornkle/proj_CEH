# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 10:04:08 2016

@author: cornkle
"""

from datetime import datetime, timedelta

# Transforms seconds of day to hours, minutes, seconds
def sec_to_time(isec):
    sec = timedelta(seconds=int(isec))
    d = datetime(1,1,1) + sec

    #print("DAYS:HOURS:MIN:SEC")
    #print("%d:%d:%d:%d" % (d.day-1, d.hour, d.minute, d.second))
    return d
    


class date_list(object):
    
    def __init__(self):
        self.years = []
        self.months = []
        self.days = []
        self.hours = []
        self.minutes = []
        self.seconds = []
        
    def add(self,y,mon,d,h,m,s):
        self.years.append(y)
        self.months.append(mon)
        self.days.append(d)
        self.hours.append(h)
        self.minutes.append(m)
        self.seconds.append(s) 
        
    def get_obj(self, ind):
        
        x = date_list()
        
        for i in ind:
           print(ind)
           x.add(self.years[i] , self.months[i] , self.days[i] , self.hours[i] , self.minutes[i] , self.seconds[i])
           
        return x 
        
    def get_str(self, ind):
        
        x = []
        
        for i in ind:
           
           x.append(str(self.years[i])+str(self.months[i]).zfill(2)+str(self.days[i]).zfill(2)+' '+ str(self.hours[i]).zfill(2)+':'+str(self.minutes[i]).zfill(2)+':'+str(self.seconds[i]).zfill(2))
           
        return x     
