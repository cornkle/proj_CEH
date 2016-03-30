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
    

# Transforms seconds of day to hours, minutes, seconds
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