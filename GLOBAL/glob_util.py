import sys
import os
import numpy as np
import xarray as xr
import glob
import ipdb
import pandas as pd


### PROJECT DOMAINS
MREGIONS = {'WAf' : [[-18,25,4,25], 'spac', 0, (1,7), (1,12)], # last is hourly offset to UCT
 'SAf' : [[20,35, -35,-15], 'spac', 2, (9,12), (1,12)],
 'india' : [[70,90, 5,30], 'asia', 5, (1,7), (1,12)],
 'china' : [[105,115,25,40], 'asia', 8 , (1,7), (1,12)],
 'australia' : [[120,140,-23, -11], 'asia', 9, (10,12), (1,12)],
 'sub_SA' : [[-68,-47, -40, -20.5], 'spac', -4, (9,12), (1,12)] ,
 #'trop_SA' : [[-75, -50, -20, -5], 'spac', -5, (1,12), (1,12), (1,12)],
 'GPlains' : [[-100,-90,32,47], 'nam', -6, (1,7), (1,12)]
}


def LT_to_UTC_hour(lt_hour, region):
 """
 Local hour to UTC hour, usable without date. Does not provide date.
 :param lt_hour: local time hour
 :param region: region of interest
 :return: UTC hour corresponding to local time hour
 """
 h = lt_hour - (MREGIONS[region])[2]

 if h >= 24:
  h = h - 24
 if h == 24:
  h = 0
 if h < 0:
  h = h + 24
 return h

def UTC_to_LT_hour(utc_hour, region):
 """
 UTC hour to LT hour, usable without date. Does not provide date.
 :param utc_hour: utc time hour
 :param region: region of interest
 :return: LT hour corresponding to utc time hour
 """
 h = utc_hour + (MREGIONS[region])[2]

 if h >= 24:
  h = h - 24
 if h == 24:
  h = 0
 if h < 0:
  h = h + 24
 return h


def LT_to_UTC_date(lt_date, region):
 """
 Local date to UTC date, uses pandas delta time.
 :param date: local date as datetime object (datetime.datetime), e.g. pd.to_datetime(date_string) or string
 :param region: region of interest
 :return: UTC date corresponding to LT time date
 """
 hourchange = (MREGIONS[region])[2]

 if isinstance(lt_date, str):
  lt_date = pd.to_datetime(lt_date)

 if hourchange < 0:
  date = lt_date + pd.Timedelta(str(np.abs(hourchange)) + ' hours')
 elif hourchange > 0:
  date = lt_date - pd.Timedelta(str(np.abs(hourchange)) + ' hours')
 else:
  date = lt_date
 return date



def UTC_to_LT_date(utc_date, region):
 """
 UTC date to LT date, uses pandas delta time.
 :param date: UTC date as datetime object (datetime.datetime), e.g. pd.to_datetime(date_string) or string "Y-m-d hh:mm:ss"
 :param region: region of interest
 :return: LT date corresponding to UTC time date
 """
 hourchange = (MREGIONS[region])[2]

 if isinstance(utc_date, str):
  utc_date = pd.to_datetime(utc_date)

 if hourchange < 0:
  date = utc_date - pd.Timedelta(str(np.abs(hourchange)) + ' hours')
 elif hourchange > 0:
  date = utc_date + pd.Timedelta(str(np.abs(hourchange)) + ' hours')
 else:
  date = utc_date
 return date




