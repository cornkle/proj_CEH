import numpy as np
import pandas as pd


### PROJECT DOMAINS
REGIONS = {
 'GPlains' : [[-110,-90,31,48], 'nam', -6, (5,9), (1,12),1], 
 'Eur' : [[-10,28,38,53], 'nam', 1, (5,9), (1,12),1], 
 'china' : [[105,120,20,45], 'asia', 8 , (5,9), (1,12),1],
 'india' : [[70,90, 5,30], 'asia', 5, (5,9), (1,12), -1], 
 'WAf' : [[-18,25,4,25], 'spac', 0, (5,9), (1,12), -1],
 'australia' : [[120,140,-25, -11], 'asia', 9,  (11,3), (1,12), -1], 
 'SAf' : [[17,33, -32,-14], 'spac', 2, (11,3), (1,12), 1], 
 'sub_SA' : [[-68,-45, -40, -20], 'spac', -4, (11,3), (1,12),1] , 
 'trop_SA' : [[-75, -45, -20, 0], 'spac', -4, (11,3), (1,12), -1],
 'CAf'     : [[13.5, 29, -14, 4], 'spac', 1, (11,3), (1,12), -1],
 'EAf'     : [[29,44,-2.5, 15.5], 'spac', 2, (3,7), (1,12), -1],
 'Atl'     : [[-46,-18,4, 25], 'spac', -2, (6,10), (1,12), -1],
 'Pcf'     : [[-135,-110,-15, 20], 'spac', -8,  (7,11), (1,12), 1],
 'InO'     : [[55,82,-15, 5], 'spac', 5, (11,3), (1,12), -1],
}

# config/plotting.py

REGION_STYLE = {

    # ======================================================
    # Tropical continental
    # warm / earthy
    # ======================================================

    "WAf": {
        "color": "firebrick",
        "ls": "-",
        "label": "West Africa",
    },

    "CAf": {
        "color": "darkorange",
        "ls": "-",
        "label": "Central Africa",
    },

    "EAf": {
        "color": "goldenrod",
        "ls": "-",
        "label": "East Africa",
    },

    "india": {
        "color": "saddlebrown",
        "ls": "-",
        "label": "India",
    },

    "trop_SA": {
        "color": "olivedrab",
        "ls": "-",
        "label": "Tropical South America",

    },

    "australia": {
        "color": "peru",
        "ls": "-",
        "label": "Australia",
    },


    # ======================================================
    # Midlatitude continental
    # blues / purples / dark neutrals
    # ======================================================

    "GPlains": {
        "color": "royalblue",
        "ls": "-",
        "label": "Great Plains",
    },

    "Eur": {
        "color": "teal",
        "ls": "-",
        "label": "Europe",
    },

    "china": {
        "color": "navy",
        "ls": "-",
        "label": "China",
    },

    "sub_SA": {
        "color": "purple",
        "ls": "-",
        "label": "Subtropical South America",
    },

    "SAf": {
        "color": "black",
        "ls": "-",
        "label": "South Africa",
    },


    # ======================================================
    # Oceanic
    # dashed + lighter/cyan family
    # ======================================================

    "Atl": {
        "color": "deepskyblue",
        "ls": "--",
        "label": "Atlantic",
    },

    "Pcf": {
        "color": "turquoise",
        "ls": "--",
        "label": "Pacific",
    },

    "InO": {
        "color": "mediumseagreen",
        "ls": "--",
        "label": "Indian Ocean",
    },


}


def LT_to_UTC_hour(lt_hour, region):
 """
 Local hour to UTC hour, usable without date. Does not provide date.
 :param lt_hour: local time hour
 :param region: region of interest
 :return: UTC hour corresponding to local time hour
 """
 h = lt_hour - (REGIONS[region])[2]

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
 h = utc_hour + (REGIONS[region])[2]

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
 hourchange = (REGIONS[region])[2]

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
 hourchange = (REGIONS[region])[2]

 if isinstance(utc_date, str):
  utc_date = pd.to_datetime(utc_date)

 if hourchange < 0:
  date = utc_date - pd.Timedelta(str(np.abs(hourchange)) + ' hours')
 elif hourchange > 0:
  date = utc_date + pd.Timedelta(str(np.abs(hourchange)) + ' hours')
 else:
  date = utc_date
 return date




