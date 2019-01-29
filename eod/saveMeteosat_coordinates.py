
from eod import msg
import pdb
from utils import constants as cnst
import pickle as pkl


ext_drive = '/media/ck/Seagate/DIR/'
local_data = ext_drive + 'mymachine/'
network_data = ext_drive + 'cornkle/'

msg_folder = network_data + 'data/OBS/meteosat_WA30'
#msg_folder
msg_folder

m = msg.ReadMsg(msg_folder, y1=2013, y2=2014, months=[8])
mdic = m.read_data(m.fpath[0], llbox=[-18, 12, 1, 17])  #[-14, 2.5, 4, 11.5]

lat = mdic.lat
lon = mdic.lon

msg_ll = {'lat': lat, 'lon' : lon}

pkl.dump(msg_ll, open(cnst.network_data + 'data/OBS/saves/VERA_msg_latlon_18W12E_1N17N.p', 'wb'))
