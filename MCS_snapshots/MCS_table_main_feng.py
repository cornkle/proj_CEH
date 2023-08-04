from MCS_snapshots import MCS_table_create_feng
import xarray as xr
import ipdb
import pandas as pd
import glob
import os


def make_table():
    """
    Start with scanning image for MCSs as defined in MCS_table_create.process_tir_image.
    :return:
    """
    lmcs = '/prj/global_water/MCS_Feng/global_v2/2d_fields/'
    out = '/prj/global_water/MCS_5000km2_tables/australia/'
    box = [120,140,-23,-11] #EAf [26,52.5,-5.6, 15.5]
    for yy in range(2001,2021):
        infiles = sorted(glob.glob(lmcs + str(yy) + '*/*.nc'))
        full_year = []
        out_dic = {}
        print('Doing', yy)
        outfile = out + str(yy)+'_MCS_5000km2_-50C_0.1degTIR-IMERG_hourly.csv'
        if os.path.isfile(outfile):
           print(outfile, ' exists, continue') 
           continue

        for infile in infiles:
            try:
                da = (xr.open_dataset(infile))
            except:
                ipdb.set_trace()
            da = da.sel(lon=slice(box[0],box[1]), lat=slice(box[2], box[3])) 
            basic_tab = MCS_table_create_feng.process_tir_image(da, 10)
            merge_tab = MCS_table_create_feng.add_environment_toTable(basic_tab, da,  envvar_take=[], rainvar_name='precipitation')
            
            merge_tab.pop('cloudMask')
            merge_tab.pop('tir')
            if len(merge_tab['date']) ==0:
               continue

            full_year.append(merge_tab)
            print('Did' , infile)
        
        for key in full_year[0].keys():
            out_dic[key] = []
        for single_tab in full_year:
            for key in single_tab.keys():
                out_dic[key].extend(single_tab[key])
                
        pd_out = pd.DataFrame.from_dict(out_dic)
        pd_out.to_csv(outfile)
        del out_dic
        del pd_out
        del merge_tab
        del basic_tab

