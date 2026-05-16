import numpy as np
import pandas as pd
import datetime,os


adj=False
hist_stats=False
hist_mask=False
hist_clim=True
hist_anom=False
copy_dirs=False

nt=3 #mask any points in grid
jobdir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/jobs/"
rundir="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/code/SEVIRI_LST/"
years=np.arange(2004,2016,1)
time_req=10
nt_back=3 #mask any points in grid
nt_for=1 #mask any points in grid

dir_old="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward/"
#dir_new="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward/"
#dir_new="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward_sig2_count50/"
#dir_new="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward_sig2_count100/"
#dir_new="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward_stats_sig2_count50_mask_sig2_count50/"
dir_new="/home/hymod/seodey/scratch/hymod/hmf/NFLICS/SEVIRI_LST/historic_mask_adj_3back1forward_stats_sig3_count50_mask_2sig_count50/"

start_time=datetime.datetime(2012,12,25,7,0)
end_time=datetime.datetime(2012,12,25,19,45)
times=pd.date_range(start_time,end_time,freq="15min")
times=[str(date).replace("-","").replace("T","").replace(":","")[-6:-2] for date in times]

if copy_dirs:
    if not os.path.exists(dir_new):
       os.makedirs(dir_new)
    for year in years:
        name="3copy_SEVIRI_LST_"+str(year)
        job_file=open(jobdir+name+'.job','w')
        job_file.writelines([  '#!/bin/sh\n',
            '#PBS -N '+name+'\n',
            '#PBS -e '+name+'.error\n',
            '#PBS -o '+name+'.output\n',
            '#PBS -q r630-shortq\n',                              
            '#PBS -l nodes=1:ppn=1\n',                
            '#PBS -l walltime=1:00:00\n',
            '. /etc/profile.d/modules.sh\n',
            'date\n',
            'cd '+rundir+'\n',
            'source activate py27 \n',
            'cd '+dir_old+'\n',
            'for f in '+str(year)+'*;do cp -r $f '+dir_new+';done\n',
            'date\n'])
        job_file.close()
 

if adj:
    for year in years:
          name="mask_seviri_adj_"+str(year)
          job_file=open(jobdir+name+'.job','w')
          job_file.writelines([  '#!/bin/sh\n',
                '#PBS -N '+name+'\n',
                '#PBS -e '+name+'.error\n',
                '#PBS -o '+name+'.output\n',
                '#PBS -q shortq\n',                              
                '#PBS -l nodes=1:ppn=1\n',                
		'#PBS -l mem=30gb\n',
                '#PBS -l walltime='+str(time_req)+':00:00\n',
                '. /etc/profile.d/modules.sh\n',
                'date\n',
                'cd '+rundir+'\n',
                'source activate py27 \n',
                'python mask_data_adj.py '+str(nt_back)+' '+str(nt_for)+' '+str(year)+'\n',
                'date\n'])
          job_file.close()

if hist_mask:
    for year in years:
          name="mask_seviri_hist_"+str(year)
          job_file=open(jobdir+name+'.job','w')
          job_file.writelines([  '#!/bin/sh\n',
                '#PBS -N '+name+'\n',
                '#PBS -e '+name+'.error\n',
                '#PBS -o '+name+'.output\n',
                '#PBS -q r630-shortq\n',                              
                '#PBS -l nodes=1:ppn=1\n',                
		        '#PBS -l mem=90gb\n',
                '#PBS -l walltime='+str(time_req)+':00:00\n',
                '. /etc/profile.d/modules.sh\n',
                'date\n',
                'cd '+rundir+'\n',
                'source activate py27 \n',
                'python mask_hist.py '+str(year)+'\n',
                'date\n'])
          job_file.close()



if hist_stats:    
    for time in times:
          name="calc_seviri_hist_"+time
          job_file=open(jobdir+name+'.job','w')
          job_file.writelines([  '#!/bin/sh\n',
                '#PBS -N '+name+'\n',
                '#PBS -e '+name+'.error\n',
                '#PBS -o '+name+'.output\n',
                '#PBS -q r630-shortq\n',                              
                '#PBS -l nodes=1:ppn=1\n',
                '#PBS -l mem=45gb\n',
                '#PBS -l walltime='+str(time_req)+':00:00\n',
                '. /etc/profile.d/modules.sh\n',
                'date\n',
                'cd '+rundir+'\n',
                'source activate py27 \n',
                'python mask_data.py '+time+'\n',
                'date\n'])
          job_file.close()

if hist_clim:    
    for time in times:
          name="calc_seviri_clim_"+time
          job_file=open(jobdir+name+'.job','w')
          job_file.writelines([  '#!/bin/sh\n',
                '#PBS -N '+name+'\n',
                '#PBS -e '+name+'.error\n',
                '#PBS -o '+name+'.output\n',
                '#PBS -q r630-shortq\n',                              
                '#PBS -l nodes=1:ppn=1\n',
                '#PBS -l mem=45gb\n',
                '#PBS -l walltime='+str(time_req)+':00:00\n',
                '. /etc/profile.d/modules.sh\n',
                'date\n',
                'cd '+rundir+'\n',
                'source activate py27 \n',
                'python calc_hist_clim.py '+time+'\n',
                'date\n'])
          job_file.close()
if hist_anom:    
    for year in years:
          name="calc_seviri_anom_"+str(year)
          job_file=open(jobdir+name+'.job','w')
          job_file.writelines([  '#!/bin/sh\n',
                '#PBS -N '+name+'\n',
                '#PBS -e '+name+'.error\n',
                '#PBS -o '+name+'.output\n',
                '#PBS -q r630-shortq\n',                              
                '#PBS -l nodes=1:ppn=1\n',                
		        '#PBS -l mem=45gb\n',
                '#PBS -l walltime='+str(time_req)+':00:00\n',
                '. /etc/profile.d/modules.sh\n',
                'date\n',
                'cd '+rundir+'\n',
                'source activate py27 \n',
                'python calc_hist_anom.py '+str(year)+'\n',
                'date\n'])
          job_file.close()


