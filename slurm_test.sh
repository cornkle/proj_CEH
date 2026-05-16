#!/bin/bash
#SBATCH --partition=short-serial 
#SBATCH -o %j.out 
#SBATCH -e %j.err
#SBATCH --time=05:00
#SBATCH --job-name=MCStest

ipython JASMIN/JASMIN_MCS_CP4_storm_box_anom_v2.py 17 'hist' 
sleep 5m
