Python code to process historical SEVIRI LST data. Created SRA Oct-Dec 2019 for SWIFT and NFLICS projects

mask_data_functions.py  - functins used by everything else
mask_data_wrap.py - creates job files to run the scripts below. JObs are created in directory "jobs"

Naming conventions for files (options set in scripts):
_adj_XbackYforward - imaged masked for missing points in X previous and Y following images
_sigX_countY - image masked for mean+- X sigma, and where < Y counts in distribution 

The scripts are run as follows. 

1.COPY historic_sorted to historic_mask_adj_3back1forward
  RUN jobs for mask_data_adj.py - Mask data with adjacent-in-time pixels masked in origioanl dat.
  (1 job per year)
  Input       historic_mask_adj_3back1forward (copy of historic_sorted)
  Output      historic_mask_adj_3back1forward (now masked)

2.RUN jobs for mask_data.py - Calculate the mean, standard deviation, and counts in the distribution of consecutive image T - (T-1)
  (1 job per time of day. MUST submit the *0.job jobs FIRST THEN the *5.job files)
  Input       historic_mask_adj_3back1forward (masked from step 1) & 
  Output      historic_stats_adj_3back1forward & historic_stats_adj_3back1forward_plots
  
3.MOVE historic_mask_adj_3back1forward to historic_mask_adj_3back1forward_sig2_count50
  RUN jobs for mask_hist.py- Mask data from historical distributions
  (1 job per year)
  Input       historic_stats_adj_3back1forward & historic_sorted & historic_mask_adj_3back1forward_sig2_count50 
  Output      historic_mask_adj_3back1forward_sig2_count50 (now masked) & historic_mask_adj_3back1forward_sig2_count50_plots

4.RUN jobs for calc_hist_clim.py - calculate the mean, standard deviation and counts for the masked historical data
  (1 job per time of day. MUST submit the *0.job jobs FIRST THEN the *5.job files)
  Input       historic_mask_adj_3back1forward_sig2_count50 (now masked from steps 1 & 3)            
  Output      historic_clim & historic_clim_plots
  (1 job per year)

5.RUN jobs for calc_hist_anom.py - calculate the anomolies between masked data and historical climitology
  Input      historic_clim & historic_mask_adj_3back1forward_sig2_count50
  Output     historic_mask_adj_3back1forward_sig2_count (now with grids of anomoly, average anomoly and associated counts) & historic_clim_anom_plots
  (1 job per year)
