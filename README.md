# `wrf_ltng_nudge`
Ken Dixon. 2 June 2014.

## Description
This repository contains the WRF code modifications necessary to nudge water
vapor toward saturation below 200 hPa where WWLLN observations are present, as in [Dixon et al. 2016](http://journals.ametsoc.org/doi/abs/10.1175/JTECH-D-15-0188.1). The repository also contains Python scripts for converting [WWLLN](http://wwlln.net) data files into LittleR (the format for WRF FDDA observations).

## Installation 

1. Make sure you can successfully compile WRF or have done so before.
2. Replace "WRFV3/phys/module_fddaobs_rtfdda.F" with the file in this repository.
3. Replace "WRFV3/share/wrf_fddaobs_in.F" with the file in this repository.
4. Recompile WRF.

## Execution

1. Obtain WWLLN ".loc" files.
2. Use "generate_obs.py" to create the LittleR format observations file required by WRF.
3. Move the resulting `OBS_DOMAIN101` file to the directory where WRF is to be run.
4. Set the namelist.
    + `obs_nudge_opt = 1`   (0=off / 1=on)
    + `max_obs = 150000`    (set this to a very high # so all WWLLN strokes are ingested)
    + `fdda_start = 0`      (# of minutes into simulation to start nudging)   
    + `fdda_end = 180`      (# of minutes into simulation to stop nudging)
    + `obs_nudge_mois = 1`  (moisture nudging on/off... also turn temp/wind off)
    + `obs_coef_mois = 3.33E-3` (5min relaxation time scale)
    + `obs_rinxy = 10`          (10km radius of influence)
    + `obsn_rinsig = 0.005`     (set so only one level is affected
    + `obs_twindo = 0.04166667` (in hours, half of the binwidth that you set in "generate_obs.py")
    + `obs_no_pbl_nudge_q = 0`  (keep this off so it isn't turned off in the PBL)
    + For more information, [WRF-ARW](http://www2.mmm.ucar.edu/wrf/users/docs/user_guide_V3/contents.html) documentation.
5. Run WRF.

## Included Files
+ Modified WRF code.
	+ `$WRF/phys/module_fddaobs_rtfdda.F`
	+ `$WRF/share/wrf_fddaobs_in.F`
+ `report.pdf`: a description of the technique and WRF modifications.
+ `namelist.input`: an example namelist.input used for 2012 derecho experiment.
+ `generate_obs.py`: script to generate WRF obs nudging input file, `OBS_DOMAIN101` (you must edit the top of this file).
+ `wwlln.py`: python functions used by `generate_obs.py` (must be in same dir where `generate_obs.py` is run OR in your python path)


## Caveats
+ The modifications were made to WRF v3.4. If you are using a later version of WRF, you will need to carefully merge differences between your version and v3.4 into the files included in this repository.
+ These modifications are not guaranteed to work when other observations are
assimilated via WRF-FDDA. 