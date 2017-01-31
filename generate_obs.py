#!/usr/bin/env python
################################################################################
# generate_obs.py 
# Ken Dixon   
# February/2013
################################################################################
#
# Pre-requisites:
#       + Fill out user-defined parameters
#       + Make sure you have run WPS (geogrid, metgrid, etc.)
#         and have the results in EXPDIR/CASENAME (usually where control run is)
#       + Set up a directory for your nudged experiment in
#              EXPDIR/CASENAME/ASSIMNAME
#         This is where "OBS_DOMAIN101" will be placed
#
# For a given assimilation period, get lightning obs and control simulation data 
# and create model appropriate obs to be used in WRF FDDA nudging
#    + Makes 100% RH profiles for locations nearest ltng strikes for 
#      assimilation times determined by the namelist.input file
#
################################################################################

import sys, os
import calendar, time
import wwlln

###########################
# User-Defined Parameters #
###########################
WRFDIR      = '/path/to/where/you/installed/wrf'         # "WRFV3" should be in this directory
EXPDIR      = '/path/to/your/simulations'                # The directory in which you run wrf should
                                                         # be in this directory ("exp" = experiment)
                                                         # and should be named CASENAME (below)
LTNGSRCDIR  = '/path/to/your/lightning/data'             # Directory containing a WWLLN ".loc" data
                                                         # e.g. "LTNGSRCDIR/CASENAME/A2012062912.loc"
CASENAME    = 'derecho'                                  # name of your experiment
                                                         # e.g. "EXPDIR/CASENAME/"
ASSIMNAME   = 'nudge_10km'                               # name of your nudging experiment (must be directory)
                                                         # e.g. "EXPDIR/CASENAME/ASSIMNAME"
LTNGPREFIX  = 'A'                                        # prefix of the WWLLN lightning files
LTNGEXT     = '.loc'                                     # extension of the WWLLN lightning files
WINDOW      = 5.0         # minutes (will get +/- window/2 for each time)
FILTERRAD   = 0.5         # filter radius in degrees (0.5)
FILTERWIN   = 60.0        # filter window in minutes (60)
PROFILETYPE = 'rh100'     # resulting vertical profile type (only 'rh100' enabled)
READWRFNC   = WRFDIR+'/WRFTOOLS/READ_WRF/read_wrf_nc'   # location of the "read_wrf_nc" tool
                                                        # that must be available
###########################

print('*****************************************')
print('**                                     **')
print('**       LTNG OBS FILE GENERATOR       **')
print('**                                     **')
print('*****************************************')
print('Project      :  '+CASENAME)
print('Obs. Nudging :  '+ASSIMNAME)
print('Profile Type :  '+PROFILETYPE)

##############
# User Check #
##############
USRCHK = raw_input('\n(0) Are these settings correct? (y=yes)  :  ')
if not USRCHK == 'y':
   print('Ok! Exiting.')
   sys.exit()
print('    SUCCESS!')

########################
# Assemble Directories #
########################
CASEDIR  = '%s/%s' %(EXPDIR,CASENAME)
ASSIMDIR = '%s/%s' %(CASEDIR,ASSIMNAME)
LTNGDIR  = '%s/%s' %(LTNGSRCDIR,CASENAME)

####################
# Change Directory #
####################
os.chdir(ASSIMDIR)

#################
# Sanity Checks #
#################
print('\n(1) Checking user parameters and file locations')
wwlln.genobs_prelim_check(CASEDIR,ASSIMDIR,LTNGDIR,"obsn2",ASSIMNAME,PROFILETYPE)
print('    SUCCESS!')

###########################
# Get asssimilation Times #
###########################
print('\n(2) Getting assimilation times from namelist.input')
nmlTDict = wwlln.get_nml_times(ASSIMDIR)
modelStartStr  = "%4.0f%02.0f%02.0f%02.0f%02.0f%02.0f" % (nmlTDict["start_year"],nmlTDict["start_month"],\
                 nmlTDict["start_day"],nmlTDict["start_hour"],nmlTDict["start_minute"],nmlTDict["start_second"])
modelStartTime = calendar.timegm(time.strptime(modelStartStr,"%Y%m%d%H%M%S"))
assimStartTime = modelStartTime + (60*nmlTDict["fdda_start"])
assimEnderTime = modelStartTime + (60*nmlTDict["fdda_end"])
print('    SUCCESS!')

##########################
# Determine Time Windows #
##########################
print('\n(3) Determining Time Windows')
numWindows = wwlln.get_num_windows(assimStartTime,assimEnderTime,WINDOW)
print('    SUCCESS!')

############################
# Get lightning boundaries #
############################
print('\n(4) Getting domain information from GEOGRID netCDF file')
wrfDomDict    = wwlln.get_wrfdomain(CASEDIR,READWRFNC)
ltngDomDict   = wwlln.get_ltngdomain(ASSIMDIR,wrfDomDict)
print('      WRF Lats     Min:  %7.2f     Max:  %7.2f' %(wrfDomDict['minLat'],wrfDomDict['maxLat']))
print('      WRF Lons    West:  %7.2f    East:  %7.2f' %(wrfDomDict['westLon'],wrfDomDict['eastLon']))
print('     LTNG Lats     Min:  %7.2f     Max:  %7.2f' %(ltngDomDict['minLat'],ltngDomDict['maxLat']))
print('     LTNG Lons    West:  %7.2f    East:  %7.2f' %(ltngDomDict['westLon'],ltngDomDict['eastLon']))
filterDomDict = wwlln.get_filterdomain(FILTERRAD,ltngDomDict)
print('   FILTER Lats     Min:  %7.2f     Max:  %7.2f' %(filterDomDict['minLat'],filterDomDict['maxLat']))
print('   FILTER Lons    West:  %7.2f    East:  %7.2f' %(filterDomDict['westLon'],filterDomDict['eastLon']))
print('    SUCCESS!')

#######################
# Get lightning files #
#######################
print('\n(5) Analyzing available WWLLN files (YYYYMMDD only)')
print('    (a) Browsing directory for WWLLN files')
ltngFileDates  = wwlln.get_wwlln_files(LTNGDIR,LTNGPREFIX,LTNGEXT)
print('    (b) Checking that assimilation times are within available files')
assimLtngDates = wwlln.get_assim_ltng_files(assimStartTime,assimEnderTime,ltngFileDates)
print('    SUCCESS!')

print('\n(6) Harvesting observations') 
flashArray = wwlln.harvest_obsn2(LTNGDIR,LTNGPREFIX,LTNGEXT,assimLtngDates,\
                                ltngDomDict,assimStartTime,assimEnderTime,\
                                WINDOW,numWindows,FILTERRAD,FILTERWIN,filterDomDict)
print('Found '+str(len(flashArray))+' strikes')
print('    SUCCESS!')

############################
# Create FDDA Observations #
############################
print('\n(7) Creating OBS_DOMAIN101')
wwlln.write_OBSDOMAIN(flashArray,PROFILETYPE,WRFDIR)
print('    SUCCESS!')
