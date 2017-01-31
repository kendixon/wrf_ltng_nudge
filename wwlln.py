####################################################################################################
# wwlln.py
# Ken Dixon
# February/2013
####################################################################################################
#
# This python file contains the functions necessary to process WWLLN lightning data
#  + Includes functions to convert WWLLN obs into moisture-nudging soundings for WRF-ARW with edits
#    to the WRF v3.4 code (Ken Dixon, April/2013). This is done by "generate_obs.py"
#
####################################################################################################
import os, sys, glob
import numpy as np
import calendar, time
from datetime import datetime
from datetime import timedelta
import operator
from math import radians, cos, sin, asin, sqrt

##############################
# generate_obs: Sanity Check #
##############################
def genobs_prelim_check(CASEDIR,ASSIMDIR,LTNGDIR,HARVEST,ASSIMNAME,PROFILETYPE):
   if os.path.exists(CASEDIR) == False:
      print('FATAL ERROR:   '+CASEDIR+' not found! Exiting.\n')
      sys.exit()
   if os.path.exists(ASSIMDIR) == False:
      print('FATAL ERROR:   '+ASSIMDIR+' not found! Exiting.\n')
      sys.exit()
   if os.path.exists(LTNGDIR) == False:
      print('FATAL ERROR:   '+LTNGDIR+' not found! Exiting.\n')
      sys.exit()
   if (HARVEST != 'obsn1') & (HARVEST != 'obsn2'):
      print('FATAL ERROR:   '+HARVEST+' not a valid harvesting method! Exiting.\n')
      sys.exit()
   if (PROFILETYPE != 'rh100'):
      print('FATAL ERROR:   '+PROFILETYPE+' not a valid vertical profile type! Exiting.\n')
      sys.exit()
   # namelist.input for assimilation
   try:
      with open('%s/namelist.input'%(ASSIMDIR)) as f: pass
   except IOError as e:
      print('FATAL ERROR:   '+ASSIMNAME+'/namelist.input not found! Exiting.\n')
      sys.exit()
   # GEOGRID file for domain data
   try:
      with open(CASEDIR+'/geo_em.d01.nc') as gg: pass
   except IOError as e:
      print('FATAL ERROR:   geo_em.d01.nc not found! Exiting.\n')
      sys.exit()
   # namelist.wps for domain data
   try:
      with open(CASEDIR+'/namelist.wps') as wps: pass
   except IOError as e:
      print('FATAL ERROR:   namelist.wps not found! Exiting.\n')
      sys.exit()

########################
# energy: Sanity Check #
########################
def energy_prelim_check(LTNGDIR,USEWRF,CASEDIR,ASSIMDIR,LTNGLAT,LTNGLON,STARTTIMESTR,ENDTIMESTR):
   if os.path.exists(LTNGDIR) == False:
      print('FATAL ERROR:   '+LTNGDIR+' not found! Exiting.\n')
      sys.exit()
   if os.path.exists('%s/images'%(LTNGDIR)) == False:
      print('ERROR      :   '+LTNGDIR+'/images not found! Exiting.\n')
      sys.exit()
   if USEWRF:
      if os.path.exists(CASEDIR) == False:
         print('FATAL ERROR:   '+CASEDIR+' not found! Exiting.\n')
         sys.exit()
      if os.path.exists(ASSIMDIR) == False:
         print('FATAL ERROR:   '+ASSIMDIR+' not found! Exiting.\n')
         sys.exit()
      # namelist.input for assimilation
      try:
         with open('%s/namelist.input'%(ASSIMDIR)) as f: pass
      except IOError as e:
         print('FATAL ERROR:   '+ASSIMNAME+'/namelist.input not found! Exiting.\n')
         sys.exit()
      # GEOGRID file for domain data
      try:
         with open(CASEDIR+'/geo_em.d01.nc') as gg: pass
      except IOError as e:
         print('FATAL ERROR:   geo_em.d01.nc not found! Exiting.\n')
         sys.exit()
      # namelist.wps for domain data
      try:
         with open(CASEDIR+'/namelist.wps') as wps: pass
      except IOError as e:
         print('FATAL ERROR:   namelist.wps not found! Exiting.\n')
         sys.exit()
   else:
      STARTTIME = calendar.timegm(time.strptime(STARTTIMESTR,"%Y%m%d%H%M"))
      ENDTIME   = calendar.timegm(time.strptime(ENDTIMESTR,"%Y%m%d%H%M"))
      if LTNGLAT[1]<=LTNGLAT[0]:
         print('FATAL ERROR:   Latitudes do not make sense')
         sys.exit()
      if STARTTIME >= ENDTIME:
         print('FATAL ERROR:   Times do not make sense')
         sys.exit()

#######################################
# Get Times from WRF Namelist(.input) #
#######################################
def get_nml_times(ASSIMDIR):
   assimTimeLines = ["fdda_start","fdda_end","start_year","start_month",\
                     "start_day","start_hour","start_minute","start_second",\
                     "end_year","end_month","end_day","end_hour","end_minute",\
                     "end_second"]
   nmlTDict  = {}   # Namelist Time Dictionary for above parameters
   f = open('%s/namelist.input' %(ASSIMDIR)) 
   for line in f:
      if line.split('=')[0].strip() in assimTimeLines:
         nmlTDict[line.split('=')[0].strip()] = float(line.split('=')[1].split(',')[0].strip())
   f.close()
   return nmlTDict

##############################
# Get number of time windows #
##############################
def get_num_windows(STARTTIME,ENDTIME,WINDOW):
   totalTimeLength = ENDTIME-STARTTIME
   if totalTimeLength%(60.0*WINDOW) == 0.0:
      numWindows = int(totalTimeLength/(60.0*WINDOW))
   else:
      print('FATAL ERROR:   Time Length not a multiple of window/bin length! Exiting.')
      sys.exit()
   return numWindows

####################################################################
# map_ltng: Get Reference Lats/Lons, Domain Size from namelist.wps #
####################################################################
def get_basemap_geogrid(CASEDIR):
   basemapDict = {}
   wpsfile = open(CASEDIR+'/namelist.wps','r')
   for line in wpsfile:
      if 'ref_lat' in line:
         basemapDict['lat_0'] = float(line.split('=')[1].split(',')[0].strip())
      if 'ref_lon' in line:
         basemapDict['lon_0'] = float(line.split('=')[1].split(',')[0].strip())
      if 'truelat1' in line:
         basemapDict['lat_1'] = float(line.split('=')[1].split(',')[0].strip())
      if 'truelat2' in line:
         basemapDict['lat_2'] = float(line.split('=')[1].split(',')[0].strip())
      if 'dx' in line:
         basemapDict['dx'] = float(line.split('=')[1].split(',')[0].strip())
      if 'dy' in line:
         basemapDict['dy'] = float(line.split('=')[1].split(',')[0].strip())
      if 'e_we' in line:
         basemapDict['gridx'] = float(line.split('=')[1].split(',')[0].strip())
      if 'e_sn' in line:
         basemapDict['gridy'] = float(line.split('=')[1].split(',')[0].strip())
   wpsfile.close()
   return basemapDict

#########################
# Get WRF Domain Bounds #
#########################
def get_wrfdomain(CASEDIR,READWRFNC):
   # set up
   wrfDomDict = {}
   # link in the GEOGRID file
   os.system('ln -sf '+CASEDIR+'/geo_em.d01.nc .')
   # Use far sides of Arakawa-C grids to include most ltng
   os.system(READWRFNC+' geo_em.d01.nc -v XLAT_V >& xlatv.out')
   os.system(READWRFNC+' geo_em.d01.nc -w XLONG_U >& xlongu.out')
   # GEOGRID Latitude Bounds
   latfile = open('xlatv.out','r')
   for line in latfile:
      if (line[1:20] == '0000-00-00_00:00:00'):
         newLine = line.replace('_',' ').replace(':',' ').replace('=',' ')
         lineArr = newLine.split()
         wrfDomDict['minLat'] = float(lineArr[5])
         wrfDomDict['maxLat'] = float(lineArr[7])
   latfile.close()
   os.system('rm xlatv.out')
   # GEOGRID Longitude Bounds
   lonfile = open('XLONG_U.out','r')
   firstLonLine = False
   longitudes = []
   for line in lonfile:
      if firstLonLine == False:
         if '0000-00-00_00:00:00' in line:
            firstLonLine = True
      else:
         lineArr = line.split()
         for elem in lineArr:
            longitudes.append(float(elem))
   lonfile.close()
   os.system('rm XLONG_U.out')
   os.system('rm xlongu.out')
   if (all((elem>=0) for elem in longitudes)) | (all((elem<=0) for elem in longitudes)):
      # Domain does not cross Prime Meridian or Int'l Date Line
      wrfDomDict['westLon'] = float(min(longitudes))
      wrfDomDict['eastLon'] = float(max(longitudes))
   else:
      # Domain does cross Prime Meridian or Int'l Date Line... hopefully not both
      # CENT_LON from namelist.wps
      wpsfile = open(CASEDIR+'/namelist.wps','r')
      for line in wpsfile:
         if 'ref_lon' in line:
           centLon = float(line.split('=')[1].split(',')[0].strip())
      wpsfile.close()
      # Check if crossing Prime Meridian or Int'l Date Line
      if (centLon < 90) & (centLon > -90):   # Crossing Prime Meridian
         wrfDomDict['westLon'] = float(min(longitudes))
         wrfDomDict['eastLon'] = float(max(longitudes))
      else:                                  # Crossing Date Line
         wrfDomDict['westLon'] = float(min([r for r in longitudes if r > 0])) # smallest positive long
         wrfDomDict['eastLon'] = float(max([r for r in longitudes if r < 0])) # largest negative long
#   lonfile = open('xlongu.out','r')
#   for line in lonfile:
#      if (line[1:20] == '0000-00-00_00:00:00'):
#         newLine = line.replace('_',' ').replace(':',' ').replace('=',' ')
#         lineArr = newLine.split()
#         lon1 = float(lineArr[5])
#         lon2 = float(lineArr[7])
#   lonfile.close()
#   os.system('rm XLONG_U.out')
#   os.system('rm xlongu.out')
#   # CENT_LON from namelist.wps
#   wpsfile = open(CASEDIR+'/namelist.wps','r')
#   for line in wpsfile:
#      if 'ref_lon' in line:
#        centLon = float(line.split('=')[1].split(',')[0].strip())
#   wpsfile.close()
#   # Determine eastLon / westLon (in case of Dateline straddling)
#   if (centLon < lon1) & (centLon < lon2):
#      wrfDomDict['eastLon'] = min([lon1,lon2])
#      wrfDomDict['westLon'] = max([lon1,lon2])
#   else:
#      wrfDomDict['eastLon'] = max([lon1,lon2])
#      wrfDomDict['westLon'] = min([lon1,lon2])
   # remove geogrid link in
   os.system('rm geo_em.d01.nc')
   return wrfDomDict

#################################################
# WRF Domain Bounds --> Lightning Domain Bounds #
#################################################
def get_ltngdomain(ASSIMDIR,wrfDomDict):
   # set up
   ltngDomDict = {}
   # RINXY from namelist.input
   f = open('%s/namelist.input' %(ASSIMDIR)) 
   for line in f:
      if 'obs_rinxy' in line:
         rinxy = float(line.split('=')[1].split(',')[0].strip())
   rinxyDeg = rinxy/112.0
   f.close()
   ### Latitude
   #### min latitude
   if wrfDomDict['minLat']-rinxyDeg <= -90:
      ltngDomDict['minLat'] = -90.0
   else:
      ltngDomDict['minLat'] = wrfDomDict['minLat']-rinxyDeg
   #### max latitude
   if wrfDomDict['maxLat']+rinxyDeg >= 90:
      ltngDomDict['maxLat'] = 90.0
   else:
      ltngDomDict['maxLat'] = wrfDomDict['maxLat']+rinxyDeg
   ### Longitude
   ltngDomDict['westLon'] = wrfDomDict['westLon']-rinxyDeg
   ltngDomDict['eastLon'] = wrfDomDict['eastLon']+rinxyDeg
   if ltngDomDict['westLon'] < -180:
      ltngDomDict['westLon'] = ltngDomDict['westLon']+360
   if ltngDomDict['eastLon'] > 180:
      ltngDomDict['eastLon'] = ltngDomDict['eastLon']-360
   return ltngDomDict

#############################################
# Get Domain Bounds for Filtering Lightning #
#############################################
def get_filterdomain(FILTERRAD,ltngDomDict):
   # Filtering is done using 0.5x0.5 degree boxes. If a lightning strike 
   # has a "buddy" within the last 60 minutes in its box, or any adjacent 
   # box, then it passes the filter. As specified in "harvest_obsn2", this 
   # filtering is not done for lone strikes detected by 6+ stations
   # + Filtering box size can be changed by user FILTERRAD = 0.5
   
   # Inverse of filter radius
   INVFR = 1.0/FILTERRAD
   # Set up Dictionary
   filterDomDict = {}
   ### Latitude
   #### min latitude
   filterDomDict['minLat'] = (np.floor(ltngDomDict['minLat']*INVFR)/INVFR)-FILTERRAD
   if filterDomDict['minLat'] <=-90:
      filterDomDict['minLat'] = -90.0
   #### max latitude
   filterDomDict['maxLat'] = (np.ceil(ltngDomDict['maxLat']*INVFR)/INVFR)+FILTERRAD
   if filterDomDict['maxLat'] >= 90:
      filterDomDict['maxLat'] = 90.0
   ### Longitude
   filterDomDict['westLon'] = (np.floor(ltngDomDict['westLon']*INVFR)/INVFR)-FILTERRAD
   filterDomDict['eastLon'] = (np.ceil(ltngDomDict['eastLon']*INVFR)/INVFR)+FILTERRAD
   if filterDomDict['westLon'] < -180:
      filterDomDict['westLon'] = filterDomDict['westLon']+360
   if filterDomDict['eastLon'] > 180:
      filterDomDict['eastLon'] = filterDomDict['eastLon']-360
   return filterDomDict

###################################
# Find WWLLN files in a directory #
###################################
def get_wwlln_files(LTNGDIR,LTNGPREFIX,LTNGEXT):
   ltngFileDates = []
   for root, dirnames, filenames in os.walk(LTNGDIR):
      for filename in glob.glob(root+'/*'+LTNGEXT):
         ltngFileDates.append(filename.split(LTNGDIR+'/')[1].split(LTNGPREFIX)[1].split(LTNGEXT)[0])
   if len(ltngFileDates) == 0:
      print('FATAL ERROR:   No Lightning files in directory! Exiting.\n')
      sys.exit()
   ltngDateTruth = [len(elem)==len(ltngFileDates[0]) for elem in ltngFileDates]
   if sum(ltngDateTruth) != len(ltngDateTruth):
      print('FATAL ERROR:   Lightning files of different datestring format! Exiting.\n')
      sys.exit()
   if (len(ltngFileDates[0]) != 8):
      print('FATAL ERROR:   Lightning files are not in YYYYMMDD format! Exiting.\n')
      sys.exit()
   return ltngFileDates

def get_wwlln_files_rt(LTNGDIR,LTNGPREFIX,LTNGEXT):
   ltngFileDates = []
   for root, dirnames, filenames in os.walk(LTNGDIR):
      for filename in glob.glob(root+'/*'+LTNGEXT):
         ltngFileDates.append(filename.split(LTNGDIR+'/')[1].split(LTNGPREFIX)[1].split(LTNGEXT)[0])
   if len(ltngFileDates) == 0:
      print('FATAL ERROR:   No Lightning files in directory! Exiting.\n')
      sys.exit()
   ltngDateTruth = [len(elem)==len(ltngFileDates[0]) for elem in ltngFileDates]
   if sum(ltngDateTruth) != len(ltngDateTruth):
      print('FATAL ERROR:   Lightning files of different datestring format! Exiting.\n')
      sys.exit()
   if (len(ltngFileDates[0]) != 8):
      print('FATAL ERROR:   Lightning files are not in YYYYMMDD format! Exiting.\n')
      sys.exit()
   return ltngFileDates

##################################
# Check if strike is in a domain #
##################################
def strike_in_domain(strikeLat,strikeLon,ltngDomDict):
   if (strikeLat <= ltngDomDict['maxLat']) & (strikeLat >= ltngDomDict['minLat']):
      if ltngDomDict['eastLon'] < ltngDomDict['westLon']:                  # dateline issue!
         if (strikeLon >= ltngDomDict['westLon']) & (strikeLon <= 180):    # inside west box
            ltngTruth = True
         elif (strikeLon <= ltngDomDict['eastLon']) & (strikeLon >= -180): # inside east box
            ltngTruth = True
         else:
            ltngTruth = False
      else:                                               # no dateline issue!
         if (strikeLon <= ltngDomDict['eastLon']) & (strikeLon >= ltngDomDict['westLon']):
            ltngTruth = True
         else:
            ltngTruth = False
   else:
      ltngTruth = False
   return ltngTruth

#################################################
# Get Ltng Files that cover assimilation period #
#################################################
def get_assim_ltng_files(assimStartTime,assimEnderTime,ltngFileDates):
   assimStartDayStr  = time.strftime("%Y%m%d",time.gmtime(assimStartTime))
   assimStartDayTime = calendar.timegm(time.strptime(assimStartDayStr,"%Y%m%d"))
   assimEnderDayStr  = time.strftime("%Y%m%d",time.gmtime(assimEnderTime))
   assimEnderDayTime = calendar.timegm(time.strptime(assimEnderDayStr,"%Y%m%d"))
   currentDayTime    = assimStartDayTime
   assimLtngDates    = []
   while currentDayTime <= assimEnderDayTime:
      currentDayStr = time.strftime("%Y%m%d",time.gmtime(currentDayTime))
      if currentDayStr not in ltngFileDates:
         print('FATAL ERROR:   Missing '+currentDayStr+' lightning file! Exiting.\n')
         sys.exit()
      else:
         assimLtngDates.append(currentDayStr)
      currentDayTime = currentDayTime + 86400
   return assimLtngDates

###############################################
# Haversine Calculation of Distance on Sphere #
###############################################
def haversine(lat1,lon1,lat2,lon2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    km = 6367 * c
    return km


#####################################
# Harvest 'obsn1' Lightning Strikes #
#  + Each strike is an 100% RH ob   #
#    at the actual time of stroke   #
#####################################
def harvest_obsn1(LTNGDIR,LTNGPREFIX,LTNGEXT,assimLtngDates,ltngDomDict,assimStartTime,assimEnderTime):
   flashArray     = []
   for day in assimLtngDates:
      f = open(LTNGDIR+'/'+LTNGPREFIX+day+LTNGEXT,'r')
      for line in f:
         newLine = line.replace('/',' ').replace(':',' ').replace(',',' ')
         lineArr = np.fromstring(newLine,sep=' ')
         sid = strike_in_domain(lineArr[6],lineArr[7],ltngDomDict)
         if sid:
            lineTimeStr = "%04.0f%02.0f%02.0f%02.0f%02.0f%02.0f" %(lineArr[0],\
                          lineArr[1],lineArr[2],lineArr[3],lineArr[4],np.floor(lineArr[5]))
            lineTime    = calendar.timegm(time.strptime(lineTimeStr,"%Y%m%d%H%M%S"))
            if (lineTime >= assimStartTime) & (lineTime < assimEnderTime):
               flashArray.append([lineArr[6],lineArr[7],lineTimeStr])
      f.close()
   flashArray.sort(key=lambda x:calendar.timegm(time.strptime(x[2],"%Y%m%d%H%M%S")))
   return flashArray

#####################################
# Harvest 'obsn2' Lightning Strikes #
#  + Each strike is a 100% RH ob    #
#    and time is center of window   #
#####################################
def harvest_obsn2(LTNGDIR,LTNGPREFIX,LTNGEXT,assimLtngDates,\
                  ltngDomDict,assimStartTime,assimEnderTime,\
                  WINDOW,numWindows,FILTERRAD,FILTERWIN,filterDomDict):

   print('    (a) Obtaining Raw Lightning Strikes')
   checkFlashArray = []
   for day in assimLtngDates:
      f = open(LTNGDIR+'/'+LTNGPREFIX+day+LTNGEXT,'r')
      for line in f:
         newLine = line.replace('/',' ').replace(':',' ').replace(',',' ')
         lineArr = np.fromstring(newLine,sep=' ')
         # Check if strike is in domain for checking flashes for filtering
         sid = strike_in_domain(lineArr[6],lineArr[7],filterDomDict)
         if sid:
            lineTimeStr = "%04.0f%02.0f%02.0f%02.0f%02.0f%02.0f" %(lineArr[0],\
                          lineArr[1],lineArr[2],lineArr[3],lineArr[4],np.floor(lineArr[5]))
            lineTime    = calendar.timegm(time.strptime(lineTimeStr,"%Y%m%d%H%M%S"))
            if (lineTime >= assimStartTime-(FILTERWIN*60.0)) &\
               (lineTime < assimEnderTime):
               checkFlashArray.append([lineArr[6],lineArr[7],lineTimeStr,lineTime,lineArr[9]])
      f.close()
   checkFlashArray.sort(key=lambda x:calendar.timegm(time.strptime(x[2],"%Y%m%d%H%M%S")))
   CFA = checkFlashArray

   print('    (b) Filtering Lightning Strikes') 
   INVFR = 1.0/FILTERRAD     # Inverse of filter radius
   flashArray = []           # This array will hold filtered flashes, in domain
   for j in range(len(CFA)):
      if (CFA[j][3] >= assimStartTime) & (CFA[j][3] < assimEnderTime):   # If in assim times...
         sid = strike_in_domain(CFA[j][0],CFA[j][1],ltngDomDict)         #       and domain...
         if sid:
            if (CFA[j][4] >= 6):                            # 6+ stations, add to flashArray
               windowIndex       = np.floor((CFA[j][3]-assimStartTime)/300.0)
               binnedLineTime    = assimStartTime+((WINDOW/2.0)*60.0)+\
                                   (windowIndex*WINDOW*60.0)
               binnedLineTimeStr = time.strftime("%Y%m%d%H%M%S",time.gmtime(binnedLineTime))
               flashArray.append([CFA[j][0],CFA[j][1],binnedLineTimeStr])
            else:                                           # 5 stations; buddy check before adding
               # SET Max/Min latitude/longitudes for filtering lightning
               tempDomDict = {}
               ### Latitude
               ##### min latitude
               tempDomDict['minLat'] = (np.floor(CFA[j][0]*INVFR)/INVFR)-FILTERRAD
               if tempDomDict['minLat'] <=-90:
                  tempDomDict['minLat'] = -90.0
               #### max latitude
               tempDomDict['maxLat'] = (np.ceil(CFA[j][0]*INVFR)/INVFR)+(2.0*FILTERRAD)
               if tempDomDict['maxLat'] >= 90:
                  tempDomDict['maxLat'] = 90.0
               ### Longitude
               tempDomDict['westLon'] = (np.floor(CFA[j][1]*INVFR)/INVFR)-FILTERRAD
               tempDomDict['eastLon'] = (np.ceil(CFA[j][1]*INVFR)/INVFR)+(2.0*FILTERRAD)
               if tempDomDict['westLon'] < -180:
                  tempDomDict['westLon'] = tempDomDict['westLon']+360
               if tempDomDict['eastLon'] > 180:
                  tempDomDict['eastLon'] = tempDomDict['eastLon']-360
               for k in range(j):
                  sidfilt = strike_in_domain(CFA[k][0],CFA[k][1],tempDomDict)
                  if (CFA[k][3] >= (CFA[j][3]-(60.0*FILTERWIN))) & \
                     (CFA[k][3] <= CFA[j][3] ) & sidfilt: # If strike in last 60m, and in +/-0.5deg
                     windowIndex       = np.floor((CFA[j][3]-assimStartTime)/300.0)
                     binnedLineTime    = assimStartTime+((WINDOW/2.0)*60.0)+\
                                         (windowIndex*WINDOW*60.0)
                     binnedLineTimeStr = time.strftime("%Y%m%d%H%M%S",time.gmtime(binnedLineTime))
                     flashArray.append([CFA[j][0],CFA[j][1],binnedLineTimeStr])
                     break  # No need to find more than one strike to verify
   return flashArray

######################################
# Harvest Realtime Lightning Strikes #
#  + Each strike counts!             #
######################################
def harvest_rt(ltngDomDict,centlat,centlon,radius,startTime,endTime):
   """
   Harvest realtime lightning given a starting and ending date.
   If radius is 0, then centlat/centlon are ignored and only ltngDomDict 
   is used to filter lightning strikes.

   """

   LTNGDIR    = "/path/to/a/ltng/directory"
   LTNGPREFIX = "A"
   LTNGEXT    = ".loc"
   ftimeArray = []    # File times used
   flashArray = []    # Flashes
   currTime   = startTime
   currFileTime = '%4d%02d%02d%02d%s0' %(currTime.year,currTime.month,\
                  currTime.day,currTime.hour,str(currTime.minute)[0])
   currFileDateTime = datetime(currTime.year,currTime.month,\
                     currTime.day,currTime.hour,int(str(currTime.minute)[0]+'0'))
   while (currTime < endTime):
      try:
         f = open(LTNGDIR+'/'+LTNGPREFIX+currFileTime+LTNGEXT,'r')
      except IOError:
         print('      Cannot open '+LTNGDIR+'/'+LTNGPREFIX+currFileTime+LTNGEXT)
         currTime     = currTime+timedelta(0,600)  # 10min realtime files
         currFileTime = '%4d%02d%02d%02d%s0' %(currTime.year,currTime.month,\
                        currTime.day,currTime.hour,str(currTime.minute)[0])
         currFileDateTime = datetime(currTime.year,currTime.month,\
                           currTime.day,currTime.hour,int(str(currTime.minute)[0]+'0'))
         continue
      ftimeArray.append(currFileDateTime)
      for line in f:
         newLine  = line.replace('/',' ').replace(':',' ').replace(',',' ')
         lineArr  = np.fromstring(newLine,sep=' ')
         sid = strike_in_domain(lineArr[6],lineArr[7],ltngDomDict)
         if sid:                 # Strike in loose domain (lat/lon box)
            if radius == 0.0:
               lineTime = datetime(int(lineArr[0]),int(lineArr[1]),int(lineArr[2]),\
                          int(lineArr[3]),int(lineArr[4]),int(np.floor(lineArr[5])))
               if (lineTime >= startTime) & (lineTime < endTime):
                  flashArray.append({'lat':lineArr[6],'lon':lineArr[7],'time':lineTime})
            else:
               dist = haversine(centlat,centlon,lineArr[6],lineArr[7])
               if dist <= radius:   # Strike in radius (km)
                  lineTime = datetime(int(lineArr[0]),int(lineArr[1]),int(lineArr[2]),\
                             int(lineArr[3]),int(lineArr[4]),int(np.floor(lineArr[5])))
                  if (lineTime >= startTime) & (lineTime < endTime):
                     flashArray.append({'lat':lineArr[6],'lon':lineArr[7],'time':lineTime})
      f.close()
      currTime     = currTime+timedelta(0,600)  # 10min realtime files
      currFileTime = '%4d%02d%02d%02d%s0' %(currTime.year,currTime.month,\
                     currTime.day,currTime.hour,str(currTime.minute)[0])
      currFileDateTime = datetime(currTime.year,currTime.month,\
                        currTime.day,currTime.hour,int(str(currTime.minute)[0]+'0'))
   returnDict = {'flashes':flashArray,'ftimes':ftimeArray}
   return returnDict


################################################
# Convert Lightning Flash Array to OBS_DOMAIN  #
################################################
def write_OBSDOMAIN(flashArray,PROFILETYPE,WRFDIR):
   wrfObsFile = open('OBS_DOMAIN101','w')
   for k in range(len(flashArray)):
      # Temporary Lat/Lon/Time
      lat  = flashArray[k][0]
      lon  = flashArray[k][1]
      time = flashArray[k][2]
      if PROFILETYPE =='rh100':
         # Determine Levels (only 2 to save space, will nudge whole column)
         plevels = list(np.arange(1000,975,-12.5)) # WRF won't take data above 80mb
         # Write Sounding Header
         wrfObsFile.write(' %14s\n' %(time))
         wrfObsFile.write('%9.2f %9.2f\n' %(lat,lon))
         wrfObsFile.write('%45s%40s\n' %('WWLLN','WWLLN'))
         # Platform, elevation, is_sounding, bogus, # of levels
         wrfObsFile.write('  %5s%38d.%6s%6s %6d\n' %('WWLLN',0,'T','T',len(plevels)))  
         # Write Sounding Data
         for m in range(len(plevels)):
            wrfObsFile.write('%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f\n' \
            %(100*plevels[m],0,-888888,-888888,-888888,-888888,-888888,-888888,-888888,-888888,100,0))
   wrfObsFile.close()

####################################################
# Form a simple Box Domain around a center lat/lon #
####################################################
def box_domain(centlat,centlon,halfwidth):
   """
   Form a domain around centlat/centlon of specified halfwidth
   Ex: centlat/lon = 40,30  ; halfwidth = 0.5
       lats = (39.5,40.5)  lons = (39.5,30.5)
   """
   # set up
   ltngDomDict = {}
   hwDeg = halfwidth  # degrees
   ### Latitude
   #### min latitude
   if centlat-hwDeg <= -90:
      ltngDomDict['minLat'] = -90.0
   else:
      ltngDomDict['minLat'] = centlat-hwDeg
   #### max latitude
   if centlat+hwDeg >= 90:
      ltngDomDict['maxLat'] = 90.0
   else:
      ltngDomDict['maxLat'] = centlat+hwDeg
   ### Longitude
   ltngDomDict['westLon'] = centlon-hwDeg
   ltngDomDict['eastLon'] = centlon+hwDeg
   if ltngDomDict['westLon'] < -180:
      ltngDomDict['westLon'] = ltngDomDict['westLon']+360
   if ltngDomDict['eastLon'] > 180:
      ltngDomDict['eastLon'] = ltngDomDict['eastLon']-360
   return ltngDomDict
