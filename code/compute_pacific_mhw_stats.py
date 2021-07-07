# JBK 2020-06-08
# process marine heatwave statistics for Pacific Islands region
# 1. read SST data from observations and CMIP6 (historical and future)
# 2. trim to required region
# 3. compute MHW statistics for region
# 4. store required data to netCDF
# OUTPUT: mhw_cats/mhw_cats.pac_is.<model>.<scenario>.<variant>.nc

# load required modules
import numpy as np
import xarray as xr
import pandas
import glob
from datetime import date
import time
import cftime
import scipy.io as io

# marine heatwave modules
import marineHeatWaves_jbk as mhw
import ecoliver as eco

# set paths and filenames
inpath = ''  # location of model source list
outpath = '' # output directory

# define which scenarios to use
scen_h = 'historical'
scen_f = 'ssp585'        # ssp126 or ssp585

# define bounds of stored region
reg_lab = 'pac_is'
reg_bnds = [120, 220, -40, 15]   # 120E-140W, 40S-15N

# define climatological baseline period
clim_b = [1995,2014]

# read model lists
# model list .txt files should have OBS path as first entry
modfile_h = open(inpath + 'cmip6_gadi_' + scen_h + '.txt', "r")
modfile_f = open(inpath + 'cmip6_gadi_' + scen_f + '.txt', "r")
modlist_h = modfile_h.read().splitlines()
modlist_f = modfile_f.read().splitlines()
modfile_h.close()
modfile_f.close()

nmods = len(modlist_h)

# range of models to process in this run
m_which = range(0,nmods)
print('PROCESS LIST')
for m in m_which:
 print('[' + str(m) + ']: ' + modlist_h[m])

# loop over all models in list
for m in m_which:
 
 # split path names to retrieve model parameters
 modcode_h = modlist_h[m].split("/")
 modcode_f = modlist_f[m].split("/")
 
 if m == 0:
  # set OBS label, e.g. NOAA_OISST.AVHRR.v2-1_modified
  modcode = '.'.join(modcode_h[4:7])
 else:
  # set model label, e.g. <source_id>.historical+ssp126.<variant_label>
  modcode = modcode_h[8] + '.' + scen_h + '+' + scen_f + '.' + modcode_h[10]
  # check the model names and experiment IDs are the same, only for CMIP6 (skip OBS)
  if not(modcode_h[8] == modcode_f[8]) or not(modcode_h[10] == modcode_f[10]):
   print('--> ERROR: problem in file paths!')
 
 print('Processing... ' + modcode)
 
 # read netcdf filenames from historical and scenario paths
 filenames_h = sorted(glob.glob(modlist_h[m] + '*.nc'))
 filenames_f = sorted(glob.glob(modlist_f[m] + '*.nc'))
 infiles = filenames_h + filenames_f
 infiles = list(dict.fromkeys(infiles))  # remove duplicates, for OBS
 
 # create output filename
 outfile = outpath + 'mhw_cats.' + reg_lab + '.' + modcode

 # ----------------------------------------------------------------------
 # Handle some CMIP6 exceptions
 #-----------------------------------------------------------------------
 # EXCEPTION: for IPSL-CM6A-LR, ssp126, there is a problem in the time array in
 # tos_Oday_IPSL-CM6A-LR_ssp126_r1i1p1f1_gn_21010101-23001231.nc, v20190903
 if modcode_f[8]=='IPSL-CM6A-LR' and modcode_f[9]=='ssp126':
  if modcode_f[10]=='r1i1p1f1' and modcode_f[14]=='v20190903':
   infiles=infiles[0:2]
 
 # EXCEPTION: for IPSL-CM6A-LR, ssp585, there is a problem in the time array in
 # tos_Oday_IPSL-CM6A-LR_ssp585_r1i1p1f1_gn_21010101-23001231.nc, v20190903
 if modcode_f[8]=='IPSL-CM6A-LR' and modcode_f[9]=='ssp585':
  if modcode_f[10]=='r1i1p1f1' and modcode_f[14]=='v20190903':
   infiles=infiles[0:2]
 
 # EXCEPTION: for MRI-ESM2-0, ssp126, there is a problem in the time array in
 # tos_Oday_MRI-ESM2-0_ssp126_r1i1p1f1_gn_22510101-23001231.nc, v20200222
 if modcode_f[8]=='MRI-ESM2-0' and modcode_f[9]=='ssp126':
  if modcode_f[10]=='r1i1p1f1' and modcode_f[14]=='v20200222':
   infiles=infiles[0:5]
 
 # EXCEPTION: for MRI-ESM2-0, ssp585, the file following file is missing
 # tos_Oday_MRI-ESM2-0_ssp585_r1i1p1f1_gn_21010101-21501231.nc, v20200120
 if modcode_f[8]=='MRI-ESM2-0' and modcode_f[9]=='ssp585':
  if modcode_f[10]=='r1i1p1f1' and modcode_f[14]=='v20200120':
   infiles=infiles[0:5]
   
 # ----------------------------------------------------------------------
 # End of CMIP6 exceptions
 #-----------------------------------------------------------------------
 
 # load data
 print('Loading data... ')
 start = time.time()
 ds = xr.open_mfdataset(infiles,combine='nested',concat_dim='time')
 end = time.time()
 print(end - start)
 
 # remove duplicate times if historial and scenario runs have overlapping elements
 # this is the case for CESM2/historical/r4i1p1f1/Oday/tos/gn/v20190308,
 # CESM2/ssp126/r4i1p1f1/Oday/tos/gn/v20200528/,
 # and CESM2/ssp585/r4i1p1f1/Oday/tos/gn/v20200528/
 _, index = np.unique(ds['time'], return_index=True)   # set of unique times
 ds = ds.isel(time=index)   # keep only unique times
 
 # rename SST variable in OBS dataset to CMIP6 standard: tos
 if m == 0:
  ds = ds.rename_vars({'sst':'tos'})
 
 # rename coords if not lat/lon
 if hasattr(ds,'longitude'):
  if hasattr(ds,'lon'):
   print('Warning: this model has both lon and longitude arrays..')
   print('Deleting longitude and latitude...')
   ds = ds.drop_vars({'latitude','longitude'})
  else:
   ds = ds.rename_vars({'longitude':'lon','latitude':'lat'})
 elif hasattr(ds,'nav_lon'):
  ds = ds.rename_vars({'nav_lon':'lon','nav_lat':'lat'})
 
 # store min and max lon, may be required for rotated lon coords
 min_lon = np.min(ds.lon.values)
 max_lon = np.max(ds.lon.values)
 
 print(min_lon)
 print(max_lon)
 
 # trim to required region
 # first case: lon is [0:360]
 if (max_lon > reg_bnds[1]):
  da = ds.where((ds.lon > reg_bnds[0]) & (ds.lon < reg_bnds[1])
    & (ds.lat > reg_bnds[2]) & (ds.lat < reg_bnds[3]), drop=True)
 # second case: lon is [-180:180]
 elif (max_lon > reg_bnds[0]) & (max_lon < reg_bnds[1]):
  da = ds.where(((ds.lon > reg_bnds[0]) | (ds.lon < reg_bnds[1]-360))
    & (ds.lat > reg_bnds[2]) & (ds.lat < reg_bnds[3]), drop=True)
 # third case: lon is [-300:60]
 elif (max_lon < reg_bnds[0]):
  da = ds.where((ds.lon > reg_bnds[0]-360) & (ds.lon < reg_bnds[1]-360)
    & (ds.lat > reg_bnds[2]) & (ds.lat < reg_bnds[3]), drop=True)
 
 # trim time, and store start and end of time array as labels
 # start at 1 Jan 1982
 # end at 31 Dec 2019 for OBS
 # end at 31 Dec 2100 for models
 # first case: time array is numpy.datetime64
 if np.issubdtype(da.time.dtype, np.datetime64):
  da = da.where(da.time >= np.datetime64('1982-01-01T00:00:00'), drop=True)
  if m == 0:    # end at 31 Dec 2019 for OBS
   da = da.where(da.time <= np.datetime64('2019-12-31T23:59:59'), drop=True)
  else:         # end at 31 Dec 2100 for models
   da = da.where(da.time <= np.datetime64('2100-12-31T23:59:59'), drop=True)
   
  tlab1 = np.datetime_as_string(da.time.values[0], unit='Y')
  tlab2 = np.datetime_as_string(da.time.values[-1], unit='Y')
  date_start = np.datetime_as_string(da.time.values[0])
  date_end = np.datetime_as_string(da.time.values[-1])
 else:
  # second case: time array is 365Day calendar
  if type(da.time.values[0]) is cftime._cftime.DatetimeNoLeap:
   da = da.where(da.time >= cftime.DatetimeNoLeap(1982,1,1,0,0,0), drop=True)
   da = da.where(da.time <= cftime.DatetimeNoLeap(2100,12,31,23,59,59), drop=True)
  # third case: time array is 360Day calendar
  elif type(da.time.values[0]) is cftime._cftime.Datetime360Day:
   da = da.where(da.time >= cftime.Datetime360Day(1982,1,1,0,0,0), drop=True)
   da = da.where(da.time <= cftime.Datetime360Day(2100,12,30,23,59,59), drop=True)
  
  tlab1 = da.time.values[0].strftime('%Y')
  tlab2 = da.time.values[-1].strftime('%Y')
  date_start = da.time.values[0].strftime()
  date_end = da.time.values[-1].strftime()
 
 # year span label
 tlab = tlab1 + '-' + tlab2
 print(tlab)
 print('Start: ' + date_start)
 print('End: ' + date_end)
 
 # clear all variables except 'tos'
 varlist = list(ds.data_vars)   # first get the variable list
 varlist.remove('tos')  # remove 'tos' from the list, which will be kept
 da = da.drop_vars(varlist)
 
 # remove zlev parameter from OBS data
 if m == 0:
  da = da.sel(zlev=0, drop=True)
 
 # load required data to memory
 print('Reading data to memory... ')
 start = time.time()
 sst = da.tos.values
 time_arr = da.time.values
 lat = da.lat.values
 lon = da.lon.values
 end = time.time()
 print(end - start)
 
 sst.shape   # check dimensions of array [time,lat,lon]
 
 # convert time array to ordinal time array for mhw.detect
 Ly_set = False   # required for mhw.detect, set to 360 for Datetime360Day
 if np.issubdtype(time_arr.dtype, np.datetime64):
  # handle time arrays that are numpy.datetime64
  yr1=pandas.to_datetime(time_arr[0]).year
  yr2=pandas.to_datetime(time_arr[-1]).year
  time_o = [0] * len(time_arr)
  for i in range(0,len(time_arr)):
   time_o[i] = pandas.Timestamp(time_arr[i]).toordinal() # convert to pandas time
 else:
  # this follows ECO's handling of time arrays for CMIP5
  yr1=time_arr[0].year
  yr2=time_arr[-1].year
  t1=[time_arr[0].year,time_arr[0].month,time_arr[0].day]
  t2=[time_arr[-1].year,time_arr[-1].month,time_arr[-1].day]
  time_o,tmp,T,year,month,day,doy = eco.timevector(t1, t2)
  if type(time_arr[0]) is cftime._cftime.DatetimeNoLeap:
   feb29s = ~((month==2) * (day==29))
   time_o = time_o[feb29s]
  elif type(time_arr[0]) is cftime._cftime.Datetime360Day:
   Ly_set = 360
   feb29s = ~((month==2) * (day==29))
   time_o = time_o[feb29s]
   month = month[feb29s]
   day = day[feb29s]
   y360 = ~((day > 30) * (month > 3))
   time_o = time_o[y360]
 
 # convert time array (list) to a numpy array (necessary for 'detect')
 time_o=np.array(time_o)
 
 # initialise variables for MHW statistics
 time_yr = list(range(yr1,yr2+1))
 cat = list(range(1,5))
 i_which = range(sst.shape[2])  # lon dimension
 j_which = range(sst.shape[1])  # lat dimension
 mhw_total_events = np.NaN*np.zeros((len(j_which), len(i_which)))
 mhw_cats_dpy = np.NaN*np.zeros((len(time_yr),len(cat),len(j_which), len(i_which)))
 sst_mean = np.NaN*np.zeros((len(time_yr),len(j_which), len(i_which)))
 
 # for observed data, initialise additional variables
 if m == 0:
  mhw_count = np.NaN*np.zeros((len(time_yr),len(j_which), len(i_which)))
  mhw_intensity = np.NaN*np.zeros((len(time_yr),len(j_which), len(i_which)))
  mhw_duration = np.NaN*np.zeros((len(time_yr),len(j_which), len(i_which)))
  mhw_count_tr = np.NaN*np.zeros((len(j_which), len(i_which)))
  mhw_intensity_tr = np.NaN*np.zeros((len(j_which), len(i_which)))
  mhw_duration_tr = np.NaN*np.zeros((len(j_which), len(i_which)))
 
 # loop over every gridpoint, and compute MHW statistics
 for i in i_which:
  start = time.time()
  print(i, 'of', len(i_which)-1)
  for j in j_which:
   # process single SST timeseries
   sst1 = sst[:,j,i]
   # skip cells with land, ice
   if np.logical_not(np.isfinite(sst1.sum())) or (sst1<-1).sum()>0:
    mhw_total_events[j,i] = 0
   else:
    # detect MHWs
    mhws, clim = mhw.detect(time_o, sst1, climatologyPeriod=clim_b, pctile=90,Ly=Ly_set)
    # perform annual averaging of statistics
    mhwBlock = mhw.blockAverage(time_o, mhws, clim, temp=sst1)
    # store total MHW counts
    mhw_total_events[j,i] = mhwBlock['count'].sum()
    # store days per year (dpy) in each MHW category
    mhw_cats_dpy[:,0,j,i] = mhwBlock['moderate_days']
    mhw_cats_dpy[:,1,j,i] = mhwBlock['strong_days']
    mhw_cats_dpy[:,2,j,i] = mhwBlock['severe_days']
    mhw_cats_dpy[:,3,j,i] = mhwBlock['extreme_days']
    # store annual mean SST
    sst_mean[:,j,i] = mhwBlock['temp_mean']
    
    # for observed data, store additional variables
    if m == 0:
     mhw_count[:,j,i] = mhwBlock['count']
     mhw_intensity[:,j,i] = mhwBlock['intensity_max']  # annual mean of max MHW intensity
     mhw_duration[:,j,i] = mhwBlock['duration']
     # mean and trend
     mean, trend, dtrend = mhw.meanTrend(mhwBlock)
     # store trend data
     mhw_count_tr[j,i] = trend['count']
     mhw_intensity_tr[j,i] = trend['intensity_max']
     mhw_duration_tr[j,i] = trend['duration']
    
  end = time.time()
  print(end - start)
 
 # create xarray dataset for processed results
 ds_out = da.drop_vars({'tos'})   # follow format of input dataset, without tos
 dim = da.tos.dims                # read dimension names
 ds_out.attrs = {}                # clear attributes
 ds_out['time'] = time_yr         # set new time coordinate
 ds_out = ds_out.assign_coords({"cat" : cat})    # set category coordinate
 
 # store new variables in dataset
 ds_out['sst_mean'] = (('time', dim[-2], dim[-1]), sst_mean)
 ds_out['mhw_total_events'] = ((dim[-2], dim[-1]), mhw_total_events)
 ds_out['mhw_cats_dpy'] = (('time', 'cat', dim[-2], dim[-1]), mhw_cats_dpy)
 
 # for observed data, store additional variables
 if m == 0:
  ds_out['mhw_count'] = (('time', dim[-2], dim[-1]), mhw_count)
  ds_out['mhw_intensity'] = (('time', dim[-2], dim[-1]), mhw_intensity)
  ds_out['mhw_duration'] = (('time', dim[-2], dim[-1]), mhw_duration)
  ds_out['mhw_count_tr'] = ((dim[-2], dim[-1]), mhw_count_tr)
  ds_out['mhw_intensity_tr'] = ((dim[-2], dim[-1]), mhw_intensity_tr)
  ds_out['mhw_duration_tr'] = ((dim[-2], dim[-1]), mhw_duration_tr)
 
 # set coverage_content_type: https://wiki.esipfed.org/Concepts_Glossary
 if m == 0:
  cct = "physicalMeasurement"
 else:
  cct = "modelResult"
 
 # set attributes of variables
 ds_out['time'] = ds_out.time.assign_attrs(units="years", standard_name="time", long_name="calendar year", axis = "T", calendar="proleptic_gregorian")
 ds_out['cat'] = ds_out.cat.assign_attrs(units="1", long_name="Marine heatwave category following Hobday et al. (2018) definition", categories="1: Moderate, 2: Strong, 3: Severe, 4: Extreme")
 ds_out['sst_mean'] = ds_out.sst_mean.assign_attrs(units=da.tos.attrs['units'], standard_name="sea_surface_temperature", long_name="Annual mean sea surface temperature", coverage_content_type=cct)
 ds_out['mhw_total_events'] = ds_out.mhw_total_events.assign_attrs(units="1", standard_name="n/a", long_name="Total number of marine heatwaves detected", coverage_content_type="auxiliaryInformation")
 ds_out['mhw_cats_dpy'] = ds_out.mhw_cats_dpy.assign_attrs(units="1", standard_name="n/a", long_name="Count of days per year (dpy) in each marine heatwave category", coverage_content_type="auxiliaryInformation")
 
 # for observed data, set attributes of additional variables
 if m == 0:
  ds_out['mhw_count'] = ds_out.mhw_count.assign_attrs(units="1", standard_name="n/a", long_name="Count of marine heatwave events in each year", coverage_content_type="auxiliaryInformation")
  ds_out['mhw_intensity'] = ds_out.mhw_intensity.assign_attrs(units=da.tos.attrs['units'], standard_name="n/a", long_name="Annual mean of maximum marine heatwave intensities in each year (as an anomaly w.r.t. seasonal climatology)", coverage_content_type="auxiliaryInformation")
  ds_out['mhw_duration'] = ds_out.mhw_duration.assign_attrs(units="1", standard_name="n/a", long_name="Mean duration (in days) of marine heatwave events in each year", coverage_content_type="auxiliaryInformation")
  ds_out['mhw_count_tr'] = ds_out.mhw_count_tr.assign_attrs(units="1", standard_name="n/a", long_name="Trend in annual counts of marine heatwave events (events/year)", coverage_content_type="auxiliaryInformation")
  ds_out['mhw_intensity_tr'] = ds_out.mhw_intensity_tr.assign_attrs(units=da.tos.attrs['units'], standard_name="n/a", long_name="Trend in annual mean of maximum of marine heatwave intensities (degree_C/year)", coverage_content_type="auxiliaryInformation")
  ds_out['mhw_duration_tr'] = ds_out.mhw_duration_tr.assign_attrs(units="1", standard_name="n/a", long_name="Trend in annual mean duration of marine heatwaves (days/year)", coverage_content_type="auxiliaryInformation")
 
 # set global attributes
 ds_out.attrs['source_code'] = "https://github.com/jbkajtar/mhw_pacific"
 ds_out.attrs['title'] = "Marine heatwave statistics for the Pacific Islands region (120E-140W, 40S-15N)"
 ds_out.attrs['summary'] = "Data generated for Holbrook et al., 'Impacts of marine heatwaves on tropical western and central Pacific island nations and their communities', Glob Planet Change, (2021)"
 ds_out.attrs['source_data'] = modcode
 ds_out.attrs['keywords'] = "marine heatwave; extreme event; impact; ocean warming; Pacific; CMIP6 projections"
 ds_out.attrs['reference climatology'] = str(clim_b[0]) + '-' + str(clim_b[1])
 ds_out.attrs['Conventions'] = "ACDD-1.3"
 
 # save dataset to netcdf file
 print('Saving data... ')
 start = time.time()
 ds_out.to_netcdf(outfile + '.nc')
 end = time.time()
 print(end - start)


