# JBK 2021-08-21
# process marine heatwave statistics for Pacific Island case study regions
# 1. read SST indices from observations and CMIP6 (historical and future)
# 2. compute MHW statistics
# 3. store required data
# OUTPUT: mhw_ts/mhw_ts.pac_is.<model>.<scenario>.<variant>.nc

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
inpath = ''  # location of SST indices
outpath = '' # output directory

# define which scenarios to use
scen_h = 'historical'
scen_f = 'ssp585'        # ssp126 or ssp585

# define bounds of stored region
reg_lab = 'pac_is'

# define climatological baseline period
clim_b = [1995,2014]

# read list of model files
scen = scen_h + '+' + scen_f
filenames = glob.glob(inpath + '*' + scen + '*.nc')
filenames_o = inpath + 'sst_indices.pac_is.NOAA_OISST.AVHRR.v2-1_modified.nc'
filenames.insert(0, filenames_o)  # insert OBS at start of list
nmods = len(filenames)

# range of models to process in this run
m_which = range(0,nmods)

# loop over all models in list
for m in m_which:
 
 # split path to retrieve filename
 modcode = filenames[m].split("/")

 # read model code
 modcode = str.replace(modcode[-1], 'sst_indices.' + reg_lab + '.', '')
 modcode = str.replace(modcode, '.nc', '')
 
 # set input filename
 infiles = filenames[m]
 
 # create output filename
 outfile = outpath + 'mhw_ts.' + reg_lab + '.' + modcode

 # load data
 print('Processing... ' + modcode)
 ds = xr.open_mfdataset(infiles,combine='nested',concat_dim='time')
 
 # load required data to memory
 print('Reading data to memory... ')
 sst = ds.tos.values
 time_arr = ds.time.values
 
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
 i_which = range(sst.shape[1])  # number of indices to process
 mhw_total_events = np.zeros(len(i_which))
 mhw_cats_dpy = np.zeros((len(time_yr),len(cat),len(i_which)))
 sst_mean = np.NaN*np.zeros((len(time_yr),len(i_which)))
 
 # loop over every gridpoint, and compute MHW statistics
 for i in i_which:
  print(i, 'of', len(i_which)-1)

  # process single SST timeseries
  sst1 = sst[:,i]

  # detect MHWs
  mhws, clim = mhw.detect(time_o, sst1, climatologyPeriod=clim_b,pctile=90,Ly=Ly_set)
  # perform annual averaging of statistics
  mhwBlock = mhw.blockAverage(time_o, mhws, clim, temp=sst1)
  # store total MHW counts
  mhw_total_events[i] = mhwBlock['count'].sum()
  # store days per year (dpy) in each MHW category
  mhw_cats_dpy[:,0,i] = mhwBlock['moderate_days']
  mhw_cats_dpy[:,1,i] = mhwBlock['strong_days']
  mhw_cats_dpy[:,2,i] = mhwBlock['severe_days']
  mhw_cats_dpy[:,3,i] = mhwBlock['extreme_days']
  # store annual mean SST
  sst_mean[:,i] = mhwBlock['temp_mean']
 
 # convert count arrays to integer arrays
 mhw_total_events = mhw_total_events.astype(int)
 mhw_cats_dpy = mhw_cats_dpy.astype(int)
 
 # create xarray dataset for processed results
 ds_out = ds.drop_vars({'tos'})   # follow format of input dataset, without tos
 dim = ds.tos.dims                # read dimension names
 ds_out['time'] = time_yr         # set new time coordinate
 ds_out = ds_out.assign_coords({"cat" : cat})    # set category coordinate
 
 # store new variables in dataset
 ds_out['sst_mean'] = (('time', 'region'), sst_mean)
 ds_out['mhw_total_events'] = ('region', mhw_total_events)
 ds_out['mhw_cats_dpy'] = (('time', 'cat', 'region'), mhw_cats_dpy)
 
 # set coverage_content_type: https://wiki.esipfed.org/Concepts_Glossary
 if m == 0:
  cct = "physicalMeasurement"
 else:
  cct = "modelResult"
 
 # set attributes of variables
 ds_out['time'] = ds_out.time.assign_attrs(units="years", standard_name="time", long_name="calendar year", axis = "T", calendar="proleptic_gregorian")
 ds_out['cat'] = ds_out.cat.assign_attrs(units="1", long_name="Marine heatwave category following Hobday et al. (2018) definition", categories="1: Moderate, 2: Strong, 3: Severe, 4: Extreme")
 ds_out['sst_mean'] = ds_out.sst_mean.assign_attrs(units=ds.tos.attrs['units'], standard_name="sea_surface_temperature", long_name="Annual mean sea surface temperature", coverage_content_type=cct)
 ds_out['mhw_total_events'] = ds_out.mhw_total_events.assign_attrs(units="1", standard_name="n/a", long_name="Total number of marine heatwaves detected", coverage_content_type="auxiliaryInformation")
 ds_out['mhw_cats_dpy'] = ds_out.mhw_cats_dpy.assign_attrs(units="1", standard_name="n/a", long_name="Count of days per year (dpy) in each marine heatwave category", coverage_content_type="auxiliaryInformation")
 
 # set global attributes
 ds_out.attrs['source_code'] = "https://github.com/jbkajtar/mhw_pacific"
 ds_out.attrs['title'] = "Marine heatwave statistics for the Pacific Islands case-study regions: Fiji, Samoa, Palau"
 ds_out.attrs['summary'] = "Data generated for Holbrook et al., 'Impacts of marine heatwaves on tropical western and central Pacific island nations and their communities', Glob Planet Change, (2021)"
 ds_out.attrs['source_data'] = modcode
 ds_out.attrs['keywords'] = "marine heatwave; extreme event; impact; ocean warming; Pacific; CMIP6 projections"
 ds_out.attrs['Conventions'] = "ACDD-1.3"
 
 # save dataset to netcdf file
 print('Saving data... ')
 ds_out.to_netcdf(outfile + '.nc')


