# JBK 2021-06-11
# store Pacific region cell areas, required for area averaging
# 1. read grid data from observations and CMIP6
# 2. trim to required region
# 3. store required data
# OUTPUT: grid/areacello.pac_is.<model>.<scenario>.<variant>.nc

# load required modules
import numpy as np
import xarray as xr
import pandas
import glob
from datetime import date
import time
import cftime
import scipy.io as io

# set paths and filenames
inpath = ''  # location of model source list
outpath = '' # output directory

# define bounds of stored region
reg_lab = 'pac_is'
reg_bnds = [120, 220, -40, 15]   # 120E-140W, 40S-15N

# read model lists
# model list .txt files should have OBS path as first entry
modfile = open(inpath + 'cmip6_gadi_areacello.txt', "r")
modlist = modfile.read().splitlines()
modfile.close()

nmods = len(modlist)

# range of models to process in this run
m_which = range(1,nmods)
print('PROCESS LIST')
for m in m_which:
 print('[' + str(m) + ']: ' + modlist[m])

# loop over all models in list
for m in m_which:
 
 # split path names to retrieve model parameters
 modcode = modlist[m].split("/")
 
 if m == 0:
  # set OBS label, e.g. NOAA_OISST.AVHRR.v2-1_modified
  modcode = '.'.join(modcode[4:7])
 else:
  # set model label, e.g. <source_id>.<experiment>.<variant_label>
  modcode = '.'.join(modcode[8:11])

 print('Processing... ' + modcode)
 
 # read netcdf filenames from historical and scenario paths
 filenames = sorted(glob.glob(modlist[m] + '*.nc'))
 infiles = filenames
 
 # create output filename
 outfile = outpath + 'areacello.' + reg_lab + '.' + modcode
 
 # load data
 print('Loading data... ')
 start = time.time()
 ds = xr.open_mfdataset(infiles,combine='nested',concat_dim='time')
 end = time.time()
 print(end - start)
 
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
 
 # clear all variables except 'areacello'
 varlist = list(ds.data_vars)  # first get the variable list
 varlist.remove('areacello')   # remove 'areacello' from the list, which will be kept
 da = da.drop_vars(varlist)
 
 # remove zlev parameter from OBS data
 if m == 0:
  da = da.sel(zlev=0, drop=True)
 
 # create xarray dataset for processed results
 ds_out = da        # follow format of input dataset
 ds_out.attrs = {}  # clear attributes
 
 # set global attributes
 ds_out.attrs['source_code'] = "https://github.com/jbkajtar/mhw_pacific"
 ds_out.attrs['title'] = "Cell area data for the Pacific Islands region (120E-140W, 40S-15N)"
 ds_out.attrs['summary'] = "Data generated for Holbrook et al., 'Impacts of marine heatwaves on tropical western and central Pacific island nations and their communities', Glob Planet Change, (2021)"
 ds_out.attrs['source_data'] = modcode
 ds_out.attrs['keywords'] = "marine heatwave; extreme event; impact; ocean warming; Pacific; CMIP6 projections"
 ds_out.attrs['Conventions'] = "ACDD-1.3"
 
 # save dataset to netcdf file
 print('Saving data... ')
 start = time.time()
 ds_out.to_netcdf(outfile + '.nc')
 end = time.time()
 print(end - start)


