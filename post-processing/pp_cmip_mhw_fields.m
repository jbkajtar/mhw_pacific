% JBK 2021-06-08
% post-process MHW statistics from CMIP6 to common grid
% use observed grid as common grid
% interpolate using nearest neighbour, to preserve model resolution
% of CMIP6 statistics, without smoothing
% read model data, pre-process, re-grid

clear all;

script_name='p42/pp_cmip_mhw_fields.m';

% region bounds
reg_lab='pac_is';
reg_bounds=[120, 220, -40, 15];
scen='historical+ssp126';

% set paths
sourcepath=''; % source path of MHW data (pac_mhw_stats)
inpath=[sourcepath 'mhw_cats/'];
outpath=[sourcepath];

% set periods of blocks
yr_start=1982;
yr_blk=[1982,2001; 2000,2019; 2020,2039; 2040,2059];
clear yr_str;
for jj=1:size(yr_blk,1)
 yr_str{jj}=[num2str(yr_blk(jj,1)) '-' num2str(yr_blk(jj,2))];
end

% list of models to process
dirlist=dir([inpath 'mhw_cats.' reg_lab '.*.' scen '*.nc']);
modcode={dirlist.name};
modcode=strrep(modcode,'mhw_cats.pac_is.','');
modcode=strrep(modcode,'.nc','');
modcode=['NOAA_OISST.AVHRR.v2-1_modified' modcode];  % append OBS dataset
mmax=length(modcode);

for m=1:mmax

 % read data
 clear ds;
 infile=[inpath 'mhw_cats.pac_is.' modcode{m} '.nc'];
 ds.lon=double(ncread(infile,'lon'));
 ds.lat=double(ncread(infile,'lat'));
 ds.mhw_cats_dpy=double(ncread(infile,'mhw_cats_dpy'));
 ds.sst_mean=double(ncread(infile,'sst_mean'));
 ds.time=double(ncread(infile,'time'));

 % set number of blocks to process: only 2 for obs
 nb=4;
 if contains(modcode{m},'NOAA_OISST')
  nb=2;
 end

 % compute the mean in each block
 for ii=1:nb
  [t1,t2]=findrange(ds.time,yr_blk(ii,1),yr_blk(ii,2));
  % compute mean SST in each block
  ds.sst_mean_blk(:,:,ii)=mean(ds.sst_mean(:,:,t1:t2),3);
  for kk=1:4
   % the moderate day counts should also include all more extreme category days
   ds.mhw_cats_blk(:,:,kk,ii)=mean(sum(ds.mhw_cats_dpy(:,:,kk:4,t1:t2),3),4);
  end
 end

 % model label
 if m==1
  model_id='NOAA_OISST';
 else
  model_id=strrep(modcode{m},['.' scen],'');
  model_id=regexprep(model_id,'.r\w*i\w*p\w*f\w*','');
 end
 disp(model_id);
 modname{m}=model_id;

 % rotate coordinates, if required
 ds.lon(ds.lon<0)=ds.lon(ds.lon<0)+360;
 
 % store observations in first position of global array
 if m==1
  lat=ds.lat;
  lon=ds.lon;
  [yq,xq]=meshgrid(lat,lon);
  mhw_cats=ds.mhw_cats_blk;
  sst_mean=ds.sst_mean_blk;
 else
  % use meshgrid to create 2D arrays if lon/lat is 1D
  if any(size(ds.lon)==1)
   [ds.lat,ds.lon]=meshgrid(ds.lat,ds.lon);
  end
  
  % regrid models onto observed grid
  for ii=1:nb
   sst_mean(:,:,ii,m)=griddata(ds.lon,ds.lat,ds.sst_mean_blk(:,:,ii),xq,yq,'nearest');
   for kk=1:4
    mhw_cats(:,:,kk,ii,m)=griddata(ds.lon,ds.lat,ds.mhw_cats_blk(:,:,kk,ii),xq,yq,'nearest');
   end
  end
  
 end
end

% additional variables for netCDF storage
model_name=char(modname);
categories=[1:4];
time_block=yr_blk;

% write to netCDF
f1=[outpath 'mhw_cats.pac_is.cmip6.' scen '.nc'];
fmt='netcdf4_classic';

% save data to netcdf
nccreate(f1,'lon', 'Dimensions',{'lon',length(lon)}, 'Datatype','single','Format',fmt);
ncwrite(f1,'lon',lon);
ncwriteatt(f1,'lon','standard_name','longitude');
ncwriteatt(f1,'lon','long_name','Longitude');
ncwriteatt(f1,'lon','units','degrees_east');
ncwriteatt(f1,'lon','axis','X');

nccreate(f1,'lat', 'Dimensions',{'lat',length(lat)}, 'Datatype','single','Format',fmt);
ncwrite(f1,'lat',lat);
ncwriteatt(f1,'lat','standard_name','latitude');
ncwriteatt(f1,'lat','long_name','Latitude');
ncwriteatt(f1,'lat','units','degrees_north');
ncwriteatt(f1,'lat','axis','Y');

nccreate(f1,'cat', 'Dimensions',{'cat',length(categories)}, 'Datatype','single','Format',fmt);
ncwrite(f1,'cat',categories);
ncwriteatt(f1,'cat','units','1');
ncwriteatt(f1,'cat','long_name','Marine heatwave category following Hobday et al. (2018) definition');
ncwriteatt(f1,'cat','categories','1: Moderate, 2: Strong, 3: Severe, 4: Extreme');

nccreate(f1,'time_block', 'Dimensions',{'time_block',size(time_block,1),'time_bnds',size(time_block,2)}, 'Datatype','single','Format',fmt);
ncwrite(f1,'time_block',time_block);
ncwriteatt(f1,'time_block','units','years')
ncwriteatt(f1,'time_block','standard_name','time');
ncwriteatt(f1,'time_block','long_name','Bounds of time averaging blocks');
ncwriteatt(f1,'time_block','axis','T');

nccreate(f1,'model_name', 'Dimensions',{'model_name',size(model_name,1),'charlen',size(model_name,2)}, 'Datatype','char', 'Format',fmt);
ncwrite(f1,'model_name',model_name);
ncwrite(f1,'model_name',model_name);
ncwriteatt(f1,'model_name','units','1');
ncwriteatt(f1,'model_name','long_name','Names of observational (NOAA_OISST) or model data sources');

nccreate(f1,'mhw_cats', 'Dimensions',{'lon',length(lon),'lat',length(lat),'cat',length(categories), 'time_block',size(time_block,1),'model_name',size(model_name,1)}, 'Datatype','single', 'Format',fmt, 'DeflateLevel',2);
ncwrite(f1,'mhw_cats',mhw_cats);
ncwriteatt(f1,'mhw_cats','units','1');
ncwriteatt(f1,'mhw_cats','standard_name','n/a');
ncwriteatt(f1,'mhw_cats','long_name','Count of days per year (dpy) in each marine heatwave category, averaged in each time block');
ncwriteatt(f1,'mhw_cats','coverage_content_type','auxiliaryInformation');

nccreate(f1,'sst_mean', 'Dimensions',{'lon',length(lon),'lat',length(lat),'time_block',size(time_block,1), 'model_name',size(model_name,1)}, 'Datatype','single', 'Format',fmt, 'DeflateLevel',2);
ncwrite(f1,'sst_mean',sst_mean);
ncwriteatt(f1,'sst_mean','units','degree_C');
ncwriteatt(f1,'sst_mean','standard_name','sea_surface_temperature');
ncwriteatt(f1,'sst_mean','long_name','Annual mean sea surface temperature, averaged in each time block');
ncwriteatt(f1,'sst_mean','coverage_content_type','modelResult');

ncwriteatt(f1,'/','source_code','https://github.com/jbkajtar/mhw_pacific');
ncwriteatt(f1,'/','title','Marine heatwave statistics for the Pacific Islands region (120E-140W, 40S-15N), averaged in 20-year blocks');
ncwriteatt(f1,'/','summary','Data generated for Holbrook et al., ''Impacts of marine heatwaves on tropical western and central Pacific island nations and their communities'', Glob Planet Change, (2021)');
ncwriteatt(f1,'/','keywords','marine heatwave; extreme event; impact; ocean warming; Pacific; CMIP6 projections');
ncwriteatt(f1,'/','Conventions','ACDD-1.3');

