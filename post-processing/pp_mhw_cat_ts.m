% JBK 2021-06-21
% post-process Pacific Islands case-study regions indices
% store MHW indices in single file, for plotting purposes
% output: mhw_cats.pac_is.cmip6.timeseries.nc

clear all;

script_name='p42/pp_mhw_cat_ts.m';

% set paths
sourcepath='';   % source path of MHW data (pac_mhw_stats)
inpath=[sourcepath 'mhw_ts/'];
outpath=[sourcepath];

% set required variables
yr_start=1982;
reg_lab='pac_is';
scen={'historical+ssp126','historical+ssp585'};

% store bounds of case study regions
% Fiji
ii=1;
idxcode{ii}='fiji';
reg_bnds(:,ii)=[174 182 -21 -14];

% Samoa
ii=2;
idxcode{ii}='samoa';
reg_bnds(:,ii)=[186 190 -15 -12];

% Palau
ii=3;
idxcode{ii}='palau';
reg_bnds(:,ii)=[130 136 2 9];
imax=length(idxcode);

% list of models to process
dirlist=dir([inpath 'mhw_ts.' reg_lab '.*.' scen{1} '*.nc']);
modcode={dirlist.name};
modcode=strrep(modcode,['mhw_ts.' reg_lab '.'],'');
modcode=strrep(modcode,'.nc','');
modcode=['NOAA_OISST.AVHRR.v2-1_modified' modcode];  % append OBS dataset
mmax=length(modcode);

% store model names
modname=strrep(modcode,['.AVHRR.v2-1_modified'],'');
modname=strrep(modname,['.' scen{1}],'');
modname=regexprep(modname,'.r\w*i\w*p\w*f\w*','');

% initialise global arrays
time=[yr_start:2100]';
for ii=1:imax
 idx.(idxcode{ii})=nan(length(time),4,mmax,2);
 sst.(idxcode{ii})=nan(length(time),mmax,2);
end

% read data
for iscen=1:2
 for m=1:mmax
  model_id=modcode{m};
  if iscen==2
   model_id=strrep(model_id,scen{1},scen{2});
  end
  disp(model_id);

  clear ds;
  infile=[inpath 'mhw_ts.' reg_lab '.' model_id '.nc'];
  ds.mhw_cats_dpy=double(ncread(infile,'mhw_cats_dpy'));
  ds.sst_mean=double(ncread(infile,'sst_mean'));
  time1=double(ncread(infile,'time'));

  % the moderate day counts should also include all more extreme days
  clear mhw_cats;
  for kk=1:4
   mhw_cats(:,kk,:)=sum(ds.mhw_cats_dpy(:,kk:4,:),2);
  end
 
  % trim to available time span
  [t1,t2]=findrange(time,time1(1),time1(end));

  % loop over all index regions
  for ir=1:imax
   
   % store annual mean SST in global array
   sst1=ds.sst_mean(ir,:);
   sst.(idxcode{ir})(t1:t2,m,iscen)=sst1;
   
   % area average MHW days in each category
   for ic=1:4
    mhw1=squeeze(mhw_cats(ir,ic,:));
    idx.(idxcode{ir})(t1:t2,ic,m,iscen)=mhw1;
   end
  end

 end
end

% additional variable changes for netCDF storage
model_name=char(modname);
scen=char(scen);
region=char(idxcode);
categories=[1:4];
clear mhw_cats sst_idx;
for ir=1:imax
 mhw_cats(:,:,:,:,ir)=idx.(idxcode{ir});
 sst_idx(:,:,:,ir)=sst.(idxcode{ir});
end

% write to netCDF
f1=[outpath 'mhw_cats.pac_is.cmip6.timeseries.nc'];
fmt='netcdf4_classic';

% save data to netcdf
nccreate(f1,'time', 'Dimensions',{'time',length(time)}, 'Datatype','single','Format',fmt);
ncwrite(f1,'time',time);
ncwriteatt(f1,'time','units','years')
ncwriteatt(f1,'time','standard_name','time');
ncwriteatt(f1,'time','long_name','calendar year');
ncwriteatt(f1,'time','axis','T');
ncwriteatt(f1,'time','calendar','proleptic_gregorian');

nccreate(f1,'cat', 'Dimensions',{'cat',length(categories)}, 'Datatype','single','Format',fmt);
ncwrite(f1,'cat',categories);
ncwriteatt(f1,'cat','units','1');
ncwriteatt(f1,'cat','long_name','Marine heatwave category following Hobday et al. (2018) definition');
ncwriteatt(f1,'cat','categories','1: Moderate, 2: Strong, 3: Severe, 4: Extreme');

nccreate(f1,'model_name', 'Dimensions',{'model_name',size(model_name,1),'charlen',size(model_name,2)}, 'Datatype','char', 'Format',fmt);
ncwrite(f1,'model_name',model_name);
ncwriteatt(f1,'model_name','units','1');
ncwriteatt(f1,'model_name','long_name','Names of observational (NOAA_OISST) or model data sources');

nccreate(f1,'scen', 'Dimensions',{'scen',size(scen,1),'charlen2',size(scen,2)}, 'Datatype','char', 'Format',fmt);
ncwrite(f1,'scen',scen);
ncwriteatt(f1,'scen','units','1');
ncwriteatt(f1,'scen','long_name','Names of scenarios, historical plus low emissions (ssp126) or high emissions (ssp585) extension');

nccreate(f1,'region', 'Dimensions',{'region',size(region,1),'charlen3',size(region,2)}, 'Datatype','char', 'Format',fmt);
ncwrite(f1,'region',region);
ncwriteatt(f1,'region','units','1');
ncwriteatt(f1,'region','long_name','Name of case study regions');

nccreate(f1,'reg_bnds', 'Dimensions',{'reg_verts',size(reg_bnds,1),'region',size(reg_bnds,2)}, 'Datatype','single', 'Format',fmt);
ncwrite(f1,'reg_bnds',reg_bnds);
ncwriteatt(f1,'reg_bnds','units','degrees');
ncwriteatt(f1,'reg_bnds','standard_name','n/a');
ncwriteatt(f1,'reg_bnds','long_name','Longitude (degrees east) and latitude (degrees north) bounds of case study regions');
ncwriteatt(f1,'reg_bnds','coverage_content_type','auxiliaryInformation');

nccreate(f1,'mhw_cats', 'Dimensions',{'time',length(time),'cat',length(categories), 'model_name',size(model_name,1),'scen',size(scen,1),'region',size(region,1)}, 'Datatype','single', 'Format',fmt, 'DeflateLevel',2);
ncwrite(f1,'mhw_cats',mhw_cats);
ncwriteatt(f1,'mhw_cats','units','1');
ncwriteatt(f1,'mhw_cats','standard_name','n/a');
ncwriteatt(f1,'mhw_cats','long_name','Count of days per year (dpy) in each marine heatwave category');
ncwriteatt(f1,'mhw_cats','coverage_content_type','auxiliaryInformation');

nccreate(f1,'sst_idx', 'Dimensions',{'time',length(time),'model_name',size(model_name,1), 'scen',size(scen,1),'region',size(region,1)}, 'Datatype','single', 'Format',fmt, 'DeflateLevel',2);
ncwrite(f1,'sst_idx',sst_idx);
ncwriteatt(f1,'sst_idx','units','degree_C');
ncwriteatt(f1,'sst_idx','standard_name','sea_surface_temperature');
ncwriteatt(f1,'sst_idx','long_name','Annual mean sea surface temperature');
ncwriteatt(f1,'sst_idx','coverage_content_type','modelResult');

ncwriteatt(f1,'/','source_code','https://github.com/jbkajtar/mhw_pacific');
ncwriteatt(f1,'/','title','Marine heatwave statistics for the Pacific Islands case-study regions: Fiji, Samoa, Palau');
ncwriteatt(f1,'/','summary','Data generated for Holbrook et al., ''Impacts of marine heatwaves on tropical western and central Pacific island nations and their communities'', Glob Planet Change, (2021)');
ncwriteatt(f1,'/','keywords','marine heatwave; extreme event; impact; ocean warming; Pacific; CMIP6 projections');
ncwriteatt(f1,'/','Conventions','ACDD-1.3');


