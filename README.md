# Tropical western and central Pacific marine heatwave data (mhw_pacific)

Code for generating marine heatwave (MHW) metrics from observational and CMIP6 data for Holbrook et al., Glob. Planet. Change (2021). A description of the data produced by this code is given in Kajtar et al., Data Br. (2021). Marine heatwaves are calculated following the Hobday et al. (2016) definition, using the marineHeatWaves python module (https://github.com/ecjoliver/marineHeatWaves).

This code has been written to generate marine heatwave metrics for the tropical western and central Pacific Ocean region (120°E-140°W, 40°S-15°N). The metrics are computed from daily sea surface temperature (SST) data, from both observations and models. The observed marine heatwave data are calculated from NOAA 0.25° daily Optimum Interpolation Sea Surface Temperature (OISST) over the period 1982-2019. The modelled marine heatwave data are from analysis of 18 model simulations as part of the Coupled Model Intercomparison Project, Phase 6 (CMIP6) over the period 1982-2100, where two future scenarios have been analysed. The marine heatwave data are provided on a grid point basis across the domain. Marine heatwave timeseries metrics are also provided for three case study regions: Fiji, Samoa, and Palau.

Further details about computing the marine heatwave metrics, and obtaining the source data, are given in Kajtar et al. (2021).

The data is publicly and freely available to download at https://doi.org/10.5281/zenodo.5069012.

# Contents

|Directory         |Description|
|------------------|-----------|
|code              |Code for computing MHW metrics|
|marineHeatWaves   |Customised MHW detection code|
|post-processing   |Code for post-processing data|
|sources           |Primary data source list|

# code

Python code for computing MHW metrics. Primary data source paths must first be specified in the files obtained from the 'sources' directory.

|File              |Description|
|------------------|-----------|
|compute_pacific_mhw_stats.py   |Compute MHW metrics across the Pacific Islands region|
|store_pacific_areacello.py     |Pre-process CMIP6 native grid data|
|compute_pacific_sst_indices.py |Compute area-average SST timeseries for case study regions|
|compute_pacific_mhw_ts.py      |Compute MHW metrics for case study regions|

# marineHeatWaves

Marine heatwaves detection code, available at: https://github.com/ecjoliver/marineHeatWaves. The code provided here has been slightly adapted, but the original module should also work.
To run the mhw_pacific code, an additional module is required, eoliver.py, available at: https://github.com/ecjoliver/TasmanSeaMHW_201516/blob/master/ecoliver.py

# post-processing

MATLAB code for post-processing of MHW metrics.

|File              |Description|
|------------------|-----------|
|findrange.m           |Function to find a range in an array|
|pp_cmip_mhw_fields.m  |Post-processing of MHW field data|
|pp_mhw_cat_ts.m       |Post-processing of MHW timeseries data, for case study regions|

# sources

Files that specify the locations of the primary source data, to be adapted by the user. The mhw_pacific reads data from these specified paths. One path per line in each file. Also provided here are the full details of CMIP6 realisations analysed in Kajtar et al. (2021), along with primary source references.

|File              |Description|
|------------------|-----------|
|cmip6_data_references.csv  |CMIP6 source list|
|cmip6_gadi_areacello.txt   |Source list for areacello|
|cmip6_gadi_historical.txt  |Source list for historical data|
|cmip6_gadi_ssp126.txt      |Source list for SSP1-2.6 scenario data|
|cmip6_gadi_ssp585.txt      |Source list for SSP5-8.5 scenario data|

# References

Hobday, A.J., L. V. Alexander, S.E. Perkins-Kirkpatrick, D.A. Smale, S.C. Straub, E.C.J. Oliver, J.A. Benthuysen, M.T. Burrows, M.G. Donat, M. Feng, N.J. Holbrook, P.J. Moore, H.A. Scannell, A. Sen Gupta, T. Wernberg, A hierarchical approach to defining marine heatwaves, Prog. Oceanogr. 141 (2016) 227–238. https://doi.org/10.1016/j.pocean.2015.12.014.

Holbrook, N.J., V. Hernaman, S. Koshiba, J. Lako, J.B. Kajtar, P. Amosa, A. Singh, Impacts of marine heatwaves on tropical western and central Pacific Island nations and their communities, Glob. Planet. Change. (2021) In Press.

Kajtar, J.B., V. Hernaman, N.J. Holbrook, P. Petrelli, Tropical western and central Pacific marine heatwave data calculated from gridded sea surface temperature observations and CMIP6. Data Br. (2021) Submitted.

# Acknowledgements

We acknowledge the World Climate Research Programme's Working Group on Coupled Modelling, which is responsible for CMIP, and we thank the climate modelling groups for producing and making available their model output. CMIP6 model outputs were made available with the assistance of resources from the National Computational Infrastructure (NCI), which is supported by the Australian Government. NOAA High Resolution SST data were provided by the NOAA National Centers for Environmental Information. We acknowledge Eric Oliver for the use of his marineHeatWaves python module, freely available at https://github.com/ecjoliver/marineHeatWaves. 
