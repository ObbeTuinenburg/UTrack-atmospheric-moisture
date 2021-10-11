# UTrack-atmospheric-moisture
Utrecht atmospheric moisture tracking model

This is the source code for atmospheric moisture tracking model to track the flows moisture from evaporation to the next precipitation location based on simulations forced with ERA5 reanalysis data.
The tracking of atmospheric moisture run as a post-processing step of reanalysis data (as is done here) comes with a number of uncertainties. These uncertainties are both numerical and on process level.
The model code here enables to do the sensitivity experiments described in the paper in Hydrology and Earth System Science (HESS): Tuinenburg, O. A. and Staal, A.: Tracking the global flows of atmospheric moisture and associated uncertainties, Hydrol. Earth Syst. Sci., 24, 2419–2435, https://doi.org/10.5194/hess-24-2419-2020, 2020.

The code is written in C and tested under a (Ubuntu) linux distribution. The model needs several atmospheric variables as its forcing. In the current version, ECMWF ERA5 reanalysis data (https://cds.climate.copernicus.eu) is used for this. Python scripts have been provided to download these. Note that you will need to create an Copernicus data store account to download these data using the Python API (see https://cds.climate.copernicus.eu/api-how-to).
