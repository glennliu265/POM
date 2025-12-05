#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Regrid REI calculated from calculate_REI_general
- Computed from [calculate_REI_general] in re-emergence folder
- Uses ACF computed by [pointwise_crosscorrelation] (also in re-emergence folder)

Created on Mon Nov 17 13:55:19 2025

@author: gliu
"""

import xesmf as xe
import numpy as np
import matplotlib.pyplot as plt
import glob
import xarray as xr
import cartopy.crs as ccrs
import tqdm

#%% User Edits

datpath = "/Users/gliu/Downloads/02_Research/01_Projects/06_POM/01_Data/proc/acfs/"
outpath = datpath

ncnames = [
    "SOM_006-0360_REI.nc",
    "MCOM_0100-0599_REI.nc",
    "FOM_1100-1599_REI.nc",
    ]

refpath = "/Users/gliu/Downloads/02_Research/01_Projects/06_POM/01_Data/proc/"
refnc   = refpath+"SOM_HBLT_clim_regrid_bilinear.nc"

#%% Function (copied from proc 2025.11.17)

def addstrtoext(name,addstr,adjust=0):
    """
    Add [addstr] to the end of a string with an extension [name.ext]
    Result should be "name+addstr+.ext"
    -4: 3 letter extension. -3: 2 letter extension
    """
    return name[:-(4+adjust)] + addstr + name[-(4+adjust):]

#%% Load Data

dsref  = xr.open_dataset(refnc).load()
ds_all = []
for nc in ncnames:
    ds = xr.open_dataset(datpath+nc).load()
    ds_all.append(ds)
    
#%% Load TLAT and TLON (unfortunately was not saved in rei)

datpath2="/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/01_Data/proc/"
ncname =datpath2+"FOM_1100-1599_lag00to60_ALL_ensALL.nc"
outname=outpath+"FOM_1100-1599_REI.nc"
tlon = xr.open_dataset(ncname).TLONG.load().squeeze()
tlat = xr.open_dataset(ncname).TLAT.load().squeeze()

#%% Set up Regridding Array
lat_out = dsref.lat.data 
lon_out = dsref.lon.data

ds_out    = xr.Dataset({'lat': (['lat'], lat_out), 'lon': (['lon'], lon_out) })

#%%



for ex in tqdm.tqdm(range(len(ds_all))):
    #
    ds = ds_all[ex]
    
    # Fix Dimensions
    dsmerge   = xr.merge([ds,tlat,tlon]) # Add tlat/tlon
    dsmerge   = dsmerge.rename(dict(lon='nlon',lat='nlat')) # Rename lat/lon to avoid conflict
    dsmerge   = dsmerge.rename({"TLONG": "lon", "TLAT": "lat"})

    regridder = xe.Regridder(dsmerge, ds_out, 'bilinear',periodic=True)
    
    dsregrid  = regridder(dsmerge)
    
    newname   = datpath + addstrtoext(ncnames[ex],"_regridded",adjust=-1)
    edict     = {'rei' : {'zlib': True}}
    dsregrid.to_netcdf(newname,encoding=edict)