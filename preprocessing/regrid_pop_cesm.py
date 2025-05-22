#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Regrid from POP to CAM Grid (CESM)
based on quick_regrid.py from reemergence module

Created on Mon May  5 16:44:43 2025

@author: gliu

"""

import time
import numpy as np
import xarray as xr
from tqdm import tqdm
import xesmf as xe
import sys

import glob

#%% Functions (Copied from proc on 2025.05.05)
def addstrtoext(name,addstr,adjust=0):
    """
    Add [addstr] to the end of a string with an extension [name.ext]
    Result should be "name+addstr+.ext"
    -4: 3 letter extension. -3: 2 letter extension
    """
    return name[:-(4+adjust)] + addstr + name[-(4+adjust):]

def make_encoding_dict(ds,encoding_type='zlib'):
    if type(ds) == xr.core.dataarray.DataArray:
        vname         = ds.name
        encoding_dict = {vname : {encoding_type:True}}
    else:
        keys          = list(ds.keys())
        values        = ({encoding_type:True},) * len(keys)
        encoding_dict = { k:v for (k,v) in zip(keys,values)}
    return encoding_dict
# -------------
#%% User Edits
# -------------

# Indicate Machine
machine       = "Astraeus"

# What to regrid?
# vname          = "HMXL"
# ncsearch       = "/stormtrack/data4/glliu/01_Data/CESM2_PiControl/proc/FCM_HMXL_0100_2000_clim.nc"

# v
vname    = "hblt"
ncsearch = "/stormtrack/data4/glliu/01_Data/CESM2_PiControl/proc/SOM_HBLT_clim.nc"


# Indicate Regridding Options
concatdim      = 'time'
method         = "bilinear"  # regridding method

# Make the output name
outpath        = "/stormtrack/data4/glliu/01_Data/CESM2_PiControl/proc/"
outname        = addstrtoext(ncsearch,"_regrid_%s" % (method),adjust=-1)

# Indicate Reference
refpath        = "/stormtrack/data4/glliu/01_Data/CESM2_PiControl/proc/"
refnc          = refpath + "FCM_1600_2000_TS_clim.nc"

print("File with be regridded to: %s" % outname)
print("Reference grid is        : %s" % (refnc))
#%% Get list of Netcdfs

# Get list of NetCDF files
nclist = glob.glob(ncsearch)
nclist.sort()
print("Found %i files" % (len(nclist)))
print(nclist)

#%% Load data

st = time.time()
if len(nclist) == 1:
    dstarg = xr.open_dataset(nclist[0])
else:
    dstarg = xr.open_mfdataset(nclist,concat_dim=concatdim)
dstarg = dstarg.load()

# Load reference dataset
dsref = xr.open_dataset(refnc).load()
print("Loaded data in %.2fs" % (time.time()-st))


# ----------
#%% Make the Grid fron the reference file
# ----------

lat     = dsref.lat
lon     = dsref.lon
xx,yy   = np.meshgrid(lon,lat)
newgrid = xr.Dataset({'lat':lat,'lon':lon})

# ----------
#%% Set some additional settings based on user input
# ----------


# Rename Latitude/Longitude to prepare for regridding
ds          = dstarg.rename({"TLONG": "lon", "TLAT": "lat"})
ds          = ds.squeeze()

#%% Set up Regridder

oldgrid     = ds
regridder   = xe.Regridder(oldgrid,newgrid,method,periodic=True)

if type(ds) == xr.core.dataset.Dataset:
    ds = ds[vname]

# Regrid
st     = time.time()
daproc = regridder(ds) # Need to input dataarray
print("Regridded in %.2fs" % (time.time()-st))

#%% Save Output

edict = make_encoding_dict(daproc)
daproc.to_netcdf(outname,
                  encoding=edict)
print("Saved output to %s" % outname)



# # Indicate the Input Variable
# varname        = "VVEL"
# ncname         = "CESM1_HTR_VVEL_NATL_scycle.nc"
# ncpath         = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/01_Data/CESM1_LE/proc/NATL/"

# # Indicate netcdf information for reference grid
# ncref       = "CESM1LE_SST_NAtl_19200101_20050101_bilinear_stdev.nc"
# ncref_path  = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/01_Data/proc/CESM1/NATL_proc/"


# # Output Information
# outname     = "CESM1_HTR_VVEL_NATL_scycle_regrid_bilinear.nc"
# outpath     = ncpath

# method      = "bilinear"  # regridding method

# #%% Open Datasets

# dstarg = xr.open_dataset(ncpath+ncname).load()
# dsref  = xr.open_dataset(ncref_path + ncref).load()

# # ---------------------------------
# #%% Load Lat/Lon Universal Variable
# # ---------------------------------

# # Set up reference lat/lon
# lat     = dsref.lat
# lon     = dsref.lon
# xx,yy   = np.meshgrid(lon,lat)

# newgrid = xr.Dataset({'lat':lat,'lon':lon})

# # ----------
# #%% Set some additional settings based on user input
# # ----------

# # Rename Latitude/Longitude to prepare for regridding
# ds          = dstarg.rename({"TLONG": "lon", "TLAT": "lat"})
# oldgrid     = ds.isel(ens=0,z_t=0) # Make it 3D


# # Set up Regridder
# regridder   = xe.Regridder(oldgrid,newgrid,method,periodic=True)

# # Regrid
# st = time.time()
# daproc = regridder(ds) # Need to input dataarray
# print("Regridded in %.2fs" % (time.time()-st))


# savename = outpath + outname
# daproc.to_netcdf(savename,
#                  encoding={varname: {'zlib': True}})

# print("Saved output to %s" % savename)
