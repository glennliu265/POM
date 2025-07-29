#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate the climatology and monthly stdev for a variable in CESM2 PiControl

Created on Mon May  5 16:51:58 2025

@author: gliu
"""


import xarray as xr
import numpy as np
import matplotlib as mpl

import matplotlib.pyplot as plt
import sys
import glob
import os

import tqdm
import time

#%% Import Stuff

amvpath = "/home/glliu/00_Scripts/01_Projects/00_Commons/" # amv module
scmpath = "/home/glliu/00_Scripts/01_Projects/01_AMV/02_stochmod/stochmod/model/" # scm module

sys.path.append(amvpath)
sys.path.append(scmpath)

from amv import proc,viz
import scm
import amv.loaders as dl

#%% User Edits

outpath = "/stormtrack/data4/glliu/01_Data/CESM2_PiControl/proc/"

# HMXL, PenOM with 3-level restoring (NOTE just grabs all available values)
vname   = "HMXL"
expname = "POM03"
outname = "%s%s_%s_0101_0499_clim.nc" % (outpath,expname,vname)

# HMXL, PiControl CESM2
vname   = "HMXL"
expname = "FCM"
outname = "%s%s_%s_0100_2000_clim.nc" % (outpath,expname,vname)

print("Output will be saved to %s" % outname)





# Indicate netcdf information for reference grid
ncref       = "CESM1LE_SST_NAtl_19200101_20050101_bilinear_stdev.nc"
ncref_path  = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/01_Data/proc/CESM1/NATL_proc/"

#%% Indicate Data Path based on input

if expname == "POM03":
    dpath   = "/stormtrack/data4/glliu/01_Data/CESM2_PiControl/POM/"
    
elif expname == "FCM":
    dpath = "/stormtrack/data4/glliu/01_Data/CESM2_PiControl/FCM/ocn/%s/" % vname
elif expname == "SOM":    
    dpath = "/stormtrack/data4/glliu/01_Data/CESM2_PiControl/SOM/"
    
#%% Get list of netcdfs

ncsearch = dpath + "*%s.*" % vname
nclist   = glob.glob(ncsearch)
nclist.sort()
print("Found %i files" % (len(nclist)))
print(nclist)

if expname == "POM" or expname == "FCM":
    print("Dropping first 100 years (%s)" % (nclist[0]))
    nclist = nclist[1:]

#%% Open netcdf
st = time.time()
keepvars = [vname,"TLONG","TLAT",'time']
ds = xr.open_mfdataset(nclist,concat_dim='time',combine='nested')
ds = proc.ds_dropvars(ds,keepvars=keepvars)
ds = ds.load()
print("Loaded data in %.2fs" % (time.time()-st))

#%% Calculate climatology

st          = time.time()
ds          = proc.fix_febstart(ds)
dsclim      = ds.groupby('time.month').mean('time')
edict       = proc.make_encoding_dict(dsclim)
dsclim.to_netcdf(outname,encoding=edict)
print("Saved in %.2fs" % (time.time()-st))

#%% Calculate Stdev

dsstd       = ds.groupby('time.month').std('time')
outname_std = outname.replace("clim","std")
print(outname_std)
dsstd.to_netcdf(outname_std,encoding=edict)
print("Saved in %.2fs" % (time.time()-st))


#%% Not loopable, but isolate and save forcing file HBLT
dpath2 = "/stormtrack/data4/glliu/01_Data/CESM2_PiControl/SOM/forcing/"
ncname = dpath2 + "pop_frc.b.e21.B1850.f09_g17.CMIP6-piControl.001_branch2.012120.nc"
ds2    = xr.open_dataset(ncname).load()

keepvars = ["time",'xc','yc','hblt']
ds2    = proc.ds_dropvars(ds2,keepvars)
ds2    = ds2.rename(dict(time='month',xc="TLONG",yc="TLAT"))
outpath = "/stormtrack/data4/glliu/01_Data/CESM2_PiControl/proc/"
outname = "%sSOM_HBLT_clim.nc" % (outpath)
ds2.to_netcdf(outname)


