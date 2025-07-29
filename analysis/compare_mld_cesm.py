#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compare Mixed-Layer Depth Climatologies in CESM2 PiControl Runs

Created on Mon May  5 17:53:25 2025

@author: gliu

"""


import xarray as xr
import numpy as np
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import scipy as sp

import matplotlib.patheffects as PathEffects

import matplotlib.pyplot as plt
import sys
import glob
import os

import tqdm
import time
import cmcrameri as cm

#%% 

# local device (currently set to run on Astraeus, customize later)
amvpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/03_Scripts/" # amv module
scmpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/03_Scripts/stochmod/model/"

sys.path.append(amvpath)
sys.path.append(scmpath)

from amv import proc,viz
import scm
import amv.loaders as dl
import cvd_utils as cvd

#%% Load the files

dpath = "/Users/gliu/Downloads/02_Research/01_Projects/06_POM/01_Data/proc/"

ncfcm = "FCM_HMXL_0100_2000_clim_regrid_bilinear.nc"
ncpom = "POM03_HMXL_0101_0499_clim_regrid_bilinear.nc"
ncsom = "SOM_HBLT_clim_regrid_bilinear.nc"

ncs   = [ncsom,ncpom,ncfcm]
names = ["SOM","POM","FCM"]

#%% Load the climatologies

dsall = []
for ii in range(3):
    ds = xr.open_dataset(dpath + ncs[ii]).load()
    if ii == 0:
        ds = ds.hblt
        
    else:
        ds = ds.HMXL
        ds = ds/100
    dsall.append(ds)
    
## Convert to meters
#dsall = [ds/100 for ds in dsall]

# Replace time
#dsall[0].month = np.arange(1,13,1)
#%% First, let's compare the mean difference relative to FCM

somdiff = dsall[0].mean('month') - dsall[-1].mean('month')
pomdiff = dsall[1].mean('month') - dsall[-1].mean('month')


#%%

ii = 1

fsz_axis   = 16
fsz_title  = 24
cints = np.arange(-100,105,5)
cproj = ccrs.PlateCarree()

fig,ax   = init_globalmap()

if ii == 0:
    comparename = "SOM - FCM"
    cnout    = "SOMvFCM"
    plotvar  = somdiff
else:
    comparename = "POM - FCM"
    cnout = "POMvFCM"
    plotvar  = pomdiff
pcm      = ax.contourf(plotvar.lon,plotvar.lat,plotvar,
                      levels=cints,cmap='cmo.balance',
                      transform=cproj,extend='both')

cl       = ax.contour(plotvar.lon,plotvar.lat,plotvar,
                      levels=cints,colors="k",
                      transform=cproj,linewidths=0.55,alpha=0.5)
ax.clabel(cl,levels=cints[::4])


cb       = viz.hcbar(pcm,ax=ax)
cb.set_label("MLD Difference [m]",fontsize=fsz_axis)
#ax.add_feature(cfeature.LAND,zorder=1,color="k")

ax.set_title(comparename,fontsize=fsz_title)


#%%

def init_globalmap(nrow=1,ncol=1,figsize=(12,8)):
    proj            = ccrs.Robinson()
    bbox            = [-180,180,-90,90]
    fig,ax          = plt.subplots(nrow,ncol,subplot_kw={'projection':proj},figsize=figsize,constrained_layout=True)
    ax.coastlines()
    ax.gridlines(ls ='dotted',draw_labels=True)
    return fig,ax




