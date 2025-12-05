#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Visualize Re-emergence Index
- Uses regridded REI from [regrid_REI_MCOM]
- Computed from [calculate_REI_general] in re-emergence folder
- Uses ACF computed by [pointwise_crosscorrelation] (also in re-emergence folder)


First
- Regrid (bilinear) to CAM5 grid...
- 

Created on Mon Nov 17 13:31:58 2025

@author: gliu

"""

import time
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import xarray as xr
import sys
import tqdm
import glob 
import scipy as sp
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.gridspec as gridspec

from scipy.io import loadmat
import matplotlib as mpl

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


#%%

acfpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/01_Data/proc/"
datpath = "/Users/gliu/Downloads/02_Research/01_Projects/06_POM/01_Data/proc/acfs/"
outpath = datpath
figpath = "/Users/gliu/Downloads/02_Research/01_Projects/06_POM/02_Figures/20251119/"
proc.makedir(figpath)

ncnames = [
    "SOM_006-0360_REI_regridded.nc",
    "MCOM_0100-0599_REI_regridded.nc",
    "FOM_1100-1599_REI_regridded.nc",
    ]

acfnames = [
    "SOM_006-0360_lag00to60_ALL_ensALL.nc",
    "MCOM_0100-0599_lag00to60_ALL_ensALL.nc",
    "FOM_1100-1599_lag00to60_ALL_ensALL.nc"
    ]

expnames=["SOM","MCOM","FOM"]


#%% Load Data

st = time.time()
# Load REI
ds_all = []
for nc in ncnames:
    ds = xr.open_dataset(datpath+nc).load()
    ds_all.append(ds)

# Load ACF
ds_acfs = []
for nc in acfnames:
    ds = xr.open_dataset(acfpath+nc).load()
    ds_acfs.append(ds)

# Squeeze TLAT/TLON and rename
renamedict = dict(lat='nlat',lon='nlon')
ds_acfs = [ds.squeeze().rename(renamedict) for ds in ds_acfs]

print("Loaded datasets in %.2fs" % (time.time()-st))    

#%% Make Some Visualizations, based on parameters from the SMIO paper

bboxplot    = [-80, 0, 20, 65]
mpl.rcParams['font.family'] = 'Avenir'
mons3       = proc.get_monstr(nletters=3)

fsz_tick    = 18
fsz_axis    = 20
fsz_title   = 16

rhocrit     = proc.ttest_rho(0.05, 2, 86)
proj        = ccrs.PlateCarree()

def init_globalmap(nrow=1,ncol=1,figsize=(12,8)):
    proj            = ccrs.Robinson()
    bbox            = [-180,180,-90,90]
    fig,ax          = plt.subplots(nrow,ncol,subplot_kw={'projection':proj},figsize=figsize,constrained_layout=True)
    ax.coastlines()
    ax.gridlines(ls ='dotted',draw_labels=True)
    ax.add_feature(cfeature.LAND,facecolor="k",zorder=2)
    
    return fig,ax

import cartopy

def init_globalmap(nrow=1,ncol=1,figsize=(12,8),centlon=-180):
    proj            = ccrs.Robinson(central_longitude=centlon)
    bbox            = [-180,180,-90,90]
    fig,ax          = plt.subplots(nrow,ncol,subplot_kw={'projection':proj},figsize=figsize,constrained_layout=True)
    
    multiax = True
    if (type(ax) == mpl.axes._axes.Axes) or (type(ax) == cartopy.mpl.geoaxes.GeoAxes):
        ax = [ax,]
        multiax = False
    
    if type(ax) == tuple or (ncol+nrow > 2):
        ax = ax.flatten()
    for a in ax:
        a.coastlines(zorder=10,lw=0.75,)#transform=proj)
        a.gridlines(ls ='dotted',draw_labels=True)
        a.add_feature(cfeature.LAND,facecolor="k",zorder=2)
    if multiax is False:
        ax = ax[0]
    else:
        ax = ax.reshape(nrow,ncol).squeeze()
    return fig,ax

#%%

ex  = 2 
iyr = 0
im  = 0



plt.rcParams.update({'font.size': 16})
cints_rei = np.arange(0,0.60,0.05)

vmax     = 0.5

for ex in range(3):
    for im in range(12):
        
        fig,ax   = init_globalmap()
        
        plotvar  = ds_all[ex].rei.isel(yr=iyr,mons=im)
        
        #pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,vmin=0,vmax=vmax,cmap='cmo.dense')
        
        pcm     = ax.contourf(plotvar.lon,plotvar.lat,plotvar,transform=proj,levels=cints_rei,cmap='cmo.dense',)#extend='both')
        cl      = ax.contour(plotvar.lon,plotvar.lat,plotvar,transform=proj,levels=cints_rei,colors='dimgray',linewidths=0.25)
        clbl    = ax.clabel(cl,levels=cints_rei[::2],fontsize=8)
        
        viz.add_fontborder(clbl,w=3)
        
        
        #ax.tick_params(labelsize=22)
        cb      = fig.colorbar(pcm,ax=ax,fraction=0.025,pad=0.01)
        ax.set_title("%s Year-%i Re-emergence Index (%s)" % (mons3[im],iyr+1,expnames[ex]))
        figname = "%sREI_%s_Y%i_mon%02i.png" % (figpath,expnames[ex],iyr+1,im+1)
        plt.savefig(figname,dpi=150,bbox_inches='tight')
        
        
#%% Make a combined plot of the JFM and JAS

iyr         = 0
#fsz_title   = 26

selmons     = [[0,1,2],[6,7,8]]


stitles     = ["Boreal Winter\n(JFM)","Austral Winter\n(JAS)"]

fig,axs     = init_globalmap(nrow=2,ncol=3,centlon=0,figsize=(22,8.5))

ii = 0
for ss in range(2):
    
    for ex in range(3):
        
        ax = axs[ss,ex]
        
        if ex == 0:
            viz.add_ylabel(stitles[ss],ax=ax,fontsize=fsz_axis)
        if ss == 0:
            ax.set_title(expnames[ex])
            
        
        # Do the plotting
        plotvar = ds_all[ex].rei.isel(yr=iyr,mons=selmons[ss]).mean('mons')
        
        pcm     = ax.contourf(plotvar.lon,plotvar.lat,plotvar,transform=proj,levels=cints_rei,cmap='cmo.dense',extend='both')#extend='both')
        #cl      = ax.contour(plotvar.lon,plotvar.lat,plotvar,transform=proj,levels=cints_rei,colors='dimgray',linewidths=0.25)
        viz.label_sp(ii,ax=ax,fontsize=fsz_axis,y=1.1)
        ii+=1
            
cb = viz.hcbar(pcm,ax=axs.flatten(),fraction=0.035,pad=0.01)
cb.set_label("Re-emergence Index [Correlation]")

figname = "%sREI_Yr0%i_Comparison_Plot_JFM_JAS.png" % (figpath,iyr+1)
plt.savefig(figname,dpi=150,bbox_inches='tight')

#for ss in range(2):
    
    
#%% Plot the *MAX* Re-emergence


vmax     = 0.5

for ex in range(3):
    
    fig,axs   = init_globalmap(1,1,centlon=0)
    
    # Plot the amplitude
    ax        = axs#[0]
    plotvar  = ds_all[ex].rei.isel(yr=iyr).max('mons')#,mons=im)
    pcm     = ax.contourf(plotvar.lon,plotvar.lat,plotvar,transform=proj,levels=cints_rei,cmap='cmo.dense',extend='both')
    cl      = ax.contour(plotvar.lon,plotvar.lat,plotvar,transform=proj,levels=cints_rei,colors='dimgray',linewidths=0.25)
    clbl    = ax.clabel(cl,levels=cints_rei[::2],fontsize=8)
    viz.add_fontborder(clbl,w=3)
    cb      = fig.colorbar(pcm,ax=ax,fraction=0.025,pad=0.01)
    ax.set_title("Max Year-%i Re-emergence Index (%s)" % (iyr+1,expnames[ex]))
    
    # #
    # ax = axs[1]
    # plotvar  = np.nanargmax(ds_all[ex].rei.isel(yr=iyr),0) + 1#,mons=im)
    # pcm     = ax.contourf(plotvar.lon,plotvar.lat,plotvar,transform=proj,levels=np.arange(1,13,1),cmap='twilight',)#extend='both')
    # cl      = ax.contour(plotvar.lon,plotvar.lat,plotvar,transform=proj,levels=np.arange(1,13,1),colors='dimgray',linewidths=0.25)
    # clbl    = ax.clabel(cl,levels=cints_rei[::2],fontsize=8)
    # viz.add_fontborder(clbl,w=3)
    # cb      = fig.colorbar(pcm,ax=ax,fraction=0.025,pad=0.01)
    # ax.set_title("Month of Maximum")
    
    # Plot the month of max
    
    #pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,vmin=0,vmax=vmax,cmap='cmo.dense')
    
    
    
    
    
    
    #ax.tick_params(labelsize=22)
    
    figname = "%sREI_%s_Y%i_MAX.png" % (figpath,expnames[ex],iyr+1,)
    plt.savefig(figname,dpi=150,bbox_inches='tight')
    

#%% Check ACF characteristics at a point
lonf           = 330
latf           = 50
locfn,loctitle = proc.make_locstring(lonf,latf)
# for ds in ds_acfs:
    
#     ds['TLAT'] = ds['TLAT'].squeeze()
#     ds['TLONG'] = ds['TLONG'].squeeze()
    

dspts   = [proc.find_tlatlon(ds,lonf,latf) for ds in ds_acfs]

#%%

fig,ax   = init_globalmap(1,1,centlon=0)
ii       = 1

plotvar = ds_all[ii].rei.isel(yr=iyr).max('mons')#,mons=im)
pcm     = ax.contourf(plotvar.lon,plotvar.lat,plotvar,transform=proj,levels=cints_rei,cmap='cmo.dense',extend='both')
cl      = ax.contour(plotvar.lon,plotvar.lat,plotvar,transform=proj,levels=cints_rei,colors='dimgray',linewidths=0.25)

# Regions with Strong Re-emergence in FOM
# plotpoints = [[-45,42],[154,40],[-30,4],[85,-2.5],[175,7.5],[5,-38],[-108,15],[130,-55]]
# pointnames = ["Gulf Stream/North Atlantic Current","Kuroshio",
#               "Tropical Atlantic","Eastern Indian Ocean","North Equatorial Pacific",
#               "Agulhas Downstream","Gap Winds","Southern Ocean"]


# Southern Ocean Checkpoints
plotpoints = [[130,-55],[111,-45],[90,-48],[32,-45],[-56,-53],[-80,-58],[174,-50],[143,-45],[-170,-44]]
pointnames = ["130E","111E","90E","Agulhas","E Malvinas","W Drake Passage","S New Zealand","S Australia","S. Pacific"]

for pp in plotpoints:
    ax.plot(pp[0],pp[1],transform=proj,markersize=25,marker="x")



#%% Selecting a Point


# fig,axs   = init_globalmap(3,1,centlon=0,figsize=(12,16))

# for ii in range(3):
    
#     if ii != 2:
#         continue
    
#     ax      = axs[ii]
#     plotvar = ds_all[ex].rei.isel(yr=iyr).max('mons')#,mons=im)
#     pcm     = ax.contourf(plotvar.lon,plotvar.lat,plotvar,transform=proj,levels=cints_rei,cmap='cmo.dense',extend='both')
#     cl      = ax.contour(plotvar.lon,plotvar.lat,plotvar,transform=proj,levels=cints_rei,colors='dimgray',linewidths=0.25)

# plotpoints = [[-50,45]]

# for pp in plotpoints:
#     ax.plot(pp[0],pp[1],transform=proj,markersize=25,marker="x")

#%% Plot points

kmonth  = 1
checkex = 2
lags    = dspts[0].isel(mons=kmonth).lags
xtks    = lags[::3]
    
for p,pp in enumerate(plotpoints):
    
    lonf,latf = pp
    if lonf < 0:
        lonf = lonf + 360
    locfn,loctitle = proc.make_locstring(lonf,latf)
    dspts   = [proc.find_tlatlon(ds,lonf,latf) for ds in ds_acfs]
    
    dscheck = ((dspts[checkex].acf)**2).sum('lags')
    kmonth  = dscheck.argmax().data.item()
    
    
    # Make the plot
    
    fig,ax = plt.subplots(1,1,constrained_layout=True,figsize=(8,4.5))
    ax,_=viz.init_acplot(kmonth,xtks,lags,ax=ax,title="")

    for ii in range(3):
        
        plotvar = dspts[ii].isel(mons=kmonth).acf
        ax.plot(plotvar.lags,plotvar,label=expnames[ii],c=expcols[ii],lw=2.5)
    ax.set_title("%s Autocorrelation Function @ %s: %s" % (mons3[kmonth],pointnames[p],loctitle))

    ax.set_xlabel("Lag from %s [months]" % mons3[kmonth])
    ax.set_ylabel("Correlation with %s Anomalies" % (mons3[kmonth]))
    ax.legend()
      
    figname = "%sACF_Example_%s.png" % (figpath,locfn)
    plt.savefig(figname,dpi=150,bbox_inches='tight')  

    
    




#%%


expcols = ["hotpink","forestgreen","midnightblue"]

kmonth  = 1
lags    = dspts[0].isel(mons=kmonth).lags
xtks    = lags[::3]

fig,ax = plt.subplots(1,1,constrained_layout=True,figsize=(8,4.5))
ax,_=viz.init_acplot(kmonth,xtks,lags,ax=ax,title="")

for ii in range(3):
    
    plotvar = dspts[ii].isel(mons=kmonth).acf
    ax.plot(plotvar.lags,plotvar,label=expnames[ii],c=expcols[ii],lw=2.5)
ax.set_title("%s Autocorrelation Function @ %s" % (mons3[kmonth],loctitle))

ax.set_xlabel("Lag from %s [months]" % mons3[kmonth])
ax.set_ylabel("Correlation with %s Anomalies" % (mons3[kmonth]))
ax.legend()
  
figname = "%sACF_Example_%s.png" % (figpath,locfn)
plt.savefig(figname,dpi=150,bbox_inches='tight')  


#%% 

