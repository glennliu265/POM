#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Exploratory analyis of Pencil Ocean Model Output

Created on Thu Apr 24 10:15:29 2025

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

#%% Plotting Inputs

# Set Plotting Options
darkmode = False
if darkmode:
    dfcol = "w"
    bgcol = np.array([15,15,15])/256
    sp_alpha = 0.05
    transparent = True
    plt.style.use('dark_background')
    mpl.rcParams['font.family']     = 'JetBrainsMono'
else:
    dfcol = "k"
    bgcol = "w"
    sp_alpha = 0.75
    transparent = False
    plt.style.use('default')

bboxplot    = [-80, 0, 20, 65]
mpl.rcParams['font.family'] = 'Avenir'
mons3       = proc.get_monstr(nletters=3)

fsz_tick    = 18
fsz_axis    = 20
fsz_title   = 16

rhocrit     = proc.ttest_rho(0.05, 2, 86)
proj        = ccrs.PlateCarree()


#%% User Edits

figpath = "/Users/gliu/Downloads/02_Research/01_Projects/06_POM/02_Figures/20250430/"
datpath = "/Users/gliu/Downloads/02_Research/01_Projects/06_POM/01_Data/proc/"
proc.makedir(figpath)


expnames = [
    "SOM_0001_0360",
    "POM03_0101-499",
    "POM05_0101-499",
    "FCM_1600_2000",
    ]

expnames_short = [
    "SOM",
    "MCOM",
    "MCOM",
    "FCM"
    ]

expnames_long = [
    "Slab Ocean (Years 1 to 360)",
    "Multi-Column Ocean, 3-level U/V Restoring (Years 101 to 499)",
    "Multi-Column Ocean, 5-level U/V Restoring (Years 101 to 499)",
    "Fully Coupled Model (Years 1600 to 2000)"
    ]


expcols = [
    "tomato",
    "deepskyblue",
    "mediumblue",
    "k",    
    ]

expmks = [
    "d",
    "^",
    "P",
    ".",
    ]

nexps   = len(expnames)

# Load Ice Mask
maskname = "cesm2_pic_limask_0.3p_0.05p_0200to2000.nc"
maskpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/01_Data/proc/model_input/masks/"
dsmask   = xr.open_dataset(maskpath + maskname).load()
dsmask180 = proc.lon360to180_xr(dsmask)

plotmask = xr.where(np.isnan(dsmask180.mask),0,1)
#%% Load the inputs

dsanoms = []
dsclims = []
for ex in tqdm.tqdm(range(nexps)):
    
    expname = expnames[ex]
    
    # Load Anomalies
    ncname  = "%s%s_TS_anom.nc" % (datpath,expname)
    dsanom  = xr.open_dataset(ncname).load()
    dsanoms.append(dsanom)
    
    # Load Climatology
    ncname  = "%s%s_TS_clim.nc" % (datpath,expname)
    dsclim  = xr.open_dataset(ncname).load()
    dsclims.append(dsclim)



dsanoms = [proc.format_ds(ds) for ds in dsanoms]
dsclims = [proc.lon360to180_xr(ds) for ds in dsclims]
stdevs  = [ds.std('time') for ds in dsanoms]


dsmasked = [ds.TS * dsmask180.mask for ds in dsanoms]
#%% Make a function

def init_globalmap(nrow=1,ncol=1,figsize=(12,8)):
    proj            = ccrs.Robinson()
    bbox            = [-180,180,-90,90]
    fig,ax          = plt.subplots(nrow,ncol,subplot_kw={'projection':proj},figsize=figsize,constrained_layout=True)
    ax.coastlines()
    ax.gridlines(ls ='dotted',draw_labels=True)
    return fig,ax

#%% Compare the stdev and mean

fsz_axis   = 16
fsz_title  = 24

ex         = 0
cproj      = ccrs.PlateCarree()
cints      = np.arange(0,1.55,0.05)
cints_mean = np.arange(270,320,2)

for ex in range(nexps):
    fig,ax   = init_globalmap()
    
    plotvar  = stdevs[ex].TS
    pcm      = ax.contourf(plotvar.lon,plotvar.lat,plotvar,
                          levels=cints,cmap='cmo.thermal',
                          transform=cproj,extend='both')
    
    plotvar  = (dsclims[ex].mean('month') * dsmask180.mask).TS
    cl       = ax.contour(plotvar.lon,plotvar.lat,plotvar,
                          levels=cints_mean,colors="w",
                          transform=cproj,linewidths=0.55,alpha=0.5)
    ax.clabel(cl)
    
    # Plot Ice
    ax.contour(plotmask.lon,plotmask.lat,plotmask,colors='cyan',
               transform=cproj,linewidths=0.66,zorder=1,linestyles='solid')
    
    cb       = viz.hcbar(pcm,ax=ax)
    cb.set_label("TS 1$\sigma$ [$\degree$C]",fontsize=fsz_axis)
    ax.add_feature(cfeature.LAND,zorder=1,color="k")
    
    ax.set_title(expnames_long[ex],fontsize=fsz_title)
    
    
    figname = "%sTS_Stdev_Mean_Global_%s.png" % (figpath,expnames[ex])
    plt.savefig(figname,dpi=150,bbox_inches='tight')



#%% How does log ratio of stdev compare to the Fully Coupled Model?

cints           = np.log(np.array([0.25, 0.5, 1, 1.25,1.5,2]))
clabs           = ["0.25x","0.5x","1x", "1.25x", "1.5x", "2x", ]
vmax            = 1.5
fsz_ticks       = 14
refvar = stdevs[-1].TS.data

lon    = dsanoms[0].lon
lat    = dsanoms[0].lat

levels = np.arange(0,2,0.2)

for ex in range(3):
    
    fig,ax   = init_globalmap()
    
    plotvar  = np.log((stdevs[ex].TS.data / refvar))
    pcm      = ax.pcolormesh(lon,lat,plotvar,
                          cmap='cmo.balance',vmin=-vmax,vmax=vmax,
                          transform=cproj)
    
    # Add Ratio Plots
    cl = ax.contour(lon,lat,plotvar * dsmask180.mask,transform=cproj,
                            levels=cints,colors="dimgray",linewidths=0.75)
    fmt= {}
    for l, s in zip(cl.levels, clabs):
        fmt[l] = s
    cl = ax.clabel(cl,fmt=fmt,fontsize=fsz_ticks)
    viz.add_fontborder(cl)
    
    # Plot Ice
    ax.contour(plotmask.lon,plotmask.lat,plotmask,colors='cyan',
               transform=cproj,linewidths=0.66,zorder=1,linestyles='solid')
    
    cb       = viz.hcbar(pcm,ax=ax)
    cb.set_label(r"Log $\frac{\sigma(%s)}{\sigma(FCM)}$" % (expnames_short[ex]),fontsize=fsz_axis)
    ax.add_feature(cfeature.LAND,zorder=1,color="k")
    
    ax.set_title(expnames_long[ex],fontsize=fsz_title)
    
    
    figname = "%sTS_Ratio_withFCM_Global_%s.png" % (figpath,expnames[ex])
    plt.savefig(figname,dpi=150,bbox_inches='tight')

#%% 2025.07.29: Remake plots for NSF Project Report
# Plot the more complex mode in the numerator

cints           = np.log(np.array([0.25, 0.5, 0.75, 1, 1.25,1.5,1.75,2]))
clabs           = ["0.25x","0.5x","0.75x","1x", "1.25x", "1.5x", "1.75","2x", ]
vmax            = 1.5
fsz_ticks       = 14

# Indicate NUmerator and denominator
inumer = 3#3#1
idenom = 0#1#0


numervar = stdevs[inumer].TS.data
denomvar = stdevs[idenom].TS.data

comparename_long = r"Log $\frac{\sigma(%s)}{\sigma(%s)}$" % (expnames_short[inumer],
                                                             expnames_short[idenom])
comparename_short = "%sto%s" % (expnames_short[inumer],
                                                             expnames_short[idenom])

fig,ax   = init_globalmap()

plotvar  = np.log((numervar / denomvar))
pcm      = ax.pcolormesh(lon,lat,plotvar,
                      cmap='cmo.balance',vmin=-vmax,vmax=vmax,
                      transform=cproj)

# Add Ratio Plots
cl = ax.contour(lon,lat,plotvar * dsmask180.mask,transform=cproj,
                        levels=cints,colors="dimgray",linewidths=0.75)
fmt= {}
for l, s in zip(cl.levels, clabs):
    fmt[l] = s
cl = ax.clabel(cl,fmt=fmt,fontsize=fsz_ticks)
viz.add_fontborder(cl)

# Plot Ice
ax.contour(plotmask.lon,plotmask.lat,plotmask,colors='cyan',
           transform=cproj,linewidths=0.66,zorder=1,linestyles='solid')

cb       = viz.hcbar(pcm,ax=ax)
cb.set_label(comparename_long,fontsize=fsz_axis)
ax.add_feature(cfeature.LAND,zorder=1,color="k")

#ax.set_title(expnames_long[ex],fontsize=fsz_title)


figname = "%sTS_Ratio_Global_%s.png" % (figpath,comparename_short)
plt.savefig(figname,dpi=150,bbox_inches='tight')

#%% Calculate R1



def calc_lagcorr(ts,lag):
    if np.any(np.isnan(ts)):
        return np.nan
    ntime = len(ts)
    r1 = np.corrcoef(ts[:(ntime-lag)],ts[lag:])[0,1]
    return r1
    

lagr1 = lambda x : calc_lagcorr(x,1)


r1s = []
for ex in range(nexps):
    #Runs in approx 160 sec
    st         = time.time()
    
    r1out = xr.apply_ufunc(
        lagr1,
        dsmasked[ex],
        input_core_dims=[['time']],
        output_core_dims=[[]],
        vectorize=True,
    )
    
    
    #r1out      = dsmasked[ex].reduce(lagr1,dim='time')
    print("Function applied in in %.2fs" % (time.time()-st))
    
    # Apply looping through basemonth, lon, lat. ('lon', 'lat', 'mon', 'rem_year')
    
    r1s.append(r1out)
    
#%% Calc Lag 6


lagr6 = lambda x : calc_lagcorr(x,6)


r6s = []
for ex in range(nexps):
    #Runs in approx 160 sec
    st         = time.time()
    
    r1out = xr.apply_ufunc(
        lagr6,
        dsmasked[ex],
        input_core_dims=[['time']],
        output_core_dims=[[]],
        vectorize=True,
    )
    
    
    #r1out      = dsmasked[ex].reduce(lagr1,dim='time')
    print("Function applied in in %.2fs" % (time.time()-st))
    
    # Apply looping through basemonth, lon, lat. ('lon', 'lat', 'mon', 'rem_year')
    
    r6s.append(r1out)

#%% Calc Lag 12


lagr12 = lambda x : calc_lagcorr(x,12)


r12s = []
for ex in range(nexps):
    #Runs in approx 160 sec
    st         = time.time()
    
    r1out = xr.apply_ufunc(
        lagr12,
        dsmasked[ex],
        input_core_dims=[['time']],
        output_core_dims=[[]],
        vectorize=True,
    )
    
    
    #r1out      = dsmasked[ex].reduce(lagr1,dim='time')
    print("Function applied in in %.2fs" % (time.time()-st))
    
    # Apply looping through basemonth, lon, lat. ('lon', 'lat', 'mon', 'rem_year')
    
    r12s.append(r1out)



#%% Plot R1/R6/R12 Globally

plotname = "R6"

fsz_axis   = 16
fsz_title  = 24

ex         = 0
cproj      = ccrs.PlateCarree()
cints      = np.arange(-1,1.05,0.05)


for ex in range(nexps):
    
    if plotname == "R1":
        plotvar  = r1s[ex]
        lag      = 1
    elif plotname == "R12":
        plotvar  = r12s[ex]
        lag       = 12
    elif plotname == "R6":
        plotvar = r6s[ex]
        lag     = 6
    
    
    fig,ax   = init_globalmap()
    
    
    pcm      = ax.contourf(plotvar.lon,plotvar.lat,plotvar,
                          levels=cints,cmap='cmo.balance',
                          transform=cproj,extend='both')
    
    
    cb       = viz.hcbar(pcm,ax=ax)
    cb.set_label("Lag %s Correlation" % lag,fontsize=fsz_axis)
    ax.add_feature(cfeature.LAND,zorder=1,color="k")
    
    ax.set_title(expnames_long[ex],fontsize=fsz_title)
    
    # Plot Ice
    ax.contour(plotmask.lon,plotmask.lat,plotmask,colors='cyan',
               transform=cproj,linewidths=0.66,zorder=1,linestyles='solid')
    
    
    figname = "%sTS_R%02i_Global_%s.png" % (figpath,lag,expnames[ex])
    plt.savefig(figname,dpi=150,bbox_inches='tight')


#%% Compute the AMV Index


amvbbox  = [-80,0,0,60]
sstin    = [ds.transpose('lon','lat','time') for ds in dsmasked]

sstin_unmasked = [ds.transpose('lon','lat','time').TS for ds in dsanoms]

amvids   = []
amvpats  = []
nasstis  = []

for ex in range(nexps):
    invar              = sstin_unmasked[ex]
    amvi,amvpat,nassti = proc.calc_AMVquick(invar.data,invar.lon.data,invar.lat.data,amvbbox,return_unsmoothed=True,
                                            mask=dsmask180.mask.data.squeeze().transpose(1,0))
    
    amvids.append(amvi)
    amvpats.append(amvpat)
    nasstis.append(nassti)


#%% Plot Global AMV Patterns

cints_amv = np.arange(-0.50,0.55,0.05)

for ex in range(nexps):
    
    if plotname == "R1":
        plotvar  = r1s[ex]
        lag      = 1
    elif plotname == "R12":
        plotvar  = r12s[ex]
        lag       = 12
    elif plotname == "R6":
        plotvar = r6s[ex]
        lag     = 6
    
    fig,ax   = init_globalmap()
    
    plotvar  = amvpats[ex].T
    pcm      = ax.contourf(lon,lat,plotvar,
                          levels=cints_amv,cmap='cmo.balance',
                          transform=cproj,extend='both')
    
    
    cb       = viz.hcbar(pcm,ax=ax)
    cb.set_label("AMV Pattern [$\degree$C per 1$\sigma$ AMV Index]",fontsize=fsz_axis)
    ax.add_feature(cfeature.LAND,zorder=1,color="k")
    
    ax.set_title(expnames_long[ex],fontsize=fsz_title)
    
    # Plot Ice
    ax.contour(plotmask.lon,plotmask.lat,plotmask,colors='cyan',
               transform=cproj,linewidths=0.66,zorder=1,linestyles='solid')
    
    
    figname = "%sTS_AMV_Pattern_Global_%s.png" % (figpath,expnames[ex])
    plt.savefig(figname,dpi=150,bbox_inches='tight')

#%% Plot AMV Indices

fig,axs = plt.subplots(4,1,constrained_layout=True,figsize=(12,8))
ylims = [-.8,.8]

for ex in range(nexps):
    
    
    ax = axs[ex]
    
    plotvar = nasstis[ex]
    trange  = np.arange(len(plotvar))
    ax.plot(trange,plotvar,c='gray',label="NASST Index")
    
    plotvar = amvids[ex]
    trange  = np.arange(len(plotvar))
    ax.plot(trange,plotvar,c='k',label="AMV Index")
    
    ax.axhline([0],ls='dotted',c='k',lw=0.55)
    
    ax.set_xlim([trange[0],trange[-1]])
    ax.set_ylim(ylims)
    
    if ex == 0:
        ax.legend(ncol=2,loc='lower center')
    ax.set_ylabel("%s SST [$\degree$C]" % expnames_short[ex])
    
    #ax.set_title("Sigma = %.2f" % (np.std(plotvar)))

figname = "%sTS_AMV_Index_Global.png" % (figpath)
plt.savefig(figname,dpi=150,bbox_inches='tight')
#%% Compare AMV Variance


rawstds = np.zeros((nexps))
lpstds  = np.zeros((nexps))

for ex in range(nexps):
    
    rawstds[ex] = np.std(nasstis[ex])
    lpstds[ex] = np.std(amvids[ex])

ratios = lpstds/rawstds


xlabs       = ["%s\n%.2f" % (expnames_short[ii],ratios[ii]*100)+"%" for ii in range(len(ratios))]

fig,ax = plt.subplots(1,1,constrained_layout=True,figsize=(6,6))

braw = ax.bar(np.arange(nexps),rawstds,color='gray')
blp  = ax.bar(np.arange(nexps),lpstds,color='k')

ax.bar_label(braw,fmt="%.04f",c='gray')
ax.bar_label(blp,fmt="%.04f",c='k')

ax.set_xticks(np.arange(nexps),labels=xlabs)
ax.set_ylabel("$\sigma$(SST Index) [$\degree$C]")

ax.set_ylim([0,1])

figname = "%sTS_AMV_NASST_Index_Stdev_Global.png" % (figpath)
plt.savefig(figname,dpi=150,bbox_inches='tight')

#%% Get AMV Variance for each pattern

aavg_amv = [proc.area_avg_cosweight(proc.sel_region_xr(amvpat * dsmask180.mask.T,amvbbox)) for amvpat in amvpats]



#%% Compute Monthly Spectra


sstreg_masked = [proc.sel_region_xr(ds * dsmask180.mask,amvbbox) for ds in sstin_unmasked]
nassti_mon    = [proc.area_avg_cosweight(ds) for ds in sstreg_masked]

nassti_mon_arr = [ds.data for ds in nassti_mon]

#%% Calculate and compute spectra

lags        = np.arange(61)
nsmooth     = 50
metrics_out = scm.compute_sm_metrics(nassti_mon_arr,nsmooth=nsmooth,lags=lags)

#%% Examine the ACF



n       = 400
conf    = 0.95
tails   = 2

kmonth = 1
xtks   = lags[::6]
fig,ax = plt.subplots(1,1,constrained_layout=True,figsize=(10,4))
ax,_   = viz.init_acplot(kmonth,xtks,lags,ax=ax)

for ex in range(nexps):
    plotvar = metrics_out['acfs'][kmonth][ex]
    cflags  = proc.calc_conflag(plotvar,conf,tails,n)
    ax.plot(lags,plotvar,marker=expmks[ex],c=expcols[ex],label=expnames_long[ex],lw=2.5,zorder=3)
    ax.fill_between(lags,cflags[:,0],cflags[:,1],color=expcols[ex],alpha=0.10,zorder=1)
    
ax.legend()
ax.set_title("%s Autocorrelation" % mons3[kmonth],fontsize=fsz_title)

savename = "%sNASST_Monthly_ACF_mon%02i.png" % (figpath,kmonth+1)
plt.savefig(savename,dpi=150,bbox_inches='tight')


#%% Plot the spectra


dtmon_fix       = 60*60*24*30
xper            = np.array([40,10,5,1,0.5])
xper_ticks      = 1 / (xper*12)

    
fig,ax          = plt.subplots(1,1,figsize=(8,4.5),constrained_layout=True)

for ex in range(nexps):
    plotspec        = metrics_out['specs'][ex] / dtmon_fix
    plotfreq        = metrics_out['freqs'][ex] * dtmon_fix
    CCs = metrics_out['CCs'][ex] / dtmon_fix

    ax.loglog(plotfreq,plotspec,lw=2.5,label=expnames_long[ex],c=expcols[ex])
    
    #ax.loglog(plotfreq,CCs[:,0],ls='dotted',lw=0.5,c=expcols[ii])
    #ax.loglog(plotfreq,CCs[:,1],ls='dashed',lw=0.9,c=expcols[ii])

ax.set_xlim([1/1000,0.5])
ax.axvline([1/(6)],label="",ls='dotted',c='gray')
ax.axvline([1/(12)],label="",ls='dotted',c='gray')
ax.axvline([1/(5*12)],label="",ls='dotted',c='gray')
ax.axvline([1/(10*12)],label="",ls='dotted',c='gray')
ax.axvline([1/(40*12)],label="",ls='dotted',c='gray')

ax.set_xlabel("Frequency (1/Month)",fontsize=14)
ax.set_ylabel("Power [$\degree C ^2 cycle \, per \, mon$]")

ax2 = ax.twiny()
ax2.set_xlim([1/1000,0.5])
ax2.set_xscale('log')
ax2.set_xticks(xper_ticks,labels=xper)
ax2.set_xlabel("Period (Years)",fontsize=14)


# Plot Confidence Interval (ERA5)
alpha           = 0.05
cloc_era        = [.5e-2,1e-1]
dof_era         = metrics_out['dofs'][-1]
cbnds_era       = proc.calc_confspec(alpha,dof_era)
proc.plot_conflog(cloc_era,cbnds_era,ax=ax,color='k',cflabel=r"95% Confidence") #+r" (dof= %.2f)" % dof_era)

ax.legend()


savename = "%sNASST_Monthly_Spectra_nsmooth%04i.png" % (figpath,nsmooth)
plt.savefig(savename,dpi=150,bbox_inches='tight')


#%% Linear Linear Plot


xlm_linear = [1/1000,0.1]
    
fig,ax          = plt.subplots(1,1,figsize=(8,4.5),constrained_layout=True)

for ex in range(nexps):
    plotspec        = metrics_out['specs'][ex] / dtmon_fix
    plotfreq        = metrics_out['freqs'][ex] * dtmon_fix
    CCs = metrics_out['CCs'][ex] / dtmon_fix

    ax.plot(plotfreq,plotspec,lw=2.5,label=expnames_long[ex],c=expcols[ex])
    


ax.set_xlim(xlm_linear)
ax.axvline([1/(6)],label="",ls='dotted',c='gray')
ax.axvline([1/(12)],label="",ls='dotted',c='gray')
ax.axvline([1/(5*12)],label="",ls='dotted',c='gray')
ax.axvline([1/(10*12)],label="",ls='dotted',c='gray')
ax.axvline([1/(40*12)],label="",ls='dotted',c='gray')

ax.set_xlabel("Frequency (1/Month)",fontsize=14)
ax.set_ylabel("Power [$\degree C ^2 cycle \, per \, mon$]")

ax2 = ax.twiny()
ax2.set_xlim(xlm_linear)
#ax2.set_xscale('log')
ax2.set_xticks(xper_ticks,labels=xper)
ax2.set_xlabel("Period (Years)",fontsize=14)
ax2.set_xlim(xlm_linear) # Set this later to not have xticks mess things up

# Plot Confidence Interval (ERA5)
alpha           = 0.05
cloc_era        = [.5e-2,1e-1]
dof_era         = metrics_out['dofs'][-1]
cbnds_era       = proc.calc_confspec(alpha,dof_era)
proc.plot_conflog(cloc_era,cbnds_era,ax=ax,color='k',cflabel=r"95% Confidence") #+r" (dof= %.2f)" % dof_era)

ax.legend()


savename = "%sNASST_Monthly_Spectra_Linear_nsmooth%04i.png" % (figpath,nsmooth)
plt.savefig(savename,dpi=150,bbox_inches='tight')



#%% Compute LP Filtered Timeseries to each map

cutoff      = 120 # In Months
order       = 6
lpfilter_ts = lambda x: proc.lp_butter(x,cutoff,order)



lpfilter_sst = []
for ex in range(nexps):
    
    #Runs in approx 160 sec
    st         = time.time()
    
    lpout = xr.apply_ufunc(
        lpfilter_ts,
        dsmasked[ex],
        input_core_dims=[['time']],
        output_core_dims=[['time']],
        vectorize=True,
    )
    
    print("Function applied in in %.2fs" % (time.time()-st))
    
    lpfilter_sst.append(lpout)

    

lpstd_sst = [ds.std('time') for ds in lpfilter_sst]
rawstd_sst = [ds.std('time') for ds in dsmasked]

#%% Check LP Filter Variance at a test point
ex   = -1
lonf = -30
latf = 50

ptraw = proc.selpt_ds(dsmasked[ex],lonf,latf)
ptlp  = proc.selpt_ds(lpfilter_sst[ex],lonf,latf)


fig,ax = plt.subplots(1,1,figsize=(12.5,3.5))
ax.plot(ptraw,label="Raw",c='gray')
ax.plot(ptlp,label="LP Filtered",c='red')


#%% Plot Raw (or LP Filtered) Variance

plot_raw    = False
ex          = -1
cproj       = ccrs.PlateCarree()
cints_raw    = np.arange(0,1.6,0.1) 

cints_lp    = np.arange(0,.90,0.05) 

for ex in range(nexps):
        
    fig,ax   = init_globalmap()
    
    plotvar  = lpstd_sst[ex] #rawstd_sst[ex]
    pcm      = ax.contourf(plotvar.lon,plotvar.lat,plotvar,
                          levels=cints_lp,cmap=cm.cm.lipari,
                          transform=cproj,extend='both')
    if plot_raw:
        plotvar  = rawstd_sst[ex]
        cl       = ax.contour(plotvar.lon,plotvar.lat,plotvar,
                              levels=cints_raw,colors="w",
                              transform=cproj,linewidths=0.55,alpha=0.5)
    else:
        cl       = ax.contour(plotvar.lon,plotvar.lat,plotvar,
                              levels=cints_lp,colors="w",
                              transform=cproj,linewidths=0.55,alpha=0.5)
        
    ax.clabel(cl)
    
    # Plot Ice
    ax.contour(plotmask.lon,plotmask.lat,plotmask,colors='cyan',
               transform=cproj,linewidths=0.66,zorder=1,linestyles='solid')
    
    cb       = viz.hcbar(pcm,ax=ax)
    cb.set_label("10-year LP Filter TS 1$\sigma$ [$\degree$C]",fontsize=fsz_axis)
    ax.add_feature(cfeature.LAND,zorder=1,color="k")
    ax.set_title(expnames_long[ex],fontsize=fsz_title)
    
    figname = "%sTS_LPfilter_Stdev_Global_%s.png" % (figpath,expnames[ex])
    plt.savefig(figname,dpi=150,bbox_inches='tight')


#%% Plot the Ratio

cmap        = 'cmo.balance'#cm.cm.roma_r#'cmo.balance'
cints_ratio = np.arange(0,1.05,.05)

for ex in range(nexps):
    fig,ax      = init_globalmap()
    
    plotvar     = lpstd_sst[ex]/rawstd_sst[ex]
    pcm         = ax.contourf(plotvar.lon,plotvar.lat,plotvar,
                          levels=cints_ratio,cmap=cmap,
                          transform=cproj,extend='both')
    
    # plotvar  = lpstd_sst[ex]#(dsclims[ex].mean('month') * dsmask180.mask).TS
    cl          = ax.contour(plotvar.lon,plotvar.lat,plotvar,
                          levels=cints_ratio[::2],colors="k",
                          transform=cproj,linewidths=0.55,alpha=0.5)
    
    
    
    ax.clabel(cl)
    
    # Plot Ice
    ax.contour(plotmask.lon,plotmask.lat,plotmask,colors='cyan',
               transform=cproj,linewidths=0.66,zorder=1,linestyles='solid')
    
    cb       = viz.hcbar(pcm,ax=ax)
    cb.set_label(r"$\sigma(TS)$ Ratio $\frac{Low  \,\, Pass}{Raw}$",fontsize=fsz_axis)
    ax.add_feature(cfeature.LAND,zorder=1,color="k")
    ax.set_title(expnames_long[ex],fontsize=fsz_title)
    
    figname = "%sTS_Stdev_Ratio_Global_%s.png" % (figpath,expnames[ex])
    plt.savefig(figname,dpi=150,bbox_inches='tight')


#%% Try taking the area average over the same SPGeastern region

bbsel        = [-40,-15,52,62]
spge_ssts    = [proc.sel_region_xr(ds,bbsel) for ds in dsmasked]
spge_aavg    = [proc.area_avg_cosweight(ds) for ds in spge_ssts]

spge_std     = [ds.std('time').data.item() for ds in spge_aavg]
spge_aavg_lp = [lpfilter_ts(ts) for ts in spge_aavg]
spge_std_lp  = [np.std(ts) for ts in spge_aavg_lp]

#%% Plot Standard Deviation

#ratios     = lpstds/rawstds

instd       = spge_std
instd_lp    = spge_std_lp

vratio      = np.array(instd_lp) / np.array(instd) * 100

xlabs       = ["%s\n%.2f" % (expnames_short[ii],vratio[ii])+"%" for ii in range(len(vratio))]
fig,ax      = plt.subplots(1,1,constrained_layout=True,figsize=(6,6))
braw        = ax.bar(np.arange(nexps),instd,color='gray')
blp         = ax.bar(np.arange(nexps),instd_lp,color='k')

ax.bar_label(braw,fmt="%.04f",c='gray')
ax.bar_label(blp,fmt="%.04f",c='k')

#ax.bar_label(vratio,fmt="%.2f",c=w,label_type='bottom')

ax.set_xticks(np.arange(nexps),labels=xlabs)
ax.set_ylabel("$\sigma$(SST) [$\degree$C]")

ax.set_ylim([0,1.0])

#figname = "%sTS_AMV_NASST_Index_Stdev_Global.png" % (figpath)
#plt.savefig(figname,dpi=150,bbox_inches='tight')

#%%


    
# figname = "%sTS_Stdev_Mean_Global_%s.png" % (figpath,expnames[ex])
# plt.savefig(figname,dpi=150,bbox_inches='tight')






#
#dsmasked = [ds.TS * dsmask180.mask for ds in dsanoms]

#def compute_amv_quick()




