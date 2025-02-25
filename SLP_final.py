# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 14:13:40 2022

@author: emrus2

This script is used to create a subplot of sea level pressure
for 3 timesteps of each day for three days, colormap and contours

bounds are [990,1040]

UPDATED ON 1/26/2023
"""

#IMPORT MODULES
import netCDF4 as nc #if this doesn't work, try to reinstall in anaconda prompt using
    #conda install netcdf4
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
from mpl_toolkits.basemap import Basemap #installed using 
    #conda install -c anaconda basemap
    
#DEFINE DATES OF ANALYSIS 
#define hour, day, month, year of analysis
year = '2020'
month = '09'
#define integers of date
year_i,month_i = (int(year[2:4]),int(month))


#create subplot for mapping multiple timesteps
fig = plt.figure(figsize=(7,7.5))
fig.suptitle('SLP',fontsize=13,fontweight="bold",y=0.9875)

#define subplot labels
sublabel = [False,'a)','e)','i)','b)','f)','j)','c)','g)','k)','d)','h)','l)']

n = 1
for day_i in np.arange(6,9,1):
    i = n
    #define string version of days
    day = str(day_i)
    if len(day) < 2:
        day = '0' + day

    #define MERRA2 location
    filename = (f'MERRA2.inst3_3d_asm_Np.{year}{month}{day}.SUB.nc')
    base_dir = 'I:\\MERRA2\\Daily_and_Subdaily\\Sea_level_pressure_3hourly'
    filepath = os.path.join(base_dir,filename)
    
    
    #COLLECT VARIABLE DATA FROM MERRA2 FILE
    #open the netcdf file in read mode
    gridfile = nc.Dataset(filepath,mode='r')
    gridlat = gridfile.variables['lat'][:]
    gridlon = gridfile.variables['lon'][:]
    time = gridfile.variables['time'][:]
    pressure = gridfile.variables['SLP'][:] #time x height x lat x lon
    gridfile.close()
    
    #CONVERT SLP INTO MB
    pressuremb = pressure/100
    
    
    #REDUCE VARIABLES TO DESIRED AREA
    #convert height to a 3D array
    pressuresm = np.squeeze(pressuremb)
    #print(np.shape(heightsm)) #output shows 8x361x576, time,lat,lon
    #define lat lon restrictions
    latmin = 29.5
    latmax = 60.5
    lonmin = -142.5
    lonmax = -101.5
    #reduce lat
    latlims = np.logical_and(gridlat > latmin, gridlat < latmax)
    latind = np.where(latlims)[0]
    gridlatreduced = gridlat[latind]
    #reduce lon
    lonlims = np.logical_and(gridlon > lonmin, gridlon < lonmax)
    lonind = np.where(lonlims)[0]
    gridlonreduced = gridlon[lonind]
    #reduce pressure
    pressurereduced = pressuresm[:,latind,:]
    pressurereduced = pressurereduced[:,:,lonind]
    
    
#%%
    #LOOP THROUGH EACH TIMESTEP
    for timestep in np.arange(0,8,2): #plotting 3,12,21Z timesteps
        zulu = str(3*timestep)
        if len(zulu) < 2:
            zulu = '0' + zulu
    
        #REDUCE VARIABLES TO DESIRED TIMESTEP
        timenew = time[timestep]
        pressurenew = pressurereduced[timestep,:,:]
    
        #MAP DESIRED VARIABLE
        #convert lat and lon into a 2D array
        lon, lat = np.meshgrid(gridlonreduced,gridlatreduced) 
        #define area threshold for basemap
        area_thresh = 1E4
        #create equidistant cylindrical projection basemap
        map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
                  urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh)
        xi, yi = map(lon,lat)
        print(i)
        ax = fig.add_subplot(4,3,i)
        #add zulu labels on top row
        if day_i == 6:
            ax.set_ylabel(f'{zulu} UTC',fontsize=8.5,fontweight="bold",labelpad=0.3)
        #add sublabels to single plots
        sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
        ax.text(0.0, 1.0, sublabel[i], transform=ax.transAxes + sublabel_loc,
            fontsize=8, verticalalignment='top', 
            bbox=dict(facecolor='1', edgecolor='none', pad=1.5))
        #add date labels on left column
        if timestep == 0:
            ax.set_title(f'{day} SEP',fontsize=10,fontweight="bold",pad=2)
        #define z-axis bounds for colormap
        lowlim = 990
        highlim = 1040
        #create colormap of MERRA2 data
        colorm = map.pcolor(xi,yi,pressurenew,shading='auto',cmap='rainbow',vmin=lowlim,vmax=highlim)
        #define border color and thickness
        border_c = '0.4'
        border_w = 0.4
        #create map features
        map.drawcoastlines(color=border_c, linewidth=border_w)
        map.drawstates(color=border_c, linewidth=border_w)
        map.drawcountries(color=border_c, linewidth=border_w)
        gridlinefont = 7
        parallels = np.arange(35.,60.,10.)
        meridians = np.arange(-134.,-109.,12.)
        if i == 3 or i == 6 or i == 9:
            map.drawparallels(parallels, labels=[0,1,0,0], fontsize=gridlinefont,color=border_c,linewidth=border_w)
            map.drawmeridians(meridians,color=border_c,linewidth=border_w)
        elif i == 10 or i == 11:
            map.drawparallels(parallels, color=border_c,linewidth=border_w)
            map.drawmeridians(meridians, labels=[0,0,0,1], fontsize=gridlinefont,color=border_c,linewidth=border_w)
        elif i == 12:
            map.drawparallels(parallels, labels=[0,1,0,0], fontsize=gridlinefont,color=border_c,linewidth=border_w)
            map.drawmeridians(meridians, labels=[0,0,0,1], fontsize=gridlinefont,color=border_c,linewidth=border_w)
        else:
            map.drawparallels(parallels, color=border_c,linewidth=border_w)
            map.drawmeridians(meridians,color=border_c,linewidth=border_w)
        #define contour color and thickness
        contour_c = '0.1'
        contour_w = 0.5
        #create contour map
        contourm = map.contour(xi,yi,pressurenew,levels=np.arange(lowlim,highlim+1,4),colors=contour_c,linewidths=contour_w)
        
        i += 3
    n += 1

#CUSTOMIZE SUBPLOT SPACING
fig.subplots_adjust(left=0.04,right=0.95,bottom=0.083, top=0.94,hspace=0.05, wspace=0.05) #bottom colorbar
#fig.add_axis([left,bottom, width,height])
cbar_ax = fig.add_axes([0.15,0.04,0.7,0.0216]) #bottom colorbar
cbar = fig.colorbar(colorm, cax=cbar_ax,ticks=np.arange(lowlim,highlim+1,10),orientation='horizontal')
cbar.ax.tick_params(labelsize=7)
cbar.set_label('hPa',fontsize=7.5,labelpad=-0.9)

#SHOW MAP
save_dir='I:\\Emma\\LaborDayWildfires\\Figures\\1FinalFigures'
os.chdir(save_dir)
#plt.savefig('F4_SLP.png',dpi=300)
plt.show()