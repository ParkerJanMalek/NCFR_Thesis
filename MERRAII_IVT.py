# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Plotting ACTUAL IVT during extreme precipitation events to compare to SOM patterns

UPDATED 6/12/2023
"""
#%% IMPORT MODULES
import netCDF4 as nc #if this doesn't work, try to reinstall in anaconda prompt using
    #conda install netcdf4
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.cm as cm
#import matplotlib.transforms as mtransforms
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
#from mpl_toolkits.basemap import Basemap #installed using 
    #conda install -c anaconda basemap
#import scipy.io
from datetime import datetime
import os
import paramiko

from matplotlib.dates import DayLocator, DateFormatter

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pyproj import Geod #this is used to convert range/azimuth to lat/lon

import matplotlib.axes as maxes

from dateutil import rrule
import datetime  as dtetme
import wget

def pull_merra(start_date1,end_date1,station_data,ts_selected,station_name,station_lon,station_lat):
    start_date = dtetme.datetime(start_date1.year,start_date1.month,start_date1.day,start_date1.hour)
    end_date = dtetme.datetime(end_date1.year,end_date1.month,end_date1.day,end_date1.hour)
    
    for dt in rrule.rrule(rrule.HOURLY, dtstart=start_date, until=end_date):
        outdir = 'G:/NCFR Thesis/NCFR_Thesis/IVT_'+station_name + '_'+str(start_date.year)+str(start_date.month)+str(start_date.day)+ str(end_date.hour)+'_'+ str(end_date.year)+ str(end_date.month)+ str(end_date.day)+ str(end_date.hour) +'\\'
        # url = 'https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXSLV.5.12.4/2017/02/MERRA2_400.tavg1_2d_int_Nx.'+str(dt.year)+str("{:02d}".format(dt.month))+str("{:02d}".format(dt.day))+'.nc4'
        # fileupload= wget.download(url,out="")
        #"scp malek@circe.rc.pdx.edu:/vol/share/climate_lab2/MERRA2/Daily_and_Subdaily/IVT_hourly/MERRA2.README.pdf D:\PSU Thesis\data\"Research Statement High-Intensity Precipitation.docx
        #"scp malek@circe.rc.pdx.edu:/vol/share/climate_lab2/Parker/Papers/"Research Statement High-Intensity Precipitation.docx" C:\Users\malekP\
        #%% IMPORT EXTREME DAYS DATA
        # change directory and import SOM data from .mat file
        #mat_dir='I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomOutput'
        # define lat, lon region of data for plotting
        latmin, latmax = (15.5,65.5)
        lonmin, lonmax = (-170.25,-105.75)
        
        #%% IMPORT MERRA2 DATA
        # define metvar
        metvars = ['SLP', '300W','Z500Anom','SLPAnom','Z850','850T','850TAnom']
        metvars = ['850TAnom','IVT']#]
        #metvar = '300W'
        for metvar in metvars:
            if metvar == 'IVT':
                filepath = "G:/NCFR Thesis/NCFR_Thesis/data/MERRA2_400.tavg1_2d_int_Nx."+str(dt.year)+str("{:02d}".format(dt.month))+str("{:02d}".format(dt.day))+".SUB.nc"#G:\\NCFR Thesis\\NCFR_Thesis\\MERRA2_400.tavg1_2d_int_Nx."+dt.20170207.SUB.nc"
            elif metvar == '850TAnom':
                filepath = "G:/NCFR Thesis/NCFR_Thesis/data/MERRA2_400.tavg1_2d_slv_Nx."+str(dt.year)+str("{:02d}".format(dt.month))+str("{:02d}".format(dt.day))+".nc4"
            #COLLECT VARIABLE DATA FROM MERRA2 FILE
            merravar = {'Z500':'H','SLP':'SLP','850T':'T','Z850':'H'}
            #open the netcdf file in read mode
            gridfile = nc.Dataset(filepath,mode='r')
            gridlat = gridfile.variables['lat'][:]
            gridlon = gridfile.variables['lon'][:]
            if metvar == 'IVT':
                Uvapor = gridfile.variables['UFLXQV'][:]
                Vvapor = gridfile.variables['VFLXQV'][:]
                merra = np.sqrt(Uvapor**2 + Vvapor**2)
            elif metvar == '850TAnom': #temperature advection
                UT = gridfile.variables['U850'][:]
                VT = gridfile.variables['V850'][:]
                T = gridfile.variables['T850'][:]
                proj=ccrs.LambertConformal(central_longitude=-90)
                lon,lat=np.meshgrid(gridlon,gridlat)
                output=proj.transform_points(ccrs.PlateCarree(),lon,lat)
                x,y=output[:,:,0],output[:,:,1]
                gradx=np.gradient(x,axis=1)
                grady=np.gradient(y,axis=0)
                merra=-(UT*(np.gradient(T,axis=1)/gradx)+VT*(np.gradient(T,axis=0)/grady))*3600
            #     merra = np.sqrt(Uvapor**2 + Vvapor**2)
            gridfile.close()
            
            merra = np.squeeze(merra)
            
           
            #%% REDUCE LAT AND LON TO DESIRED AREA
            
            #REDUCE VARIABLES TO DESIRED AREA
            #reduce lat
            latlims = np.logical_and(gridlat > latmin, gridlat < latmax)
            latind = np.where(latlims)[0]
            gridlatreduced = gridlat[latind]
            #reduce lon
            lonlims = np.logical_and(gridlon > lonmin, gridlon < lonmax)
            lonind = np.where(lonlims)[0]
            gridlonreduced = gridlon[lonind]
            #reduce pressure
            merrareduced = merra[:,latind,:]
            merrareduced = merrareduced[:,:,lonind]
            
            #print(np.amin(merrareduced),np.amax(merrareduced))
            
            
            #%% CREATE ANOMALY MAP 
            #GENERATE CUSTOM COLORMAP
            def center_colormap(lowlim, highlim, center=0):
                dv = max(-lowlim, highlim) * 2
                N = int(256 * dv / (highlim-lowlim))
                bwr = cm.get_cmap('seismic', N)
                newcolors = bwr(np.linspace(0, 1, N))
                beg = int((dv / 2 + lowlim)*N / dv)
                end = N - int((dv / 2 - highlim)*N / dv)
                newmap = ListedColormap(newcolors[beg:end])
                return newmap
            
            #%% DEFINE PLOTTING VARIABLES
            minmax = 0.002 *3600 #s to hour
            lowanom, highanom = (-minmax, minmax)
            newmap = center_colormap(lowanom, highanom, center=0)
            lowlims = {'Z500':2850,'SLP':985,'IVT':0,'300W':0,'850T':252,'Z500Anom':lowanom,'Z850':1187,'SLPAnom':lowanom,'850TAnom':-minmax}
            highlims = {'Z500':5700,'SLP':1022,'IVT':1700,'300W':56,'850T':293,'Z500Anom':highanom,'Z850':1548,'SLPAnom':highanom,'850TAnom':minmax}
            
            contourstart = {'Z500':3000,'SLP':990,'IVT':0,'300W':5,'850T':250,'Z500Anom':-1.75,'Z850':1190,'SLPAnom':-2.25,'850TAnom':-minmax}
            contourint = {'Z500':200,'SLP':4,'IVT':100,'300W':5,'850T':2.5,'Z500Anom':0.25,'Z850':30,'SLPAnom':0.25,'850TAnom':minmax/6}
            
            cbarstart = {'Z500':3000,'SLP':990,'IVT':0,'300W':0,'850T':250,'Z500Anom':-2.0,'Z850':1200,'SLPAnom':-2.4,'850TAnom':-minmax}
            cbarint = {'Z500':500,'SLP':5,'IVT':150,'300W':10,'850T':5,'Z500Anom':0.5,'Z850':50,'SLPAnom':0.4,'850TAnom':minmax/6}
            
            colormap = {'Z500':'jet','SLP':'rainbow','IVT':'gnuplot2_r','300W':'hot_r','850T':'turbo','Z500Anom':newmap,'Z850':'turbo','SLPAnom':newmap,'850TAnom':'coolwarm'}
            cbarlabs = {'Z500':'m','SLP':'hPa','IVT':'kg $\mathregular{m^{-1}}$ $\mathregular{s^{-1}}$','300W':'m/s','850T':'K','Z500Anom':r'$\mathbf{\sigma}$','Z850':'m','SLPAnom':r'$\mathbf{\sigma}$','850TAnom':'Degrees/hr'}
            plottitle = {'Z500':'Z500','SLP':'SLP','IVT':'IVT','300W':'300 hPa Wind','850T':'850 hPa Temperature Advection','Z500Anom':'Z500 Anomaly','Z850':'Z850','SLPAnom':'SLP Anomaly','850TAnom':'850 hPa Temperature Advection'}
            #%% PLOT NODES from MATLAB
            
            #create subplot for mapping multiple timesteps
            #MAP DESIRED VARIABLE
            # define date of plot
        
            # Read in NetCDF4 file (add a directory path if necessary):
        
            #Start Plotting Data
            
            # Plot the data using matplotlib and cartopy
            #n=1 #testing
            #for dt in rrule.rrule(rrule.HOURLY, dtstart=start_date, until=end_date):
            #for n in np.arange(1,24):
                
                
            datetitle =  plottitle[metvar]+" on " + str(dt.year)+"-"+str(dt.month)+"-"+str(dt.day)+ "-" + str(dt.hour) +":00 UTC"
            figtitle = plottitle[metvar]+" on " + str(dt.year)+"-"+str(dt.month)+"-"+str(dt.day)+ "-" + str(dt.hour) 
            # Set the figure size, projection, and extent
            fig = plt.figure(figsize=(20, 20))
            
            
            
            ax1 = fig.add_subplot(211)
            #fig, ax = plt.subplots(figsize=(20, 20))
            # Plot the data!
            date_filter = station_data[ts_selected['start']:ts_selected['end']]
            ax1.bar(date_filter.index,date_filter,width=0.01)
            ax1.axvline(x=dt,linewidth=4, color='r')
            plt.ylabel('Precipiation (mm)', fontsize=25)
            plt.xticks(fontsize=15,rotation=40)
            plt.yticks(fontsize=15)
            ax1.grid()
            date_form = DateFormatter("%m-%d-%Y-%H")
            ax1.xaxis.set_major_formatter(date_form)
            
            
            ax1.set_title(datetitle, fontsize=40)
            
            
            
            
            usemap_proj = ccrs.PlateCarree(central_longitude=180)
            usemap_proj._threshold /= 20.  # to make greatcircle smooth
            
            ax2 = fig.add_subplot(212,projection=usemap_proj)
            ax2.set_global()
            
            border_c = '0.4'
            border_w = 12
            ax2.coastlines(resolution="110m",linewidth=1)
            ax2.gridlines(linestyle='--',color='black')
            ax2.add_feature(cfeature.LAND)
            ax2.add_feature(cfeature.OCEAN,color="white")
            ax2.add_feature(cfeature.COASTLINE)
            ax2.add_feature(cfeature.BORDERS, linestyle=':', zorder=3)
            ax2.add_feature(cfeature.STATES, linestyle=':', zorder=3)
            gl = ax2.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), color='gray', linewidth=0.3)
            gl.top_labels = False
            gl.right_labels = False
            gl.xlabel_style = {'size': 25}
            gl.ylabel_style = {'size': 25}
            # set appropriate extents: (lon_min, lon_max, lat_min, lat_max)
            ax2.set_extent([lonmin, lonmax, latmin, latmax], crs=ccrs.PlateCarree())
            ax2.plot(station_lon,station_lat,marker='o', color='red', markersize=15, transform=ccrs.PlateCarree())
            
            #define area threshold for basemap
            area_thresh = 1E4
            #create equidistant cylindrical projection basemap
            
            # usemap_proj = ccrs.PlateCarree(central_longitude=180)
            # usemap_proj._threshold /= 20.  # to make greatcircle smooth
            # ax = fig.add_subplot(212,projection=usemap_proj)
            # ax.axis('off')
            #ax = plt.axes(projection=usemap_proj)
            # set appropriate extents: (lon_min, lon_max, lat_min, lat_max)
            # Set contour levels, then draw the plot and a colorbar
            arr = merrareduced[dt.hour,:,:]
            contour_c = '0.1'
            contour_w = 0.7
            lon, lat = np.meshgrid(gridlonreduced,gridlatreduced) 
            plt.scatter(-120.9,39.5,color='r',marker='*',linewidths=5)
            colorm = plt.pcolor(lon,lat,arr,shading='auto',cmap=colormap[metvar],vmin=lowlims[metvar],vmax=highlims[metvar])
            mp =  plt.contourf(lon, lat, arr, np.arange(contourstart[metvar],highlims[metvar],contourint[metvar]), transform=ccrs.PlateCarree(),cmap=colormap[metvar])
            mp2 = plt.contour(lon,lat,arr,colors=contour_c,linewidths=contour_w,levels=np.arange(contourstart[metvar],highlims[metvar]+1,contourint[metvar]))
            
            cbar = plt.colorbar(mp,ticks=np.arange(cbarstart[metvar],highlims[metvar],cbarint[metvar]),orientation='vertical',pad=0.08)
            cbar.set_label(cbarlabs[metvar],fontsize=20,labelpad=0.5,fontweight='bold')
            cbar.ax.tick_params(labelsize=20)
            fig.tight_layout()
            plt.show()
            fig.savefig(outdir+figtitle+".png")
            #plt.close()
                
    #contourm = map.contour(xi,yi,arr,colors=contour_c,linewidths=contour_w,levels=np.arange(contourstart['IVT'],highlims['IVT']+1,contourint['IVT']),zorder=2)
    # Save the plot as a PNG image
    
    #fig.savefig('MERRA2_t2m.png', format='png', dpi=360)
     
    
    
    
    
    
    















   
    # for n in np.arange(1,24):
    #     fig = plt.figure()
    #     datetitle =  "IVT on 2017-02-07 - " + str(i+1) +":00 UTC"
    #     # reduce merra to desired day
    #     arr = merrareduced[1,:,:]
    #     #convert lat and lon into a 2D array
    #     lon, lat = np.meshgrid(gridlonreduced,gridlatreduced) 
      
    #     geodetic = ccrs.Geodetic()
    #     plate_carree = ccrs.PlateCarree(central_longitude=180)
         
        
        
    #     map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
    #               urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh)
    #     xi, yi = map(lon,lat)
    #     ax = fig.add_subplot()
    #     ax.set_title(datetitle,pad=4,fontsize=12)
    #     #ax.set_title('{:%d %b}'.format(datetitle),pad=4,fontsize=12)
    #     # sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
    #     # ax.text(0.0, 1.0, extremeevent[i], transform=ax.transAxes + sublabel_loc,
    #     #     fontsize=9, fontweight='bold', verticalalignment='top', 
    #     #     bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
    #     #create colormap of MERRA2 data
    #     colorm = map.pcolor(xi,yi,arr,shading='auto',cmap=colormap['IVT'],vmin=lowlims['IVT'],vmax=highlims['IVT'],zorder=1)
        
    #     #define border color and thickness
    #     border_c = '0.4'
    #     border_w = 0.4
    #     #create map features
    #     map.drawcoastlines(color=border_c, linewidth=border_w)
    #     map.drawstates(color=border_c, linewidth=border_w)
    #     map.drawcountries(color=border_c, linewidth=border_w)
    #     gridlinefont = 8.5
    #     parallels = np.arange(20.,71.,20.)
    #     meridians = np.arange(-160.,-109.,20.)
    #     map.drawparallels(parallels, labels=[1,0,0,0], fontsize=gridlinefont,color=border_c,linewidth=border_w)
    #     map.drawmeridians(meridians, labels=[0,0,0,1], fontsize=gridlinefont,color=border_c,linewidth=border_w)
    #       #define contour color and thickness
    #     contour_c = '0.1'
    #     contour_w = 0.7
    #     #create contour map
    #     contourm = map.contour(xi,yi,arr,colors=contour_c,linewidths=contour_w,levels=np.arange(contourstart['IVT'],highlims['IVT']+1,contourint['IVT']),zorder=2)
    #     plt.clabel(contourm,levels=contourm.levels[::2],fontsize=6,inline_spacing=1,colors='k',zorder=2,manual=False)
            
    #     #add yuba shape
    #     #map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed,linewidth=0.8,color='r')
    #     plt.scatter(-120.9,39.5,color='r',marker='*',linewidths=0.7,zorder=4)
    #     #cbar_ax = fig.add_axes([0.9,0.05,0.025,0.88]) #bottom colorbar
    #     cbar = fig.colorbar(colorm,ticks=np.arange(cbarstart['IVT'],highlims['IVT']+1,cbarint['IVT']),orientation='vertical')
    #     cbar.ax.tick_params(labelsize=8)
    #     cbar.set_label(cbarlabs['IVT'],fontsize=8.5,labelpad=0.5,fontweight='bold')
            
    #     # #CUSTOMIZE SUBPLOT SPACING
    #     # fig.subplots_adjust(left=0.05,right=0.89,bottom=0.021, top=0.955,hspace=0.05, wspace=0.05) #bottom colorbar
    #     # #fig.add_axis([left,bottom, width,height])
    #     # cbar_ax = fig.add_axes([0.9,0.05,0.025,0.88]) #bottom colorbar
    #     # cbar = fig.colorbar(colorm, cax=cbar_ax,ticks=np.arange(cbarstart['IVT'],highlims['IVT']+1,cbarint['IVT']),orientation='vertical')
    #     # cbar.ax.tick_params(labelsize=8)
    #     # cbar.set_label(cbarlabs['IVT'],fontsize=8.5,labelpad=0.5,fontweight='bold')
        
            
    #     #SHOW MAP
    #     fig.savefig("G:\\NCFR Thesis\\NCFR_Thesis\\IVT-2017-02-03 "+str(i+1)+".png",dpi=300)
    #     i = i+1
    #     plt.show()