# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 18:10:10 2023

Pulls data based from MRMS database based on start and end dates

@author: parke
"""
import datetime as dtetme
import wget
from dateutil import rrule
from datetime import datetime, timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt 
import gzip
import os
import urllib
import re
import numpy as np
import cartopy.crs as ccrs
import xarray as xr
import pandas as pd
import cartopy
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.axes as maxes
import shutil
from metpy.io import Level2File
from metpy.plots import add_timestamp, ctables
from matplotlib.dates import DayLocator, DateFormatter

def map_event(start_date,end_date,station_lon,station_lat,station_name,station_data,ts_selected):

    
    #bounding box California
    min_lon = 360+station_lon - 4
    max_lon = 360+station_lon + 4
    min_lat = station_lat - 4
    max_lat = station_lat + 4
    
    
    #bounding box Russia River
    #min_lon = 238
    #max_lon = 243
    #min_lat = 35
    #max_lat = 40
    
    #bounding box Yuba Feather
    #min_lon = 238
    #max_lon = 243
    #min_lat = 35
    #max_lat = 40
    
    #bounding box Santa Ana
    #min_lon = 238
    #max_lon = 243
    #min_lat = 30
    #max_lat = 45

    
    
    # if bound=="ca":
    #     min_lon = 245
    #     max_lon = 250
    #     min_lat = 33
    #     max_lat = 35
    # elif bound=="yfrr":
    #     min_lon = 235
    #     max_lon = 243
    #     min_lat = 37
    #     max_lat = 43
    # elif bound=="sa":
    #     min_lon = 238
    #     max_lon = 243
    #     min_lat = 33
    #     max_lat = 35
    # else:
    #    min_lon = 230.005
    #    max_lon = 299.994998
    #    min_lat = 20.005
    #    max_lat = 54.995 
       
    
    
    
    file_var = 'GaugeCorr_QPE_01H_00' #PrecipRate_00, GaugeCorr_QPE_01H_00,RadarOnly_QPE_01H_00
    file_var2 = 'GaugeCorr_QPE_01H' #GaugeCorr_QPE_01H,PrecipRate,RadarOnly_QPE_01H
    
    outdir = 'G:\\NCFR Thesis\\NCFR_Thesis\\'+station_name+'_'+file_var2 + '_'+str(start_date.year)+str(start_date.month)+str(start_date.day)+ str(end_date.hour)+'_'+ str(end_date.year)+ str(end_date.month)+ str(end_date.day)+ str(end_date.hour) +'\\'
    
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.mkdir(outdir)
    
    for dt in rrule.rrule(rrule.HOURLY, dtstart=start_date, until=end_date):
        month = str(dt.month).zfill(2)
        day = str(dt.day).zfill(2)
        hour = str(dt.hour).zfill(2)
        
        filename = file_var+'.00_'+str(dt.year)+str(month)+str(day)+'-'+str(hour)+'0000.grib2.gz'
        url = 'https://mtarchive.geol.iastate.edu/'+str(dt.year)+'/'+str(month)+'/'+str(day)+'/mrms/ncep/'+file_var2+'/'+filename
        print(url)
        if os.path.exists(filename):
            os.remove(filename) 
            
            
        fileupload= wget.download(url, out="")
        
        with gzip.open(fileupload, 'rb') as f:
            uncompressed_data = f.read()
            
            temp_file_path_name = re.sub(".gz","",f.name)
            
        if os.path.exists(temp_file_path_name):
            os.remove(temp_file_path_name)
        
        with open(temp_file_path_name, 'wb') as temp_file:
            temp_file.write(uncompressed_data)
            
        print(temp_file_path_name)
            
         
            
            # draw filled contours.
        clevs = np.linspace(0, 45, 21)
       # draw filled contours.
        #clevs = [0, 1, 2.5, 5, 7.5, 10, 15, 20, 30, 40,
         #        50, 70, 100, 150, 200, 250, 300, 400, 500, 600, 750]
        # In future MetPy
        # norm, cmap = ctables.registry.get_with_boundaries('precipitation', clevs)
        cmap_data = [(0.0, 0.0, 0.0,0.0),
                     (0.3137255012989044, 0.8156862854957581, 0.8156862854957581),
                     (0.0, 1.0, 1.0),
                     (0.0, 0.8784313797950745, 0.501960813999176),
                     (0.0, 0.7529411911964417, 0.0),
                     (0.501960813999176, 0.8784313797950745, 0.0),
                     (1.0, 1.0, 0.0),
                     (1.0, 0.6274510025978088, 0.0),
                     (1.0, 0.0, 0.0),
                     (1.0, 0.125490203499794, 0.501960813999176),
                     (0.9411764740943909, 0.250980406999588, 1.0),
                     (0.501960813999176, 0.125490203499794, 1.0),
                     (0.250980406999588, 0.250980406999588, 1.0),
                     (0.125490203499794, 0.125490203499794, 0.501960813999176),
                     (0.125490203499794, 0.125490203499794, 0.125490203499794),
                     (0.501960813999176, 0.501960813999176, 0.501960813999176),
                     (0.8784313797950745, 0.8784313797950745, 0.8784313797950745),
                     (0.9333333373069763, 0.8313725590705872, 0.7372549176216125),
                     (0.8549019694328308, 0.6509804129600525, 0.47058823704719543),
                     (0.6274510025978088, 0.42352941632270813, 0.23529411852359772),
                     (0.4000000059604645, 0.20000000298023224, 0.0)]
        # cmap = mcolors.ListedColormap(cmap_data, 'precipitation')
        # norm = mcolors.BoundaryNorm(clevs, cmap.N)
        norm, cmap = ctables.registry.get_with_steps('NWSReflectivity', 5, 5)
        
        # Number of colors for the smoothed color map
        num_colors = len(cmap_data)
    
        # Create a new color map with interpolated colors
        cmap_smooth = mcolors.LinearSegmentedColormap.from_list('smooth_cmap', cmap_data, num_colors)
    
            
        grbs = xr.open_dataset(temp_file_path_name,engine='cfgrib')
        # grt=grbs[1]
        # value = grt.values
        
        lat_ind = np.where((np.array(grbs.latitude) >=np.array(min_lat)) & (np.array(grbs.latitude) <=np.array(max_lat)))[0]
        #lon_min_ind = np.where(np.array(grbs.longitude) >=np.array(min_lon))[0]
        lon_ind = np.where((np.array(grbs.longitude) >=np.array(min_lon)) & (np.array(grbs.longitude) <=np.array(max_lon)))[0]
        ca_values = grbs.unknown[lat_ind,lon_ind]
        
        
        norm = Normalize(vmin=ca_values.min(),vmax=ca_values.max())
        
        
        #     #---Convert latitude/longitude 1D arrays to 2D.
        lat = grbs.latitude[lat_ind]
        lon = grbs.longitude[lon_ind]
        lon2d, lat2d = np.meshgrid(lon,lat)
        
        fig = plt.figure(figsize=(60, 60))
        ax1 = fig.add_subplot(211)
        date_filter = station_data[ts_selected['start']:ts_selected['end']]
        ax1.bar(date_filter.index,date_filter,width=0.01)
        ax1.axvline(x=dt,linewidth=4, color='r')
        fig.suptitle('test' , fontsize=60)
        plt.ylabel('Precipiation (mm)', fontsize=50)
        plt.xlabel('Date', fontsize=50)
        plt.xticks(fontsize=30,rotation=40)
        plt.yticks(fontsize=30)
        ax1.grid()
        date_form = DateFormatter("%m-%d-%Y-%H")
        ax1.xaxis.set_major_formatter(date_form)
        
        ax1.axis('off')
        
        usemap_proj = ccrs.PlateCarree(central_longitude=180)
        usemap_proj._threshold /= 20.  # to make greatcircle smooth
         
        ax = fig.add_subplot(212,projection=usemap_proj)
        # set appropriate extents: (lon_min, lon_max, lat_min, lat_max)
        ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())
         
        geodetic = ccrs.Geodetic()
        plate_carree = ccrs.PlateCarree(central_longitude=180)
        ax.spines['right'].set_visible(False)
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.OCEAN,color="white")
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.BORDERS, linestyle=':', zorder=2)
        ax.add_feature(cfeature.STATES, linestyle=':', zorder=2)
        # plot grid lines
        gl=ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), color='gray', linewidth=0.3)
        gl.xlabel_style = {'size': 25}
        gl.ylabel_style = {'size': 25}
        
        cm = ax.contourf(lon2d, lat2d, ca_values,clevs,cmap=cmap_smooth,
                         transform=ccrs.PlateCarree(), zorder=1)
        ax.plot(station_lon,station_lat,marker='o', color='red', markersize=20, transform=ccrs.PlateCarree())
        
        ax.tick_params(axis='x', labelsize=35)
        ax.tick_params(axis='y', labelsize=35)
        norm1 = mcolors.Normalize(vmin=0, vmax=1)
        # colorbar and labels
        # divider = make_axes_locatable(ax)
         
        # cax = divider.append_axes("right", size="5%", axes_class=maxes.Axes, pad=0.05)
        # cbar = plt.colorbar(cm, cax=cax,cmap=cmap_smooth, orientation='vertical')
        # cbar.set_label(label='mm',size=35)
        # cbar.ax.tick_params(labelsize=35)
        divider = make_axes_locatable(ax)
            
        cax = divider.append_axes("right", size="5%", axes_class=maxes.Axes, pad=0.05)
        cbar = plt.colorbar(cm, cax=cax, orientation='vertical')
        cbar.set_label(label='mm',size=35)
        cbar.ax.tick_params(labelsize=35)
        plt.setp(ax.get_xticklabels(), fontsize=35)
        plt.setp(ax.get_yticklabels(), fontsize=35)
        ax.tick_params(axis='x', labelsize=35)
        ax.tick_params(axis='y', labelsize=35)
        
        plt.suptitle('Gauge-Corrected QPE '+str(dt.year)+str(month)+str(day)+" - "+str(hour) + ":00 UTC", fontsize=50)
        plt.tight_layout()
        plt.show()
        
        #cb = plt.colorbar(cm,orientation="horizontal",cmap=cmap_smooth)
        #ax.set_title('Gauge-Corrected QPE '+str(dt.year)+str(month)+str(day)+'-'+str(hour),fontsize=35);
        #cb.ax.tick_params(labelsize=25)
        # Add a label to the color bar
        #cb.set_label('mm',fontsize=25)
        
        #cb  = ax.colorbar(cf,"bottom", size="7%", pad="10%",fig=fig,ax=ax)
        
        #plt.show()
        
     
        # # Create the basemap object
        # bm = Basemap(projection="cyl",
        #               llcrnrlat=min_lat-1,
        #               urcrnrlat=max_lat+1,
        #               llcrnrlon=min_lon-1,
        #               urcrnrlon=max_lon+1,
        #               resolution='l')
    
        # bm.shadedrelief() 
        # bm.drawcoastlines()
        # bm.drawstates()
        # bm.drawcountries()
        # cbartiks = np.arange(-3,25,3)
        # cf  = bm.contourf(lon2d,lat2d,ca_values)
        # cb  = bm.colorbar(cf,"bottom", size="7%", pad="10%",fig=fig,ax=ax)
        fig.savefig(outdir+file_var2+'_'+str(dt.year)+str(dt.month)+str(dt.day)+str(hour)+'.png')
        # fig, ax = plt.subplots(figsize=(8,8))
        # im = ax.imshow(ca_values, extent=(min_lon,max_lon, min_lat, max_lat))
        # ax.set_title(grt.name)
        # cbar = plt.colorbar(im, ax=ax)
        # cbar.ax.set_ylabel(grt.units)
        # plt.show()
       # grbs.close()
        
        
        if os.path.exists(temp_file_path_name):
            os.remove(temp_file_path_name)
            os.remove(filename) 
            os.remove(temp_file_path_name+".923a8.idx")
        


        