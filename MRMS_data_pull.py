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

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import Normalize
import shutil

def map_event(start_date,end_date,station_lon,station_lat,station_name):

    
    #bounding box California
    min_lon = 360+station_lon - 2
    max_lon = 360+station_lon + 2
    min_lat = station_lat - 2
    max_lat = station_lat + 2
    
    
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
       
    
    
    
    file_var = 'RadarOnly_QPE_01H_00' #PrecipRate_00, GaugeCorr_QPE_01H_00,RadarOnly_QPE_01H_00
    file_var2 = 'RadarOnly_QPE_01H' #GaugeCorr_QPE_01H,PrecipRate,RadarOnly_QPE_01H
    
    outdir = 'D:\\PSU Thesis\\data\\'+station_name+'_'+file_var2 + '_'+str(start_date.year)+str(start_date.month)+str(start_date.day)+ str(end_date.hour)+'_'+ str(end_date.year)+ str(end_date.month)+ str(end_date.day)+ str(end_date.hour) +'\\'
    
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
        cmap = mcolors.ListedColormap(cmap_data, 'precipitation')
        norm = mcolors.BoundaryNorm(clevs, cmap.N)
        
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
        
        
        
        fig = plt.figure(figsize=(12,7))
    
        # this declares a recentered projection for Pacific areas
        usemap_proj = ccrs.PlateCarree(central_longitude=180)
        usemap_proj._threshold /= 20.  # to make greatcircle smooth
    
        ax = plt.axes(projection=usemap_proj)
        # set appropriate extents: (lon_min, lon_max, lat_min, lat_max)
        ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())
    
        geodetic = ccrs.Geodetic()
        plate_carree = ccrs.PlateCarree(central_longitude=180)
    
    
    
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.OCEAN,color="white")
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.BORDERS, linestyle=':', zorder=2)
        ax.add_feature(cfeature.STATES, linestyle=':', zorder=2)
        # plot grid lines
        ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), color='gray', linewidth=0.3)
        
        
        cm = ax.contourf(lon2d, lat2d, ca_values,clevs,cmap=cmap_smooth,
                         transform=ccrs.PlateCarree(), zorder=1)
        ax.plot(station_lon,station_lat,marker='o', color='red', markersize=10, transform=ccrs.PlateCarree())
        
        ax.tick_params(axis='both', which='both', labelsize=8, direction='out')
        norm1 = mcolors.Normalize(vmin=0, vmax=1)
        # colorbar and labels
        cb = plt.colorbar(cm,orientation="horizontal",cmap=cmap_smooth)
        ax.set_title('Gauge-Corrected QPE');
        
        # Add a label to the color bar
        cb.set_label('mm')
        
        #cb  = ax.colorbar(cf,"bottom", size="7%", pad="10%",fig=fig,ax=ax)
        
        plt.show()
        
     
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
        


        