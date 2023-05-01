# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 18:10:10 2023

Pulls data based from MRMS database based on start and end dates

@author: parke
"""
import datetime as dtetme
import wget
import pygrib
from dateutil import rrule
from datetime import datetime, timedelta
import matplotlib.pyplot as plt 
import gzip
import pygrib
import os
import urllib
import re
import numpy as np
from mpl_toolkits.basemap import Basemap
import nio
import xarray as xr



start_date = dtetme.datetime(2016,1,9,13)
end_date = dtetme.datetime(2016,1,9,13)

bound_ca = True
#bounding box California

if (bound_ca):
    min_lon = 230
    max_lon = 250
    min_lat = 30
    max_lat = 45
else:
   min_lon = 230.005
   max_lon = 299.994998
   min_lat = 20.005
   max_lat = 54.995 
   



file_var = 'GaugeCorr_QPE_01H_00' #GaugeCorr_QPE_01H_00
file_var2 = 'GaugeCorr_QPE_01H' #GaugeCorr_QPE_01H

for dt in rrule.rrule(rrule.HOURLY, dtstart=start_date, until=end_date):
    month = str(dt.month).zfill(2)
    day = str(dt.day).zfill(2)
    hour = str(dt.hour).zfill(2)
    
    filename = file_var+'.00_'+str(dt.year)+str(month)+str(day)+'-'+str(hour)+'0000.grib2.gz'
    url = 'https://mtarchive.geol.iastate.edu/'+str(dt.year)+'/'+str(month)+'/'+str(day)+'/mrms/ncep/'+file_var2+'/'+filename
    if os.path.exists(filename):
        os.remove(filename) 
        
        
    fileupload= wget.download(url)
    
    with gzip.open(fileupload, 'rb') as f:
        uncompressed_data = f.read()
        
        temp_file_path_name = re.sub(".gz","",f.name)
        
    if os.path.exists(temp_file_path_name):
        os.remove(temp_file_path_name)
    
    with open(temp_file_path_name, 'wb') as temp_file:
        temp_file.write(uncompressed_data)
        
    print(temp_file_path_name)
        
        
    grbs = xr.open_dataset(temp_file_path_name,engine='cfgrib')
    # grt=grbs[1]
    # value = grt.values
    lat_min_ind = np.where(np.array(grbs.latitude) ==np.array(grbs.latitude.min()))[0][0]
    lat_max_ind =  np.where(np.array(grbs.latitude) ==np.array(grbs.latitude.max()))[0][0]
    lon_min_ind = np.where(np.array(grbs.longitude) ==np.array(grbs.longitude.min()))[0][0]
    lon_max_ind = np.where(np.array(grbs.longitude) ==np.array(grbs.longitude.max()))[0][0]
    ca_values = grbs.unknown[lat_max_ind:lat_min_ind+1,lon_min_ind:lon_max_ind+1]
    #     #---Convert latitude/longitude 1D arrays to 2D.
    lat = grbs.latitude
    lon = grbs.longitude
    lon2d, lat2d = np.meshgrid(lon,lat)
    
    
    fig = plt.figure(figsize=(8,8))
    ax  = fig.add_axes([0.1,0.1,0.8,0.9])

    # Define and plot the meridians and parallels
    min_lat = min_lat
    max_lat = max_lat
    min_lon = min_lon
    max_lon = max_lon
    

    # Create the basemap object
    bm = Basemap(projection="cyl",
                  llcrnrlat=min_lat-1,
                  urcrnrlat=max_lat+1,
                  llcrnrlon=min_lon-1,
                  urcrnrlon=max_lon+1,
                  resolution='l')

    
    bm.drawcoastlines()
    bm.drawstates()
    bm.drawcountries()
    cf  = bm.contourf(lon2d,lat2d,ca_values)
    cb  = bm.colorbar(cf,"bottom", size="7%", pad="10%")
    # fig, ax = plt.subplots(figsize=(8,8))
    # im = ax.imshow(ca_values, extent=(min_lon,max_lon, min_lat, max_lat))
    # ax.set_title(grt.name)
    # cbar = plt.colorbar(im, ax=ax)
    # cbar.ax.set_ylabel(grt.units)
    # plt.show()
    grbs.close()
    
    
    if os.path.exists(temp_file_path_name):
        os.remove(temp_file_path_name)
        os.remove(filename) 
        