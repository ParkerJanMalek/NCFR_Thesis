# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 15:12:26 2023

@author: parke
"""
import xarray as xr
import netCDF4 as nc
import numpy as np
import scipy.io as sc
import mat73
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import geopandas
from shapely.geometry import Point


fds = xr.open_dataset("D:/PSU Thesis/data/globalARcatalog_MERRA2_1980-2020_v3.1.nc",chunks='auto')

# defin 
land_lat=np.squeeze(fds.lflat)
land_lon=np.squeeze(fds.lflon)
lathead = np.squeeze(fds.hlat)
lonhead = np.squeeze(fds.hlon)
kid = np.squeeze(fds.kid)
kstatus = np.squeeze(fds.kstatus)

#xmin xmax ymin ymax
latlonbox = [360-124.409591, 360-114.131211, 32.5, 42]
#latlonbox = [180+120, 180+110, 30, 45]
#open AR time series
with open('D:\\PSU Thesis\\data\\AR_California_Landfall.pickle', 'rb') as data:
    date_AR = pickle.load(data)

if(not date_AR):
    # loop through landfall coordinates, pick ones that occur within range of bounding box described above
    # ((lat_no_nan[j] >= latlonbox[2]) & (lat_no_nan[j] <= latlonbox[3])) & 
    kid_no_nan = []
    kstatus_no_nan = []
    lat_no_nan= []
    lon_no_nan = []
    lath_no_nan= []
    lonh_no_nan = []
    
    date_AR = []
    for i in range(0,len(land_lat[:,0])):
        print(i)
        lat = np.array(land_lat[i,:])
        lon = np.array(land_lon[i,:])
        lath = np.array(lathead[i,:])
        lonh = np.array(lonhead[i,:])
        kid_interm = np.array(kid[i,:])
        kstatus_interm = np.array(kstatus[i,:])
        
        
        lat_no_nan = lat[~np.isnan(lat)]
        lon_no_nan = lon[~np.isnan(lon)]
        
        for j in range(0,len(lon_no_nan)):
            if  ((lat_no_nan[j] >= latlonbox[2]) & (lat_no_nan[j] <= latlonbox[3])) & ((lon_no_nan[j] >= latlonbox[0]) & (lon_no_nan[j] <= latlonbox[1])):
                date_AR.append(np.array([lat_no_nan[j],lon_no_nan[j],np.array(fds.kstatus.time)[i].astype('datetime64[h]')]))
    with open("AR_California_Landfall.pickle", "wb") as f:
        pickle.dump(date_AR, f, protocol=pickle.HIGHEST_PROTOCOL)
AR_point = []

for i in range(0,len(date_AR)):
   latlon = [date_AR[i][0],date_AR[i][1]-360]
   AR_point.append(latlon)
AR_pd = pd.DataFrame(AR_point,columns=["latitude","longitude"])

    
    
#plots of landfall points on map


california = geopandas.read_file('D:\\PSU Thesis\\data\\CA_State_TIGER2016.shp')
california =  california.to_crs(4326)
fig,ax = plt.subplots(figsize = (15,15))
california.plot(ax=ax)
geometry = [Point(xy) for xy in zip(AR_pd.longitude,AR_pd.latitude)]
geo_df = geopandas.GeoDataFrame(geometry = geometry)
california.crs = {'init':"epsg:4326"}
geo_df.crs = {'init':"epsg:4326"}
g = geo_df.plot(ax = ax, markersize = 20, color = 'red',marker = '*')
plt.show()


