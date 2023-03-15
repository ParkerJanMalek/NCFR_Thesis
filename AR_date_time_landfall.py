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
import os
from shapely.geometry import Point


fds = xr.open_dataset("D:/PSU Thesis/data/globalARcatalog_MERRA2_1980-2020_v3.1.nc",chunks='auto')

# defin 
#tlat=np.squeeze(fds.tlat)
#tlon=np.squeeze(fds.tlon)
#lat=np.squeeze(fds.lat)
#lon=np.squeeze(fds.lon)
month = np.squeeze(fds.month)
day = np.squeeze(fds.day)
year = np.squeeze(fds.year)
hour = np.squeeze(fds.hour)
#kid = np.squeeze(fds.kid)
#kstatus = np.squeeze(fds.kstatus)

land_loc=np.squeeze(fds.lfloc)
lathead = np.squeeze(fds.lflat)
lonhead = np.squeeze(fds.lflon)
#shape = np.squeeze(fds.shape)
kid = np.squeeze(fds.kid)
kstatus = np.squeeze(fds.kstatus)
kidmap = np.squeeze(fds.kidmap)
#xmin xmax ymin ymax
#latlonbox = [360-124.409591, 360-114.131211, 32.5, 42]
latlonbox = [230, 250, 32.5, 41]
#latlonbox = [100,260,0,60]

#latlon_where = np.where(((lathead >= latlonbox[2]) & (lathead <= latlonbox[3])) & ((lonhead >= latlonbox[0]) & (lonhead <= latlonbox[1])) )
#latlon_where = np.where(((land_loc >= latlonbox[2]) & (land_loc <= latlonbox[3])) & ((land_loc >= latlonbox[0]) & (land_loc <= latlonbox[1])) )
for i in range(0,len(kid[:,0])):
    t = np.where(kid[i,:] == '255856')
    if(t):
        print(t)
# lat_CA = []
# lon_CA = []
# time = []
# for i in range(0,len(land_loc[:,0,0])):
#     ind = np.where(~np.isnan(land_loc[i,:,:]))
#     for j in range(0,len(ind[0])):
#         if((ind[0][j] >= latlonbox[2]) & (ind[0][j] >= latlonbox[3]):
#             lat_CA.append(ind[0][j])
#         if((ind[0][j] >= latlonbox[0]) & (ind[0][j] >= latlonbox[1]):
#             lon_CA.append(ind[0][j])
#     latlon_CA.append(np.array([lathead[ind],lonhead[ind]]))
#     time.append(np.array([year[ind],month[ind],day[ind],hour[ind]]))
    
# latlon_where = np.where(((lathead >= latlonbox[2]) & (lathead <= latlonbox[3])) & ((lonhead >= latlonbox[0]) & (lonhead <= latlonbox[1])) )


latlon_CA = []
time = []
for i in range(0,len(latlon_where[0])):
    ind = latlon_where[0][i],latlon_where[1][i]
    latlon_CA.append(np.array([lathead[ind],lonhead[ind]]))
    time.append(np.array([year[ind],month[ind],day[ind],hour[ind]]))
                  
    
# landf_loc = land_loc[:,np.where((land_loc.lat >= latlonbox[2]) & (land_loc.lat <= latlonbox[3]))[0],np.where((land_loc.lon >= latlonbox[0]) & (land_loc.lon <= latlonbox[1]))[0]]
# saved_AR_path = 'D:\\PSU Thesis\\data\\AR_California_Landfall.pickle'
# #latlonbox = [180+120, 180+110, 30, 45]
# #open AR time series
# with open(saved_AR_path, 'rb') as data:
#     date_AR = pickle.load(data)

# if(not os.path.exists(saved_AR_path)):
#     # loop through landfall coordinates, pick ones that occur within range of bounding box described above
#     # ((lat_no_nan[j] >= latlonbox[2]) & (lat_no_nan[j] <= latlonbox[3])) & 
#     kid_no_nan = []
#     kstatus_no_nan = []
#     lat_no_nan= []
#     lon_no_nan = []
#     lath_no_nan= []
#     lonh_no_nan = []
    
#     date_AR = []
#     for i in range(0,len(land_lat[:,0])):
#         print(i)
#         lat = np.array(land_lat[i,:])
#         lon = np.array(land_lon[i,:])
#         lath = np.array(lathead[i,:])
#         lonh = np.array(lonhead[i,:])
#         kid_interm = np.array(kid[i,:])
#         kstatus_interm = np.array(kstatus[i,:])
        
        
#         lat_no_nan = lat[~np.isnan(lat)]
#         lon_no_nan = lon[~np.isnan(lon)]
        
#         for j in range(0,len(lon_no_nan)):
#             if  ((lat_no_nan[j] >= latlonbox[2]) & (lat_no_nan[j] <= latlonbox[3])) & ((lon_no_nan[j] >= latlonbox[0]) & (lon_no_nan[j] <= latlonbox[1])):
#                 date_AR.append(np.array([lat_no_nan[j],lon_no_nan[j],np.array(fds.kstatus.time)[i].astype('datetime64[h]')]))
#     with open("AR_California_Landfall.pickle", "wb") as f:
#         pickle.dump(date_AR, f, protocol=pickle.HIGHEST_PROTOCOL)
        
AR_point = []
date_AR = latlon_CA
for i in range(0,len(date_AR)):
    latlon = [date_AR[i][0],date_AR[i][1]]
    AR_point.append(latlon)
AR_pd = pd.DataFrame(AR_point,columns=["latitude","longitude"])

    
    
#plots of landfall points on map
# world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))


# california = geopandas.read_file('D:\\PSU Thesis\\data\\CA_State_TIGER2016.shp')
# california =  california.to_crs(4326)
# fig,ax = plt.subplots(figsize = (30,30))
# california.plot(ax=ax)
# geometry = [Point(xy) for xy in zip(AR_pd.longitude,AR_pd.latitude)]
# geo_df = geopandas.GeoDataFrame(geometry = geometry)
# california.crs = {'init':"epsg:4326"}
# geo_df.crs = {'init':"epsg:4326"}
# g = geo_df.plot(ax = ax, markersize = 200, color = 'red',marker = '*',label="AR Landfall Locations")
# ax.legend(loc="upper left",fontsize=40)
# ax.set_title("Californai AR Landfall Locations")
# plt.show()
# fig.savefig('D:\\PSU Thesis\\data\\California_AR_Landfall.png')


