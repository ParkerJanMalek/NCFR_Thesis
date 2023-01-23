# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 15:12:26 2023

@author: parke
"""
import xarray as xr
import netCDF4 as nc
import numpy as np
import scipy.io
import mat73

fds = xr.open_dataset("D:/PSU Thesis/data/globalARcatalog_MERRA2_1980-2020_v3.1.nc",chunks='auto')



land_lat=np.squeeze(fds.lflat)
land_lon=np.squeeze(fds.lflon)
lathead = np.squeeze(fds.hlat)
lonhead = np.squeeze(fds.hlon)
kid = np.squeeze(fds.kid)
kstatus = np.squeeze(fds.kstatus)

#xmin xmax ymin ymax
latlonbox = [360-124.409591, 360-114.131211, 32.5, 42]
#latlonbox = [180+120, 180+110, 30, 45]

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
            date_AR.append(np.array([lat_no_nan[j],lon_no_nan[j],np.array(fds.kstatus.time)[i].astype('M8[D]')]))
               
        
        


