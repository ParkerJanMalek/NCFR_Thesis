"""
Plotting AWS-hosted NEXRAD Level 2 Data
=======================================================

Access NEXRAD radar data via Amazon Web Services and plot with MetPy

Accessing data remotely is a powerful tool for big data, such as NEXRAD radar data.
By accessing it in the cloud, you can save time and space from downloading the data locally.

"""

######################################################################
# Access the data in the AWS cloud. In this example, we're plotting data
# from the Evansville, IN radar, which had convection within its
# domain on 06/26/2019.
# def pull_radar(start_date1, end_date1):
import netCDF4 as nc
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
import boto3
import botocore
from botocore.client import Config
import matplotlib.pyplot as plt
import matplotlib
import metpy.calc as mpcalc
from metpy.units import units
from metpy.io import Level2File
from metpy.plots import add_timestamp, ctables
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from dateutil import rrule
from datetime import datetime, timedelta
import datetime  as dtetme
from matplotlib.dates import DayLocator, DateFormatter
import pickle
import xarray
import os
from matplotlib.animation import FuncAnimation
from IPython import display
import imageio

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pyproj import Geod #this is used to convert range/azimuth to lat/lon

import matplotlib.axes as maxes

def pull_radar(start_date1,end_date1,station_data,station_data_temp,station_data_wind,station_data_direction,station_data_dew,ts_selected,station_name,radar_name,slat,slon):
    start_date = dtetme.datetime(start_date1.year,start_date1.month,start_date1.day,start_date1.hour)
    end_date = dtetme.datetime(end_date1.year,end_date1.month,end_date1.day,end_date1.hour)
    images = []
    s3 = boto3.resource('s3', config=Config(signature_version=botocore.UNSIGNED,
                                            user_agent_extra='Resource'))
    bucket = s3.Bucket('noaa-nexrad-level2')
    outdir = 'G:\\NCFR Thesis\\NCFR_Thesis\\combined_'+station_name + '_'+str(start_date.year)+str(start_date.month)+str(start_date.day)+ str(end_date.hour)+'_'+ str(end_date.year)+ str(end_date.month)+ str(end_date.day)+ str(end_date.hour) +'\\'
    for dt in rrule.rrule(rrule.HOURLY, dtstart=start_date, until=end_date):
        fig_creation = ['presentation']
        paper_dates = ts_selected['event_rainfall'].index[np.argsort(ts_selected['event_rainfall'])][::-1][0:4].sort_values()
        
        
        for ii in fig_creation:
            if (ii == "paper" and ((dt.year in paper_dates.year) and (dt.month in paper_dates.month) and (dt.day in paper_dates.day) and (dt.hour in paper_dates.hour))) or (ii == "presentation"):
                
                #if ii == "paper" and ((dt.year in paper_dates.year) and (dt.month in paper_dates.month) and (dt.day in paper_dates.day) and (dt.hour in paper_dates.hour)):
                month = str(dt.month).zfill(2)
                day = str(dt.day).zfill(2)
                hour = str(dt.hour).zfill(2)
                #print(str(dt.year) + '/' + month + '/' + day +'/'+ str(dt.year) + month + day + '_'+hour)
                # construct a list of days/hours to loop through.
        
                radar_object1 = bucket.objects.filter(Prefix=str(dt.year) + '/' + month + '/' + day + '/KBBX/KBBX'+ str(dt.year) + month + day + '_'+hour)
                #print(len(list(radar_object1)))
                if len(list(radar_object1)) == 0:
                    radar_object1 = bucket.objects.filter(Prefix=str(dt.year) + '/' + month + '/' + day + '/KDAX/KDAX'+ str(dt.year) + month + day + '_'+hour)
                #fig = plt.figure(figsize=(40, 20))
                p_flag=0
                for obj in radar_object1:
                   #if p_flag == 0:
                    #%% IMPORT MERRA2 DATA
                    # define metvar
                    metvars = ['SLP', '300W','Z500Anom','SLPAnom','Z850','850T','850TAdv']
                    metvars = ['IVT','850T','SLP','850TAdv']#]
                    #metvar = '300W'
                    if ii == "paper":
                        fig = plt.figure(figsize=(10, 30))
                        fig.tight_layout(rect=[0.05, 0.05, 1, 0.93])   
                        fig.set_size_inches([60,60])
                    else:
                        fig = plt.figure(figsize=(60, 20))
                     
                    var=0
                    plot_text = ['(b)','(c)','(d)','(e)']
                    for metvar in metvars:
                        pt = plot_text[var]
                        savestr = obj.key.split("/")[-1]
                        
                        
                        
                        if ii == "presentation":
                            ax1 = fig.add_subplot(2,1,2)
                            # Filter data
                            date_filter = station_data[ts_selected['start']:ts_selected['end']]
                            date_filter_temp = station_data_temp[ts_selected['start']:ts_selected['end']]
                            date_filter_wind = station_data_wind[ts_selected['start']:ts_selected['end']]
                            date_filter_direction = station_data_direction[ts_selected['start']:ts_selected['end']]
                            date_filter_dew = station_data_dew[ts_selected['start']:ts_selected['end']]
                            
                            # Plot precipitation as bars
                            bar1 = ax1.bar(date_filter.index, date_filter, width=0.03, label="Precipitation")
                            ax1.axvline(x=dt, linewidth=4, color='m', label='Threshold')
                            ax1.set_ylabel('Precipitation (mm)', fontsize=20)
                            ax1.set_xlabel('Hour', fontsize=20)
                            ax1.tick_params(axis='both', which='major', labelsize=17)
                            ax1.set_ylim([0, max(date_filter.max() * 1.1, 4.5)])
                            ax1.grid()
                            
                            # Add temperature on secondary y-axis
                            ax_temp = ax1.twinx()
                            ax_temp.plot(date_filter_temp.index, date_filter_temp, linestyle="-", color='red', linewidth=6, label="Surface Temperature")
                            ax_temp.set_ylabel('Temperature ($^\circ$C)', fontsize=20)
                            ax_temp.tick_params(axis='y', labelsize=15)
                            dew_mask = np.isfinite(date_filter_dew)
                            # Add dew point on the same secondary y-axis
                            ax_temp.plot(date_filter_dew.index[dew_mask], date_filter_dew[dew_mask], linestyle="-", color='blue', linewidth=6, label="Dew Point Temperature")
                            ax_temp.legend(loc="upper right")
                            #plot precip
                            # bar1 = ax1.bar(date_filter.index,date_filter,width=0.01)
                            # ax1.axvline(x=dt,linewidth=4, color='m')
                            # ax1.set_ylabel('Precipiation (mm)', fontsize=55)
                            # ax1.set_xlabel('Hour', fontsize=50)
                            # ax1.tick_params(axis='both', which='major', labelsize=55)
                            # ax1.set_ylim([0,4.5])
                            # ax1.grid()
                            
                            
                            date_filter_direction_text = date_filter_direction.copy()
                            
                            
                           
                            date_filter_direction_text[(date_filter_direction>=157.5) & (date_filter_direction <202.5)] = 'S'
                            date_filter_direction_text[(date_filter_direction>=202.5) & (date_filter_direction <247.5)] = 'SW'
                            date_filter_direction_text[(date_filter_direction>=247.5) & (date_filter_direction <292.5)] = 'W'
                            date_filter_direction_text[(date_filter_direction>=292.5) & (date_filter_direction <337.5)] = 'NW'
                            date_filter_direction_text[(date_filter_direction>=337.5) | (date_filter_direction <22.5)] = 'N'
                            date_filter_direction_text[(date_filter_direction>=22.5) & (date_filter_direction <67.5)] = 'NE'
                            date_filter_direction_text[(date_filter_direction>=67.5) & (date_filter_direction <112.5)] = 'E'
                            date_filter_direction_text[(date_filter_direction>=112.5) & (date_filter_direction <157.5)] = 'SE'
                            for rect in np.arange(0,len(bar1)):
                                height = bar1[rect].get_height()
                                ax1.text(bar1[rect].get_x() + bar1[rect].get_width() / 2.0, height,date_filter_direction_text[rect], ha='center', va='bottom',fontsize=35,fontweight='bold')
 
                            date_form = DateFormatter("%H:%M")
                            ax1.xaxis.set_major_formatter(date_form)
                            
                                     
                         # Use MetPy to read the file
                        f = Level2File(obj.get()['Body'])
                        sweep = 0
                      # First item in ray is header, which has azimuth angle
                        az = np.array([ray[0].az_angle for ray in f.sweeps[sweep]])
         
                        ref_hdr = f.sweeps[sweep][0][4][b'REF'][0]
                        ref_range = np.arange(ref_hdr.num_gates) * ref_hdr.gate_width + ref_hdr.first_gate
                        ref = np.array([ray[4][b'REF'][1] for ray in f.sweeps[sweep]])
         
                        rho_hdr = f.sweeps[sweep][0][4][b'RHO'][0]
                        rho_range = (np.arange(rho_hdr.num_gates + 1) - 0.5) * rho_hdr.gate_width + rho_hdr.first_gate
                        rho = np.array([ray[4][b'RHO'][1] for ray in f.sweeps[sweep]])
         
                        phi_hdr = f.sweeps[sweep][0][4][b'PHI'][0]
                        phi_range = (np.arange(phi_hdr.num_gates + 1) - 0.5) * phi_hdr.gate_width + phi_hdr.first_gate
                        phi = np.array([ray[4][b'PHI'][1] for ray in f.sweeps[sweep]])
         
                        zdr_hdr = f.sweeps[sweep][0][4][b'ZDR'][0]
                        zdr_range = (np.arange(zdr_hdr.num_gates + 1) - 0.5) * zdr_hdr.gate_width + zdr_hdr.first_gate
                        zdr = np.array([ray[4][b'ZDR'][1] for ray in f.sweeps[sweep]])
         
                        rLAT = f.sweeps[0][0][1].lat
                        rLON = f.sweeps[0][0][1].lon
         
         
         
                        bot_left_lon = 360+rLON - 12
                        top_right_lon = 360+rLON + 12
                        bot_left_lat = rLAT - 12
                        top_right_lat = rLAT + 12
         
                          ######################################################################
                          # Plot the data
                          # -------------
                          #
                          # Use MetPy and Matplotlib to plot the data
                          #
                          # KMUX - SF
                          # KBBX - Oroville
                          # KSOX - Santa Ana
                          # Get the NWS reflectivity colortable from MetPy
                        ref_norm, ref_cmap = ctables.registry.get_with_steps('NWSReflectivity', 5, 5)
         
         
         
                      # this declares a recentered projection for Pacific areas
                        usemap_proj = ccrs.PlateCarree(central_longitude=180)
                        usemap_proj._threshold /= 20.  # to make greatcircle smooth
                        if ii == "paper" and ((dt.year in paper_dates.year) and (dt.month in paper_dates.month) and (dt.day in paper_dates.day) and (dt.hour in paper_dates.hour)):
                            ax = fig.add_subplot(5,1,1,projection=usemap_proj)
                        else:
                            ax = fig.add_subplot(2,5,1,projection=usemap_proj)
                        ax.axis('off')
                        ax.set_extent([bot_left_lon+8, top_right_lon-8, bot_left_lat+8, top_right_lat-8], crs=ccrs.PlateCarree())
                  
                        geodetic = ccrs.Geodetic()
                        plate_carree = ccrs.PlateCarree(central_longitude=180)
         
         
         
                        ax.coastlines(resolution="110m",linewidth=1)
                        ax.gridlines(linestyle='--',color='black')
                        ax.add_feature(cfeature.LAND)
                        ax.add_feature(cfeature.OCEAN,color="white")
                        ax.add_feature(cfeature.COASTLINE)
                        ax.add_feature(cfeature.BORDERS, linestyle=':', zorder=3)
                        ax.add_feature(cfeature.STATES, linestyle=':', zorder=3)
                        
                        ax.text(0.05, 0.05, '(a)',fontsize=35,fontweight ='bold', transform=ax.transAxes)
                        
                      
                        
                         # plot grid lines
                        gl = ax.gridlines(draw_labels=False, crs=ccrs.PlateCarree(), color='gray', linewidth=0.3)
                        gl.top_labels = False
                        gl.right_labels = False
                        gl.xlabel_style = {'size': 25}
                        gl.ylabel_style = {'size': 25}
         
                        var_range = ref_range
                        colors = ref_cmap
                        var_data = ref
                        lbl = "NWSReflectivity"
                      # Turn into an array, then mask
                        data = np.ma.array(var_data)
                        data[np.isnan(data)] = np.ma.masked
         
                       # Convert az, range to a lat/lon
                        g = Geod(ellps='clrk66') # This is the type of ellipse the earth is projected on.
                                                # There are other types of ellipses you can use,
                                                # but the differences will be small
                        center_lat = np.ones([len(az),len(ref_range)])*rLAT
                        center_lon = np.ones([len(az),len(var_range)])*(360+rLON)
                        az2D = np.ones_like(center_lat)*az[:,None]
                        rng2D = np.ones_like(center_lat)*np.transpose(var_range[:,None])*1000
                        lon,lat,back=g.fwd(center_lon,center_lat,az2D,rng2D)
         
                      # Convert az,range to x,y
                        xlocs = var_range * np.sin(np.deg2rad(az[:, np.newaxis]))
                        ylocs = var_range * np.cos(np.deg2rad(az[:, np.newaxis]))
         
                      # Define norm for reflectivity
                        norm = ref_norm if colors == ref_cmap else None
         
         
                      # Plot the data
                        a = ax.pcolormesh(360+lon,lat, data, cmap=colors, norm=norm,transform=ccrs.PlateCarree())
                        ax.plot(360+slon,slat,marker='o', color='red', markersize=15, transform=ccrs.PlateCarree())
         
                        
                        
                        if (ii == "paper" and ((dt.year == list(paper_dates.year)[3]) and (dt.month == list(paper_dates.month)[3]) and (dt.day == list(paper_dates.day)[3]) and (dt.hour == list(paper_dates.hour)[3]))):
                            cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
                            #divider = make_axes_locatable(ax)
                            #cbar_ax = divider.append_axes("right", size="5%", axes_class=maxes.Axes, pad=0.01)
                           # cbar_ax = fig.add_axes([0, 0])  # Adjust these values as needed
                            cbar = fig.colorbar(a, cax = cax,orientation='vertical')
                            cbar.set_label(label=lbl,size=35,labelpad=0.5,fontweight='bold')
                            cbar.ax.tick_params(labelsize=30)
                        if ii=="presentation":
                            #cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
                            divider = make_axes_locatable(ax)
                            cbar_ax = divider.append_axes("right", size="5%", axes_class=maxes.Axes, pad=0.01)
                           # cbar_ax = fig.add_axes([0, 0])  # Adjust these values as needed
                            cbar = fig.colorbar(a, cax = cbar_ax,orientation='vertical')
                            cbar.set_label(label=lbl,size=35,labelpad=0.5,fontweight='bold')
                            cbar.ax.tick_params(labelsize=30)
                        #plt.setp(ax.get_xticklabels(), fontsize=40)
                        #plt.setp(ax.get_yticklabels(), fontsize=40)
                        #ax.tick_params(axis='x', labelsize=40)
                        #ax.tick_params(axis='y', labelsize=40)
                        #add_timestamp(ax, f.dt, y=0.02, high_contrast=False)
                        ax.set_title(savestr.split("_")[0][4:8]+"-"+savestr.split("_")[0][8:10]+"-"+savestr.split("_")[0][10:]+ " " +savestr.split("_")[1][0:2]+":"+savestr.split("_")[1][2:4], fontsize=45,pad=10,fontweight='bold')
                     
         
         
         
                    
                      # define lat, lon region of data for plotting
                        latmin, latmax = (15.5,65.5)
                        lonmin, lonmax = (-170.25,-105.75)
         
                        if metvar == 'IVT':
                            filepath = "G:/NCFR Thesis/NCFR_Thesis/data/MERRA2_400.tavg1_2d_int_Nx."+str(dt.year)+str("{:02d}".format(dt.month))+str("{:02d}".format(dt.day))+".SUB.nc"#G:\\NCFR Thesis\\NCFR_Thesis\\MERRA2_400.tavg1_2d_int_Nx."+dt.20170207.SUB.nc"
                            filepath2 = "G:/NCFR Thesis/NCFR_Thesis/data/MERRA2_400.tavg1_2d_slv_Nx."+str(dt.year)+str("{:02d}".format(dt.month))+str("{:02d}".format(dt.day))+".nc4"
                        elif metvar == '850TAdv' or metvar == '850T' or metvar == 'SLP':
                            filepath = "G:/NCFR Thesis/NCFR_Thesis/data/MERRA2_400.tavg1_2d_slv_Nx."+str(dt.year)+str("{:02d}".format(dt.month))+str("{:02d}".format(dt.day))+".nc4"
                       
                      #COLLECT VARIABLE DATA FROM MERRA2 FILE
                        merravar = {'Z500':'H','SLP':'SLP','850T':'T','Z850':'H'}
                      #open the netcdf file in read mode
                        gridfile = nc.Dataset(filepath,mode='r')
                        gridlat = gridfile.variables['lat'][:]
                        gridlon = gridfile.variables['lon'][:]
                        if metvar == 'IVT':
                            gridfile2 = nc.Dataset(filepath2,mode='r')
                            Uvapor = gridfile.variables['UFLXQV'][:]
                            Vvapor = gridfile.variables['VFLXQV'][:]
                            merra = np.sqrt(Uvapor**2 + Vvapor**2)
                            
                            
                            
                            merra_upper_windU = gridfile2.variables['U250'][:]
                            merra_upper_windV = gridfile2.variables['V250'][:]
                            
                            windspeed_upper = (merra_upper_windU ** 2 + merra_upper_windV**2)**0.5
                            merra_slp = gridfile2.variables['SLP'][:]/100
                            
                            
                            gridlat2 = gridfile2.variables['lat'][:]
                            gridlon2 = gridfile2.variables['lon'][:]
                        elif metvar == '850TAdv': #temperature advection
                            gftadv = xarray.open_dataset(filepath)
                            UT = gftadv.U850
                            VT = gftadv.V850
                            T = gftadv.T850
                            #dx, dy = mpcalc.lat_lon_grid_deltas(gridlon, gridlat)
                            # proj=ccrs.LambertConformal(central_longitude=-90)
                            # lon,lat=np.meshgrid(gridlon,gridlat)
                            # output=proj.transform_points(ccrs.PlateCarree(central_longitude=180),lon,lat)
                            # x,y=output[:,:,0],output[:,:,1]
                            # gradx=np.gradient(x,axis=1)
                            # grady=np.gradient(y,axis=0)
                            # merra=-(UT*(np.gradient(T,axis=1)/gradx)+VT*(np.gradient(T,axis=0)/grady))*3600
                           
                            merra = np.array(mpcalc.advection(T, UT, VT)) * 3600
                           
                             #merra = mpcalc.advection(T, u=UT, v=VT)
                        elif metvar == 'SLP': 
                            merra = gridfile.variables['SLP'][:]/100
                            
                        elif metvar == '850T':
                            merra = gridfile.variables['T850'][:]
                        
                        
                        
                        if metvar != 'IVT' and metvar != '850TAdv':
                            uwind = gridfile.variables['U10M'][:]
                            vwind = gridfile.variables['V10M'][:]
                            
                            windspeed = (uwind ** 2 + vwind**2)**0.5
                        
                        
                            latlims = np.logical_and(gridlat2 > latmin, gridlat2 < latmax)
                            latind = np.where(latlims)[0]
                            gridlatreduced = gridlat[latind]
                          #reduce lon
                            lonlims = np.logical_and(gridlon2 > lonmin, gridlon2 < lonmax)
                            lonind = np.where(lonlims)[0]
                            gridlonreduced = gridlon2[lonind]
                          #reduce pressure
                          
                            merrareducedu = uwind[:,latind,:]
                            merrareduced_u = merrareducedu[:,:,lonind]
                            merrareducedv = vwind[:,latind,:]
                            merrareduced_v = merrareducedv[:,:,lonind]
                            merrareducedw = windspeed[:,latind,:]
                            merrareduced_w = merrareducedw[:,:,lonind]
                      
                            merrareduced_u = np.squeeze(merrareduced_u)
                            merrareduced_v = np.squeeze(merrareduced_v)
                            merrareduced_w = np.squeeze(merrareduced_w)
                            
                            merrareduced_u_h = merrareduced_u[dt.hour,:,:]
                            merrareduced_v_h = merrareduced_v[dt.hour,:,:]
                            
                            
                        
         
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
                        
                        merra = np.squeeze(merra)
                        
                        ##################################SLP#################################
                        #reduce lat
                        if metvar == 'IVT':
                            latlims = np.logical_and(gridlat2 > latmin, gridlat2 < latmax)
                            latind = np.where(latlims)[0]
                            gridlatreduced = gridlat2[latind]
                            #reduce lon
                            lonlims = np.logical_and(gridlon2 > lonmin, gridlon2 < lonmax)
                            lonind = np.where(lonlims)[0]
                            gridlonreduced = gridlon2[lonind]
                            
                            
                            
                            latlims = np.logical_and(gridlat > latmin, gridlat < latmax)
                            latind1 = np.where(latlims)[0]
                            #gridlatreduced = gridlat[latind]
                            #reduce lon
                            lonlims = np.logical_and(gridlon > lonmin, gridlon < lonmax)
                            lonind1 = np.where(lonlims)[0]
                            #gridlonreduced = gridlon[lonind]
                              
                            Uvapor = Uvapor[:,latind1,:]
                            Uvapor = Uvapor[:,:,lonind1]
                            Uvapor_hour = Uvapor[dt.hour,:,:]
                            mupperU = merra_upper_windU[:,latind,:]
                            mupperU = mupperU[:,:,lonind]
                            mupperU_hour = mupperU[dt.hour,:,:]
                            Vvapor = Vvapor[:,latind1,:]
                            Vvapor = Vvapor[:,:,lonind1]
                            Vvapor_hour = Vvapor[dt.hour,:,:]
                            mupperV = merra_upper_windV[:,latind,:]
                            mupperV = mupperV[:,:,lonind]
                            mupperV_hour = mupperV[dt.hour,:,:]
                            
                            
                            merrareduced_slp = merra_slp[:,latind,:]
                            merrareduced_slp = merrareduced_slp[:,:,lonind]
                            arr_slp = merrareduced_slp[dt.hour,:,:]
                        
                        #isolate to hour
                        
                        arr = merrareduced[dt.hour,:,:]
                        
         
         
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
                        mini = round(np.nanmin(arr),5)#0.002 *3600 #s to hour
                        maxi = round(np.nanmax(arr),5)
                        # if metvar=="850TAdv":
                        #     print(np.nanmin(arr))
                        #     print(np.nanmax(arr))
                        lowanom, highanom = (mini, maxi)
                        #newmap = center_colormap(lowanom, highanom, center=0)
                        lowlims = {'Z500':2850,'SLP':975,'IVT':0,'300W':0,'850T':240,'Z500Anom':lowanom,'Z850':1187,'SLPAnom':lowanom,'850TAdv':-10}
                        highlims = {'Z500':5700,'SLP':1050,'IVT':1700,'300W':56,'850T':293,'Z500Anom':highanom,'Z850':1548,'SLPAnom':highanom,'850TAdv':10}
         
                        contourstart = {'Z500':3000,'SLP':975,'IVT':250,'300W':5,'850T':240,'Z500Anom':-1.75,'Z850':1190,'SLPAnom':-2.25,'850TAdv':-10}
                        contourint = {'Z500':200,'SLP':4,'IVT':100,'300W':5,'850T':2.5,'Z500Anom':0.25,'Z850':30,'SLPAnom':0.25,'850TAdv':1.5}
         
                        cbarstart = {'Z500':3000,'SLP':975,'IVT':250,'300W':0,'850T':240,'Z500Anom':-2.0,'Z850':1200,'SLPAnom':-2.4,'850TAdv':-10}
                        cbarint = {'Z500':500,'SLP':5,'IVT':150,'300W':10,'850T':5,'Z500Anom':0.5,'Z850':50,'SLPAnom':0.4,'850TAdv':1.5}
         
                        colormap = {'Z500':'jet','SLP':'rainbow','IVT':'gnuplot2_r','300W':'hot_r','850T':'turbo','Z500Anom':'turbo','Z850':'turbo','SLPAnom':'turbo','850TAdv':'coolwarm'}
                        cbarlabs = {'Z500':'m','SLP':'hPa','IVT':'kg $\mathregular{m^{-1}}$ $\mathregular{s^{-1}}$','300W':'m/s','850T':'K','Z500Anom':r'$\mathbf{\sigma}$','Z850':'m','SLPAnom':r'$\mathbf{\sigma}$','850TAdv':'Degrees/hr'}
                        plottitle = {'Z500':'Z500','SLP':'SLP','IVT':'IVT','300W':'300 hPa Wind','850T':'850 hPa Temperature','Z500Anom':'Z500 Anomaly','Z850':'Z850','SLPAnom':'SLP Anomaly','850TAdv':'850 hPa \n Temperature Advection'}
                      #%% PLOT NODES from MATLAB
         
                      #create subplot for mapping multiple timesteps
                      #MAP DESIRED VARIABLE
                      # define date of plot
         
                      # Read in NetCDF4 file (add a directory path if necessary):
         
                      #Start Plotting Data
         
                      # Plot the data using matplotlib and cartopy
                     
         
                        datetitle =  plottitle[metvar]#+" on " + str(dt.year)+"-"+str(dt.month)+"-"+str(dt.day)+ "-" + str(dt.hour) +":00 UTC"
                        figtitle = plottitle[metvar]+" on " + str(dt.year)+"-"+str(dt.month)+"-"+str(dt.day)+ "-" + str(dt.hour)
         
         
         
         
                        usemap_proj = ccrs.PlateCarree(central_longitude=180)
                        usemap_proj._threshold /= 20.  # to make greatcircle smooth
                        
                        
                        if ii == "paper" and ((dt.year in paper_dates.year) and (dt.month in paper_dates.month) and (dt.day in paper_dates.day) and (dt.hour in paper_dates.hour)):
                            ax2 = fig.add_subplot(5,1,2+var,projection=usemap_proj)
                        else:
                            ax2 = fig.add_subplot(2,5,2+var,projection=usemap_proj)
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
                        
                        ax2.text(0.05, 0.05, pt,fontsize=35,fontweight ='bold', transform=ax2.transAxes)
                        
                        gl = ax2.gridlines(draw_labels=False, crs=ccrs.PlateCarree(), color='gray', linewidth=0.3)
                        gl.top_labels = False
                        gl.right_labels = False
                        gl.xlabel_style = {'size': 15}
                        gl.ylabel_style = {'size': 15}
                      # set appropriate extents: (lon_min, lon_max, lat_min, lat_max)
                       
                        lon, lat = np.meshgrid(gridlonreduced,gridlatreduced)
                        if metvar=="850TAdv" or metvar=="SLP" or metvar == "850T":
                            ax2.set_extent([bot_left_lon, top_right_lon, bot_left_lat, top_right_lat], crs=ccrs.PlateCarree())
                        else:
                            ax2.set_extent([lonmin, lonmax, latmin, latmax], crs=ccrs.PlateCarree())
                        ax2.plot(360+slon,slat,marker='o', color='red', markersize=15, transform=ccrs.PlateCarree())
                        
                        
                        
                      #define area threshold for basemap
                        area_thresh = 1E4
                      #create equidistant cylindrical projection basemap
         
                      #Set contour levels, then draw the plot and a colorbar
                        
                        contour_c = '0.1'
                        contour_w = 0.7
                        
                        if metvar != "IVT":
                            mp =  ax2.contourf(lon, lat, arr, np.arange(contourstart[metvar],highlims[metvar],contourint[metvar]), transform=ccrs.PlateCarree(),cmap=colormap[metvar])
                        
                        #if metvar=="850TAdv":
                        #    mp =  ax2.contourf(lon, lat, arr, np.linspace(-4,4,5), transform=ccrs.PlateCarree(),cmap=colormap[metvar])
                            
                        #    cbar = plt.colorbar(mp,ticks=np.linspace(np.nanmin(arr),np.nanmax(arr),5),orientation='vertical',pad=0.08)
                        #else:
                        if (ii == "paper" and ((dt.year == list(paper_dates.year)[3]) and (dt.month == list(paper_dates.month)[3]) and (dt.day == list(paper_dates.day)[3]) and (dt.hour == list(paper_dates.hour)[3]))) and metvar !="IVT":
                            cax = fig.add_axes([ax2.get_position().x1+0.01,ax2.get_position().y0,0.02,ax2.get_position().height])
                            #divider = make_axes_locatable(ax2)
                            #cbar_ax = divider.append_axes("right", size="5%", axes_class=maxes.Axes, pad=0.01)
                            cbar = plt.colorbar(mp,ticks=np.arange(cbarstart[metvar],highlims[metvar],cbarint[metvar]),cax=cax,orientation='vertical')
                            cbar.set_label(cbarlabs[metvar],fontsize=35,labelpad=1,fontweight='bold')
                            cbar.ax.tick_params(labelsize=30)
                        if ii == "presentation" and metvar != "IVT":
                            #cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
                            divider = make_axes_locatable(ax2)
                            cbar_ax = divider.append_axes("right", size="5%", axes_class=maxes.Axes, pad=0.01)
                           # cbar_ax = fig.add_axes([0, 0])  # Adjust these values as needed
                            cbar = plt.colorbar(mp,ticks=np.arange(cbarstart[metvar],highlims[metvar],cbarint[metvar]),cax=cbar_ax,orientation='vertical')
                            cbar.set_label(label=cbarlabs[metvar],size=35,labelpad=0.5,fontweight='bold')
                            cbar.ax.tick_params(labelsize=30)
                        if ii == "paper":    
                            ax2.text(.5,.95,datetitle,
                                       horizontalalignment='center',
                                       fontweight="bold",
                                       transform=ax2.transAxes,fontsize=38, bbox=dict(facecolor='white', edgecolor='black', pad=10.0))
                        else:
                            ax2.text(.5,.95,datetitle,
                                       horizontalalignment='center',
                                       fontweight="bold",
                                       transform=ax2.transAxes,fontsize=34, bbox=dict(facecolor='white', edgecolor='black', pad=10.0))
                        
                        if metvar == "IVT":
                            windspeed_upper = (mupperU_hour ** 2 + mupperV_hour**2)**0.5
                            if  ii == "paper":
                                mp =  ax2.contourf(lon, lat, windspeed_upper, np.arange(0,120,10), transform=ccrs.PlateCarree(),cmap='gray_r')
                                mp2 =  ax2.contour(lon, lat, arr, np.arange(contourstart[metvar],highlims[metvar],contourint[metvar]), transform=ccrs.PlateCarree(),cmap=colormap[metvar])
                            if ((ii == "paper" and ((dt.year == list(paper_dates.year)[3]) and (dt.month == list(paper_dates.month)[3]) and (dt.day == list(paper_dates.day)[3]) and (dt.hour == list(paper_dates.hour)[3])))) :
                            
                                
                                
                                mp =  ax2.contourf(lon, lat, windspeed_upper, np.arange(0,120,10), transform=ccrs.PlateCarree(),cmap='gray_r')
                                #if ii == "paper" and ((dt.year in paper_dates.year) and (dt.month in paper_dates.month) and (dt.day in paper_dates.day) and (dt.hour in paper_dates.hour)):
                                divider = make_axes_locatable(ax2)
                                cbar_ax = divider.append_axes("right", size="5%", axes_class=maxes.Axes, pad=0.1)
                               # cbar_ax = fig.add_axes([0, 0])  # Adjust these values as needed
                                cbar = plt.colorbar(mp,ticks=np.arange(0,120,10),cax=cbar_ax,orientation='vertical')
                                cbar.set_label(label='250 hPa Winds (m/s)',size=35,labelpad=0.7,fontweight='bold')
                                cbar.ax.tick_params(labelsize=30)
                                #ax2.quiver(lon,lat,Uvapor_hour,Vvapor_hour,transform=ccrs.PlateCarree(),regrid_shape=20,color="lightblue")
                                mp2 =  ax2.contour(lon, lat, arr, np.arange(contourstart[metvar],highlims[metvar],contourint[metvar]), transform=ccrs.PlateCarree(),cmap=colormap[metvar],alpha=0.7)
                                divider2 = make_axes_locatable(ax2)
                                cbar_ax2 = divider2.append_axes("bottom", size="5%", axes_class=maxes.Axes, pad=0.01)
                               # cbar_ax = fig.add_axes([0, 0])  # Adjust these values as needed
                                cbar2 = plt.colorbar(mp2,ticks=np.arange(contourstart[metvar],highlims[metvar],contourint[metvar]),cax=cbar_ax2,orientation='horizontal')
                                cbar2.set_label(label=cbarlabs[metvar],size=35,labelpad=0.5,fontweight='bold')
                                cbar2.ax.tick_params(labelsize=19)
                            if (ii == "presentation") :
                            
                                
                                
                                mp =  ax2.contourf(lon, lat, windspeed_upper, np.arange(0,120,10), transform=ccrs.PlateCarree(),cmap='gray_r')
                                #if ii == "paper" and ((dt.year in paper_dates.year) and (dt.month in paper_dates.month) and (dt.day in paper_dates.day) and (dt.hour in paper_dates.hour)):
                                divider = make_axes_locatable(ax2)
                                cbar_ax = divider.append_axes("right", size="5%", axes_class=maxes.Axes, pad=0.1)
                               # cbar_ax = fig.add_axes([0, 0])  # Adjust these values as needed
                                cbar = plt.colorbar(mp,ticks=np.arange(0,120,10),cax=cbar_ax,orientation='vertical')
                                cbar.set_label(label='250 hPa Winds (m/s)',size=35,labelpad=0.7,fontweight='bold')
                                cbar.ax.tick_params(labelsize=30)
                                #ax2.quiver(lon,lat,Uvapor_hour,Vvapor_hour,transform=ccrs.PlateCarree(),regrid_shape=20,color="lightblue")
                                mp2 =  ax2.contour(lon, lat, arr, np.arange(contourstart[metvar],highlims[metvar],contourint[metvar]), transform=ccrs.PlateCarree(),cmap=colormap[metvar])
                                divider2 = make_axes_locatable(ax2)
                                cbar_ax2 = divider2.append_axes("bottom", size="5%", axes_class=maxes.Axes, pad=0.01)
                               # cbar_ax = fig.add_axes([0, 0])  # Adjust these values as needed
                                cbar2 = plt.colorbar(mp2,ticks=np.arange(contourstart[metvar],highlims[metvar],contourint[metvar]),cax=cbar_ax2,orientation='horizontal')
                                cbar2.set_label(label=cbarlabs[metvar],size=35,labelpad=0.5,fontweight='bold')
                                cbar2.ax.tick_params(labelsize=16)
                             #ax2.quiver(lon,lat,mupperU_hour,mupperV_hour,transform=ccrs.PlateCarree(),regrid_shape=20,color="gray") 
                        if metvar == "850T":
                            #if ii == "paper" and ((dt.year in paper_dates.year) and (dt.month in paper_dates.month) and (dt.day in paper_dates.day) and (dt.hour in paper_dates.hour)):
                            ax2.quiver(lon,lat,merrareduced_u_h,merrareduced_v_h,transform=ccrs.PlateCarree(),regrid_shape=20) 
                            # else:
                            #     ax2.quiver(lon,lat,merrareduced_u_h,merrareduced_v_h,transform=ccrs.PlateCarree(),regrid_shape=20)
                         # else:
                           #     ax2.quiver(lon,lat,merrareduced_u_h,merrareduced_v_h,transform=ccrs.PlateCarree(),regrid_shape=20)
                      
                        var=var+1
                         # if ii == "paper" and ((dt.year in paper_dates.year) and (dt.month in paper_dates.month) and (dt.day in paper_dates.day) and (dt.hour in paper_dates.hour)):
                         #     #fig.suptitle('Precipitation Pulse ' + str(dt.year)+"-"+str(dt.month)+"-"+str(dt.day)+ "-" + str(dt.hour) +":00 UTC" , fontsize=60,fontweight="bold",y=0.98)
                         # else:
                         #     fig.suptitle('Precipitation Pulse Tracker', fontsize=60,fontweight="bold",y=0.98)
                            
                        
                    #plt.show()
                    if ii == "paper" and ((dt.year in paper_dates.year) and (dt.month in paper_dates.month) and (dt.day in paper_dates.day) and (dt.hour in paper_dates.hour)):
                         #fig.savefig(outdir+savestr+"_"+ii+".png")
                         for ax in fig.axes:
                             ax.set_anchor('C')
                        # fig.subplots_adjust(bottom=0, top=1, left=0, right=1, hspace=0, wspace=0)
                         if not os.path.exists(outdir+'paper\\'):
                              os.mkdir(outdir+'paper\\')
                         # with open(outdir+'data\\'+savestr+"_"+ii+".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                         #     #pickle.dump(fig, f) 
                         #     pickle.dump(fig,f)
                         fig.savefig(outdir+'paper\\'+savestr+"_"+ii+".jpg",bbox_inches='tight')
                    p_flag=p_flag+1             
                        
                    if ii == "presentation":
                        fig.savefig(outdir+savestr+"_"+ii+".png")
                        #png_dir = outdir+savestr+"_"+ii+".png"
                        #images.append(imageio.imread(png_dir))
                       
                        
                        # Make it pause at the end so that the viewers can ponder
                        # for _ in range(10):
                        #     images.append(imageio.imread(file_path))

        
                        
                    plt.close('all')
                    p_flag = p_flag+1
                    
    #imageio.mimsave(outdir+savestr+"movie_"+ii+".gif", images)





