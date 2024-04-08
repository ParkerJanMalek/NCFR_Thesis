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
from metpy.io import Level2File
from metpy.plots import add_timestamp, ctables
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from dateutil import rrule
from datetime import datetime, timedelta
import datetime  as dtetme
from matplotlib.dates import DayLocator, DateFormatter

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pyproj import Geod #this is used to convert range/azimuth to lat/lon

import matplotlib.axes as maxes

def pull_radar(start_date1,end_date1,station_data,ts_selected,station_name,radar_name,slat,slon):
    start_date = dtetme.datetime(start_date1.year,start_date1.month,start_date1.day,start_date1.hour)
    end_date = dtetme.datetime(end_date1.year,end_date1.month,end_date1.day,end_date1.hour)
    print(start_date)
    print(end_date)
    s3 = boto3.resource('s3', config=Config(signature_version=botocore.UNSIGNED,
                                            user_agent_extra='Resource'))
    bucket = s3.Bucket('noaa-nexrad-level2')
    outdir = 'G:\\NCFR Thesis\\NCFR_Thesis\\combined_'+station_name + '_'+str(start_date.year)+str(start_date.month)+str(start_date.day)+ str(end_date.hour)+'_'+ str(end_date.year)+ str(end_date.month)+ str(end_date.day)+ str(end_date.hour) +'\\'
    for dt in rrule.rrule(rrule.HOURLY, dtstart=start_date, until=end_date):
        month = str(dt.month).zfill(2)
        day = str(dt.day).zfill(2)
        hour = str(dt.hour).zfill(2)
        print(str(dt.year) + '/' + month + '/' + day + '/'+radar_name+'/'+radar_name+ str(dt.year) + month + day + '_'+hour)
        # construct a list of days/hours to loop through.

        radar_object1 = bucket.objects.filter(Prefix=str(dt.year) + '/' + month + '/' + day + '/KBBX/KBBX'+ str(dt.year) + month + day + '_'+hour)
        #print(len(list(radar_object1)))
        if len(list(radar_object1)) == 0:
            radar_object1 = bucket.objects.filter(Prefix=str(dt.year) + '/' + month + '/' + day + '/KDAX/KDAX'+ str(dt.year) + month + day + '_'+hour)
        fig = plt.figure(figsize=(40, 20))
        for obj in radar_object1:


           #%% IMPORT MERRA2 DATA
           # define metvar
           metvars = ['SLP', '300W','Z500Anom','SLPAnom','Z850','850T','850TAdv']
           metvars = ['IVT','850T','SLP']#]
           #metvar = '300W'
           fig = plt.figure(figsize=(60, 20))
           var=0
           for metvar in metvars:
               savestr = obj.key.split("/")[-1]
               print(savestr)
               ax1 = fig.add_subplot(2,1,2)
               #fig, ax = plt.subplots(figsize=(20, 20))
               # Plot the data!
               date_filter = station_data[ts_selected['start']:ts_selected['end']]
               ax1.bar(date_filter.index,date_filter,width=0.01)
               ax1.axvline(x=dt,linewidth=4, color='r')
               plt.ylabel('Precipiation (mm)', fontsize=30)
               plt.xlabel('Hour', fontsize=30)
               plt.xticks(fontsize=30,rotation=40)
               plt.yticks(fontsize=30)
               plt.ylim([0,4.5])
               ax1.grid()
               date_form = DateFormatter("%H:%M")
               ax1.xaxis.set_major_formatter(date_form)
               #ax1.set_title('Precipitation Pulse Tracker',fontsize=40,pad=10)

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



               bot_left_lon = 360+rLON - 4
               top_right_lon = 360+rLON + 4
               bot_left_lat = rLAT - 4
               top_right_lat = rLAT + 4

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



             #ax.set_xticks([])
             #ax.set_yticks([])
             # this declares a recentered projection for Pacific areas
               usemap_proj = ccrs.PlateCarree(central_longitude=180)
               usemap_proj._threshold /= 20.  # to make greatcircle smooth
               ax = fig.add_subplot(2,4,1,projection=usemap_proj)
               ax.axis('off')
             #ax = plt.axes(projection=usemap_proj)
             # set appropriate extents: (lon_min, lon_max, lat_min, lat_max)
               ax.set_extent([bot_left_lon, top_right_lon, bot_left_lat, top_right_lat], crs=ccrs.PlateCarree())

               geodetic = ccrs.Geodetic()
               plate_carree = ccrs.PlateCarree(central_longitude=180)



               ax.add_feature(cfeature.LAND)
               ax.add_feature(cfeature.OCEAN,color="white")
               ax.add_feature(cfeature.COASTLINE)
               ax.add_feature(cfeature.BORDERS, linestyle=':', zorder=2)
               ax.add_feature(cfeature.STATES, linestyle=':', zorder=2)
                # plot grid lines
               gl = ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), color='gray', linewidth=0.3)
               gl.top_labels = False
               gl.right_labels = False
               gl.xlabel_style = {'size': 25}
               gl.ylabel_style = {'size': 25}
             # for var_data, colors, lbl in zip((ref),
             #                                                 ([ref_cmap]),
             #                                                 ("NWSReflectivity")):

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
               ax.plot(360+slon,slat,marker='o', color='red', markersize=20, transform=ccrs.PlateCarree())

               divider = make_axes_locatable(ax)

               cax = divider.append_axes("right", size="5%", axes_class=maxes.Axes, pad=1.05)
               cbar = fig.colorbar(a, cax=cax, orientation='vertical',pad=0.08)
               cbar.set_label(label=lbl,size=20,labelpad=0.5)
               cbar.ax.tick_params(labelsize=20)
               plt.setp(ax.get_xticklabels(), fontsize=20)
               plt.setp(ax.get_yticklabels(), fontsize=20)
               ax.tick_params(axis='x', labelsize=20)
               ax.tick_params(axis='y', labelsize=20)
               #ax.set_aspect('equal', 'datalim')
               # ax.set_xlim(-100, 100)
               # ax.set_ylim(-100, 100)
               add_timestamp(ax, f.dt, y=0.02, high_contrast=False)
               # plt.suptitle(savestr[0:4]+' Level 2 Data '+ savestr.split("_")[0][4:]+"_"+savestr.split("_")[1], fontsize=50)
               # plt.tight_layout()
               # plt.show()
               # fig.savefig(savestr+".png")
               # i=i+1
               ax.set_title('Radar Data on '+ savestr.split("_")[0][4:8]+"-"+savestr.split("_")[0][8:10]+"-"+savestr.split("_")[0][10:]+ " " +savestr.split("_")[1][0:2]+":"+savestr.split("_")[1][2:4], fontsize=40,pad=10)
               plt.tight_layout()
              # plt.show()
              # return(fig)



             #############################################################################
             #outdir2 = 'G:/NCFR Thesis/NCFR_Thesis/IVT_'+station_name + '_'+str(start_date.year)+str(start_date.month)+str(start_date.day)+ str(end_date.hour)+'_'+ str(end_date.year)+ str(end_date.month)+ str(end_date.day)+ str(end_date.hour) +'\\'
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


               if metvar == 'IVT':
                   filepath = "G:/NCFR Thesis/NCFR_Thesis/data/MERRA2_400.tavg1_2d_int_Nx."+str(dt.year)+str("{:02d}".format(dt.month))+str("{:02d}".format(dt.day))+".SUB.nc"#G:\\NCFR Thesis\\NCFR_Thesis\\MERRA2_400.tavg1_2d_int_Nx."+dt.20170207.SUB.nc"
               elif metvar == '850TAdv' or metvar == '850T' or metvar == 'SLP':
                   filepath = "G:/NCFR Thesis/NCFR_Thesis/data/MERRA2_400.tavg1_2d_slv_Nx."+str(dt.year)+str("{:02d}".format(dt.month))+str("{:02d}".format(dt.day))+".nc4"
               windv = "G:/NCFR Thesis/NCFR_Thesis/era5_10m_v_component_of_wind_2017_hourly_165E-80W_25N-80N.nc" 
               windu = "G:/NCFR Thesis/NCFR_Thesis/era5_10m_u_component_of_wind_2017_hourly_165E-80W_25N-80N.nc"
               wv = nc.Dataset(windv,mode='r')
               wu = nc.Dataset(windu,mode='r')
               
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
               elif metvar == '850TAdv': #temperature advection
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
               elif metvar == 'SLP':
                   merra = gridfile.variables['SLP'][:]/100
               elif metvar == '850T':
                   merra = gridfile.variables['T850'][:]
                    
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
               #isolate to hour
               arr = merrareduced[dt.hour,:,:]

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
               mini = round(np.min(arr),1)#0.002 *3600 #s to hour
               maxi = round(np.max(arr),1)
               lowanom, highanom = (mini, maxi)
               newmap = center_colormap(lowanom, highanom, center=0)
               lowlims = {'Z500':2850,'SLP':970,'IVT':0,'300W':0,'850T':252,'Z500Anom':lowanom,'Z850':1187,'SLPAnom':lowanom,'850TAdv':mini}
               highlims = {'Z500':5700,'SLP':1022,'IVT':1700,'300W':56,'850T':293,'Z500Anom':highanom,'Z850':1548,'SLPAnom':highanom,'850TAdv':maxi}

               contourstart = {'Z500':3000,'SLP':980,'IVT':0,'300W':5,'850T':250,'Z500Anom':-1.75,'Z850':1190,'SLPAnom':-2.25,'850TAdv':mini}
               contourint = {'Z500':200,'SLP':4,'IVT':100,'300W':5,'850T':2.5,'Z500Anom':0.25,'Z850':30,'SLPAnom':0.25,'850TAdv':maxi/6}

               cbarstart = {'Z500':3000,'SLP':980,'IVT':0,'300W':0,'850T':250,'Z500Anom':-2.0,'Z850':1200,'SLPAnom':-2.4,'850TAdv':mini}
               cbarint = {'Z500':500,'SLP':5,'IVT':150,'300W':10,'850T':5,'Z500Anom':0.5,'Z850':50,'SLPAnom':0.4,'850TAdv':maxi/6}

               colormap = {'Z500':'jet','SLP':'rainbow','IVT':'gnuplot2_r','300W':'hot_r','850T':'turbo','Z500Anom':newmap,'Z850':'turbo','SLPAnom':newmap,'850TAdv':'coolwarm'}
               cbarlabs = {'Z500':'m','SLP':'hPa','IVT':'kg $\mathregular{m^{-1}}$ $\mathregular{s^{-1}}$','300W':'m/s','850T':'K','Z500Anom':r'$\mathbf{\sigma}$','Z850':'m','SLPAnom':r'$\mathbf{\sigma}$','850TAdv':'Degrees/hr'}
               plottitle = {'Z500':'Z500','SLP':'SLP','IVT':'IVT','300W':'300 hPa Wind','850T':'850 hPa Temperature','Z500Anom':'Z500 Anomaly','Z850':'Z850','SLPAnom':'SLP Anomaly','850TAdv':'850 hPa Temperature Advection'}
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
             #fig.add_subplot



             # ax1 = fig.add_subplot(211)
             # #fig, ax = plt.subplots(figsize=(20, 20))
             # # Plot the data!
             # date_filter = station_data[ts_selected['start']:ts_selected['end']]
             # ax1.bar(date_filter.index,date_filter,width=0.01)
             # ax1.axvline(x=dt,linewidth=4, color='r')
             # plt.ylabel('Precipiation (mm)', fontsize=25)
             # plt.xticks(fontsize=15,rotation=40)
             # plt.yticks(fontsize=15)
             # ax1.grid()
             # date_form = DateFormatter("%m-%d-%Y-%H")
             # ax1.xaxis.set_major_formatter(date_form)


             # ax1.set_title(datetitle, fontsize=40)




               usemap_proj = ccrs.PlateCarree(central_longitude=180)
               usemap_proj._threshold /= 20.  # to make greatcircle smooth
  
               ax2 = fig.add_subplot(2,4,2+var,projection=usemap_proj)
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
               if metvar=="850TAdv" or metvar=="SLP" or metvar == "850T":
                   ax2.set_extent([bot_left_lon, top_right_lon, bot_left_lat, top_right_lat], crs=ccrs.PlateCarree())
               else:
                   ax2.set_extent([lonmin, lonmax, latmin, latmax], crs=ccrs.PlateCarree())
               ax2.plot(360+slon,slat,marker='o', color='red', markersize=15, transform=ccrs.PlateCarree())
  
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
               
               contour_c = '0.1'
               contour_w = 0.7
               lon, lat = np.meshgrid(gridlonreduced,gridlatreduced)
             #ax2.scatter(-120.9,39.5,color='r',marker='*',linewidths=5)
               colorm = ax2.pcolor(lon,lat,arr,shading='auto',cmap=colormap[metvar],vmin=lowlims[metvar],vmax=highlims[metvar])
               mp =  ax2.contourf(lon, lat, arr, np.arange(contourstart[metvar],highlims[metvar],contourint[metvar]), transform=ccrs.PlateCarree(),cmap=colormap[metvar])
               mp2 = ax2.contour(lon,lat,arr,colors=contour_c,linewidths=contour_w,levels=np.arange(contourstart[metvar],highlims[metvar]+1,contourint[metvar]))
  
               cbar = plt.colorbar(mp,ticks=np.arange(cbarstart[metvar],highlims[metvar],cbarint[metvar]),orientation='vertical',pad=0.08)
               cbar.set_label(cbarlabs[metvar],fontsize=20,labelpad=0.5,fontweight='bold')
               cbar.ax.tick_params(labelsize=20)
               ax2.set_title(figtitle,fontsize=40,pad=10)
               fig.tight_layout()
               var=var+1
               fig.suptitle('Precipitation Pulse Tracker', fontsize=60,fontweight="bold",y=0.98)
           #plt.show()
           fig.savefig(outdir+savestr+".png")
           plt.close('all')





