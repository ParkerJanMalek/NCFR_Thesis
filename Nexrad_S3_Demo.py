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

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pyproj import Geod #this is used to convert range/azimuth to lat/lon

import matplotlib.axes as maxes

def pull_radar(start_date1,end_date1):
    start_date = dtetme.datetime(start_date1.year,start_date1.month,start_date1.day,start_date1.hour)
    end_date = dtetme.datetime(end_date1.year,end_date1.month,end_date1.day,end_date1.hour)
    print(start_date)
    print(end_date)
    s3 = boto3.resource('s3', config=Config(signature_version=botocore.UNSIGNED,
                                            user_agent_extra='Resource'))
    bucket = s3.Bucket('noaa-nexrad-level2')
    
    for dt in rrule.rrule(rrule.HOURLY, dtstart=start_date, until=end_date):
        print(rrule.rrule(rrule.HOURLY, dtstart=start_date, until=end_date))
        month = str(dt.month).zfill(2)
        day = str(dt.day).zfill(2)
        hour = str(dt.hour).zfill(2)
        print(str(dt.year) + '/' + month + '/' + day + '/KBBX/KBBX'+ str(dt.year) + month + day + '_'+hour)
        # construct a list of days/hours to loop through.
        
        radar_object1 = bucket.objects.filter(Prefix=str(dt.year) + '/' + month + '/' + day + '/KBBX/KBBX'+ str(dt.year) + month + day + '_'+hour) 
       
        for obj in radar_object1:
            fig, ax = plt.subplots(figsize=(20, 20))
            # Plot the data!
            savestr = obj.key.split("/")[-1]
            print(savestr)
            
            
            # Use MetPy to read the file
            f = Level2File(obj.get()['Body'])
            
            
            ######################################################################
            # Subset Data
            # -----------
            #
            # With the file comes a lot of data, including multiple elevations and products.
            # In the next block, we'll pull out the specific data we want to plot.
            #
            
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
            
            
            ax.axis('off')
            #ax.set_xticks([])
            #ax.set_yticks([])
            # this declares a recentered projection for Pacific areas
            usemap_proj = ccrs.PlateCarree(central_longitude=180)
            usemap_proj._threshold /= 20.  # to make greatcircle smooth
             
            ax = plt.axes(projection=usemap_proj)
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
            ax.plot(360+rLON,rLAT,marker='o', color='red', markersize=20, transform=ccrs.PlateCarree())
             
            divider = make_axes_locatable(ax)
             
            cax = divider.append_axes("right", size="5%", axes_class=maxes.Axes, pad=0.05)
            cbar = fig.colorbar(a, cax=cax, orientation='vertical')
            cbar.set_label(label=lbl,size=35)
            cbar.ax.tick_params(labelsize=35)
            plt.setp(ax.get_xticklabels(), fontsize=35)
            plt.setp(ax.get_yticklabels(), fontsize=35)
            ax.tick_params(axis='x', labelsize=35)
            ax.tick_params(axis='y', labelsize=35)
            #ax.set_aspect('equal', 'datalim')
            # ax.set_xlim(-100, 100)
            # ax.set_ylim(-100, 100)
            add_timestamp(ax, f.dt, y=0.02, high_contrast=False)
            # plt.suptitle(savestr[0:4]+' Level 2 Data '+ savestr.split("_")[0][4:]+"_"+savestr.split("_")[1], fontsize=50)
            # plt.tight_layout()
            # plt.show()
            # fig.savefig(savestr+".png")
            # i=i+1
            ax.spines['right'].set_visible(False)
            plt.suptitle(savestr[0:4]+' Level 2 Data '+ savestr.split("_")[0][4:]+" "+savestr.split("_")[1][0:2]+":"+savestr.split("_")[1][2:4] +" UTC", fontsize=50)
            plt.tight_layout()
            plt.show()
           # return(fig)
            fig.savefig(savestr+".png")
             
            
        
    
         
              
