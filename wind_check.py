fig = plt.figure(figsize=(20, 20))

ax2 = fig.add_subplot(1,1,1,projection=usemap_proj)
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

ax2.set_extent([bot_left_lon, top_right_lon, bot_left_lat, top_right_lat], crs=ccrs.PlateCarree())
 
mp =  ax2.contourf(wlon_filter,wlat_filter, wind, np.arange(np.min(wind),np.max(wind),2.5), transform=ccrs.PlateCarree(),cmap=colormap['SLP'])

ax2.quiver(wlon_filter,wlat_filter,windtotu,windtotv,transform=ccrs.PlateCarree(),regrid_shape=20) #sizes=dict(emptybarb=0.001, spacing=0.1, height=0.5)
mp2 = ax2.contour(lon,lat,arr,colors=contour_c,linewidths=contour_w,levels=np.arange(contourstart[metvar],highlims[metvar]+1,contourint[metvar]))

cbar = plt.colorbar(mp,ticks=np.arange(cbarstart[metvar],highlims[metvar],cbarint[metvar]),orientation='vertical',pad=0.08)
cbar.set_label(cbarlabs[metvar],fontsize=20,labelpad=0.5,fontweight='bold')
cbar.ax.tick_params(labelsize=20)