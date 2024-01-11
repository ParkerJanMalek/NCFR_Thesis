import urllib
import json
import pandas as pd
import matplotlib.pyplot as plt
import datetime  as dtetme
from datetime import datetime
import scipy.io
import numpy as np
import pickle
from matplotlib.dates import DayLocator, DateFormatter
import shutil
import os
import cartopy
import cartopy.crs as ccrs


#radar


######


import MRMS_data_pull as MRMS
import Nexrad_S3_Demo as radar
import MERRAII_IVT as merra


# stage four quanitiative precipitation (hourly radar product)
# find the ARs that are colocated
# Specify request parameters (as strings)
token = 'aa31874e86fb42d9b2ea6b293f1bb004' # use your own token

################################FUNCTIONS############################################

def count_and_sum_events(time_series,hours_between):
    event_start = []
    event_end = []
    event_acum = []
    event_id = []
    for i in np.arange(0,len(time_series)-1):
        if time_series[i] != 0 and ((np.all(time_series[i-hours_between:i] == 0) and not time_series[i-hours_between:i].empty) or 
                                    (i < hours_between and np.all(time_series[:i] == 0))):
           
            event_start.append(time_series.index[i])
            start_i = i
            j = i
        if time_series[i] != 0 and np.all(time_series[i+1:i+(hours_between+1)] == 0):
            event_end.append(time_series.index[i])
            end_i = i
            event_total = np.sum(time_series[start_i:end_i+1])
            event_acum.append(np.sum(time_series[start_i:end_i+1]))
            k = i
    event_id.append({'start':np.array(event_start),'end':np.array(event_end),'accum':np.array(event_acum)})     
                
    return(event_id)

def find_nearest_value(target, values):
    nearest = None
    min_difference = float('inf')

    for value in values:
        difference = abs(value - target)
        if difference < min_difference:
            min_difference = difference
            nearest = value

    return nearest

def myround(x, prec=2, base=.5):
  return round(base * round(float(x)/base),prec)

def data_availability_check(station_name,token,var):
    if responseDict_por['SUMMARY']['RESPONSE_MESSAGE'] == 'OK':
        return(True)
    else:
        return(False)
    
####################################################################################

station_name_list =['kove']

var = 'precip_accum_one_hour'
unit = 'precip|mm'
plot_full_record = False
output_5_perc_data = 'G:\\NCFR Thesis\\NCFR_Thesis\\precip_threshold_dates\\'

#define cumsum dates
start_date = dtetme.datetime(2017,2,7,0)
end_date = dtetme.datetime(2017,2,7,23)



for i in station_name_list:
    print(i)
    from datetime import datetime
    #bound = station_window[i]
    station_name = i
    args_por = {
        'obtimezone':'UTC',
        'vars':var,
        'stids': station_name,
        'units':unit,
        'token':token
    }

    apiString_por = urllib.parse.urlencode(args_por)
    url_por = "https://api.synopticdata.com/v2/stations/nearesttime"
    fullUrl_por = '{}?{}'.format(url_por,apiString_por)
    response_por = urllib.request.urlopen(fullUrl_por)
    responseDict_por = json.loads(response_por.read())
    if data_availability_check(station_name,token,var):
        por_start = responseDict_por['STATION'][0]['PERIOD_OF_RECORD']['start']
        por_end = responseDict_por['STATION'][0]['PERIOD_OF_RECORD']['end']
        station_lon = float(responseDict_por['STATION'][0]['LONGITUDE'])
        station_lat = float(responseDict_por['STATION'][0]['LATITUDE'])
        
        s = datetime.strptime(por_start[0:10],"%Y-%m-%d")
        e = datetime.strptime(por_end[0:10],"%Y-%m-%d")
        por_start_fmt = datetime.strftime(datetime.strptime(por_start[0:10],"%Y-%m-%d"),"%Y%m%d%H%S")
        por_end_fmt = datetime.strftime(datetime.strptime(por_end[0:10],"%Y-%m-%d"),"%Y%m%d%H%S")
        print(e.year-s.year)
        args = {
            'start':por_start_fmt,
            'end':por_end_fmt,
            'obtimezone':'UTC',
            'vars': var,
            'stids': station_name,
                #['C3BCC','C3BVS','C3CAT','C3DLA','C3DRW','C3FRC','C3GPO','C3HDC','C3HRD','C3NBB','C3NCM','C3POR','C3PVN','C3SKI','C3SKY','C3SOD','C3WDG','C3WPO'],
            'units':unit,
            'token':token
        }
        
        
        apiString = urllib.parse.urlencode(args)
        url = "https://api.synopticdata.com/v2/stations/timeseries"
        fullUrl = '{}?{}'.format(url,apiString)
        
        # Open the URL and convert the returned JSON into a dictionary
        response = urllib.request.urlopen(fullUrl)
        responseDict = json.loads(response.read())
        
        # Isolate the time and temperature from the response dictionary
        dateTime = responseDict['STATION'][0]['OBSERVATIONS']['date_time']
        value_pull = responseDict['STATION'][0]['OBSERVATIONS'][args['vars']+'_set_1']
        station_data_out = pd.Series(value_pull,index=pd.to_datetime(dateTime))
        station_data_out.to_csv('raw_rainfalls_'+station_name+'.csv')
        
        # resample to only the hourly observations. For inter-hourly observations, the mean of these values is taken
        station_data_hourly = station_data_out.resample('60min').mean()
        station_data_hourly.to_csv('hourly_rainfalls_'+station_name+'.csv')
        station_data_hourly[station_data_hourly.isnull()] = 0
        #station_data_hourly[station_data_hourly <=0.02500] = 0
        
        #isolate events to the 
        s_e_sum = count_and_sum_events(station_data_hourly,12)
        
        sum_events = s_e_sum[0]['accum']
        
        uq_start = s_e_sum[0]['start']
        uq_end = s_e_sum[0]['end']
        
        time_spans = []
        
        for i in np.arange(0,len(sum_events)-1):
            time_spans.append({'start':uq_start[i],'end':uq_end[i],'event_total':sum_events[i]})
        
        event_id = []
        event_total = []
        
        
        for r in time_spans:
            event_percentage = 100 * station_data_hourly[r['start']:r['end']]/r['event_total']
            dates_of_event = pd.DataFrame(np.unique(pd.to_datetime(station_data_hourly[r['start']:r['end']].index.date)))
            event_id.append({'start':r['start'],'end':r['end'],'event_rainfall':station_data_hourly[r['start']:r['end']],'event_total':r['event_total'] ,'event_perc':event_percentage,'multimodal':False})
            event_total.append(r['event_total'])

        
        #select Oroville events with greater than 3mm in any hour
        ts_selected  = [x for x in event_id if x['start'].date().year == 2017 and x['start'].date().month == 2] #selecte all events in February 2017
            
            
        ##############################cumulative precip sums (not used but saved)################################
        cum_sum_precip = station_data_hourly.groupby([station_data_hourly.index.year,
                                                 station_data_hourly.index.month,
                                                 station_data_hourly.index.day]).cumsum()
        
        
        cum_sum_precip_t = cum_sum_precip.reset_index()
        
        cols= ["index"]
        cum_sum_precip_t['Date'] = cum_sum_precip_t[cols].apply(lambda x: '-'.join(x.values.astype(str)), axis="columns")
        cum_sum_precip_t['Date'] = pd.to_datetime(cum_sum_precip_t["Date"])
        
        cum_precip = station_data_hourly.groupby([station_data_hourly.index.year,
                                                 station_data_hourly.index.month,
                                                 station_data_hourly.index.day]).sum()
        
        
        t=cum_sum_precip_t

        t=t.reset_index()


        ar_cols = ['year','month','day']
        
        merged_cumsum = t
        ###############################################################################################################

        ###########################################Plotting############################################################
        
        file_var2 = 'GaugeCorr_QPE_01H'
        #outdirck = 'G:\\NCFR Thesis\\NCFR_Thesis\\'+station_name+'_'+file_var2 + '_'+str(start_date.year)+str(start_date.month)+str(start_date.day)+ str(end_date.hour)+'_'+ str(end_date.year)+ str(end_date.month)+ str(end_date.day)+ str(end_date.hour) +'\\'
        # if not os.path.exists(outdirck):
        #     MRMS.map_event(start_date, end_date, station_lon, station_lat,station_name)
            
        multimodal= 1
        ts_selected[1]['start'] = pd.Timestamp('2017-2-6-22',tz='UTC')
        ts_selected_isolate = ts_selected[0:5]
        tt=[ts_selected_isolate[1]]
        for i in ts_selected_isolate:
            start_date = i['start']
            end_date = i['end']
            outdirck = 'G:\\NCFR Thesis\\NCFR_Thesis\\IVT_'+station_name + '_'+str(start_date.year)+str(start_date.month)+str(start_date.day)+ str(end_date.hour)+'_'+ str(end_date.year)+ str(end_date.month)+ str(end_date.day)+ str(end_date.hour) +'\\'
            if os.path.exists(outdirck):
                shutil.rmtree(outdirck)
            os.mkdir(outdirck)
            merra.pull_merra(i['start'], i['end'],station_data_hourly,i,station_name,station_lon,station_lat)
            outdirck = 'G:\\NCFR Thesis\\NCFR_Thesis\\Radar_'+station_name + '_'+str(start_date.year)+str(start_date.month)+str(start_date.day)+ str(end_date.hour)+'_'+ str(end_date.year)+ str(end_date.month)+ str(end_date.day)+ str(end_date.hour) +'\\'
            if os.path.exists(outdirck):
                shutil.rmtree(outdirck)
            os.mkdir(outdirck)
            radar.pull_radar(i['start'], i['end'],station_data_hourly,i,station_name)
            MRMS.map_event(start_date, end_date, station_lon, station_lat,station_name,station_data_hourly,i)
            # print(i)
            # date_filter = station_data_hourly[i['start']:i['end']]
            # date_filter_cumsum = date_filter.cumsum()
            # fig,ax = plt.subplots(figsize=(40, 20))
            # ax.bar(date_filter.index,date_filter,width=0.01)
            # ax.set_title('Synoptic Rainfall Pulse Event at Oroville Municipal Airport' , fontsize=60)
            # plt.ylabel('Precipiation (mm)', fontsize=50)
            # plt.xlabel('Date', fontsize=50)
            # plt.xticks(fontsize=35,rotation=40)
            # plt.yticks(fontsize=35)
            # ax.grid()
            # date_form = DateFormatter("%m-%d-%Y-%H")
            # ax.xaxis.set_major_formatter(date_form)
            # fig.tight_layout()
            # fig.savefig(station_name + "_"+args['vars'] +str(multimodal)+ '.jpeg')  
            # multimodal = multimodal + 1
        # fig1,ax1 = plt.subplots(figsize=(20, 12.5))
        # ax1.plot(date_filter.index,date_filter)
        # fig1.suptitle('Hourly Precipitation at Oroville Airport ('+station_name.upper() + ") \n in February 2017", fontsize=30)
        # plt.ylabel('Precipiation (mm)', fontsize=25)
        # plt.xlabel('Date', fontsize=25,labelpad=-0.5)
        # plt.xticks(fontsize=15,rotation=45)
        # plt.yticks(fontsize=15)
        # ax1.grid()
        # date_form = DateFormatter("%m-%d-%Y-%H")
        # ax1.xaxis.set_major_formatter(date_form)
        # ax2 = fig1.gca()
        # ax2.set_ylim([0, None])
    else:
        print("No data available for " + station_name)
        
        