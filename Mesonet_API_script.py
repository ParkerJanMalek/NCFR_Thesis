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
import combine_figures as cb


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

plot_full_record = False
output_5_perc_data = 'G:\\NCFR Thesis\\NCFR_Thesis\\precip_threshold_dates\\'




#define cumsum dates
start_date = dtetme.datetime(2017,2,7,0)
end_date = dtetme.datetime(2017,2,7,23)

unit = ['degrees','temp|C','speed|mph','precip|mm']
var = ['wind_direction','air_temp','wind_speed','precip_accum_one_hour']

for i in station_name_list:
    print(i)
    from datetime import datetime
    #bound = station_window[i]
    
    station_name = i
    station_data_out = {}
    for j in np.arange(0,len(unit)):
        print(j)
        if unit[j] == 'degrees':
            
            args_por = {
                'obtimezone':'UTC',
                'vars':var[j],
                'stids': station_name,
                'token':token
            }
        else:
            args_por = {
                'obtimezone':'UTC',
                'vars':var[j],
                'stids': station_name,
                'units':unit[j],
                'token':token
            }
        apiString_por = urllib.parse.urlencode(args_por)
        url_por = "https://api.synopticdata.com/v2/stations/nearesttime"
        fullUrl_por = '{}?{}'.format(url_por,apiString_por)
        response_por = urllib.request.urlopen(fullUrl_por)
        responseDict_por = json.loads(response_por.read())
        por_start = responseDict_por['STATION'][0]['PERIOD_OF_RECORD']['start']
        por_end = responseDict_por['STATION'][0]['PERIOD_OF_RECORD']['end']
        station_lon = float(responseDict_por['STATION'][0]['LONGITUDE'])
        station_lat = float(responseDict_por['STATION'][0]['LATITUDE'])
        s = datetime.strptime(por_start[0:10],"%Y-%m-%d")
        e = datetime.strptime(por_end[0:10],"%Y-%m-%d")
        por_start_fmt = datetime.strftime(datetime.strptime('2017-02-01',"%Y-%m-%d"),"%Y%m%d%H%S")
        por_end_fmt = datetime.strftime(datetime.strptime('2017-02-28',"%Y-%m-%d"),"%Y%m%d%H%S")
        print(e.year-s.year)
        
        if unit[j] == 'degrees':
            
            args = {
                'start':por_start_fmt,
                'end':por_end_fmt,
                'obtimezone':'UTC',
                'vars': var[j],
                'stids': station_name,
                'token':token
            }
        else:
           args = {
               'start':por_start_fmt,
               'end':por_end_fmt,
               'obtimezone':'UTC',
               'vars': var[j],
               'stids': station_name,
               'units':unit[j],
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
        
        station_data_out[unit[j]] = pd.Series(value_pull,index=pd.to_datetime(dateTime))
    
    station_data_precip = station_data_out[unit[3]]
    station_data_wind = station_data_out[unit[2]] 
    station_data_temp = station_data_out[unit[1]]
    station_data_direction = station_data_out[unit[0]]
        
    station_data_precip.to_csv('raw_rainfalls_'+station_name+'.csv')
        
        


    # resample to only the hourly observations. For inter-hourly observations, the mean of these values is taken
    station_data_hourly = station_data_precip.resample('60min').mean()
    station_data_hourly_t = station_data_temp.resample('15min').mean()
    station_data_hourly_w = station_data_wind.resample('60min').mean()
    station_data_hourly_d = station_data_direction.resample('60min').first()
    station_data_hourly.to_csv('hourly_rainfalls_'+station_name+'.csv')
    station_data_hourly[station_data_hourly.isnull()] = 0
   # station_data_hourly_t = station_data_hourly_t[station_data_hourly_t.isnull()]
    #station_data_hourly_t[station_data_hourly_t.isnull()] = 0
    #station_data_hourly_w[station_data_hourly_w.isnull()] = 0
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
    #ts_selected[1]['start'] = pd.Timestamp('2017-2-6-22',tz='UTC')
    ts_selected_isolate = ts_selected[0:6]
    
    #split into pulses
    pulse_split = ts_selected.copy()
    def adjust_events(event):
        event['event_rainfall'] = event['event_rainfall'][event['start']:event['end'] ]
        event['event_perc'] = event['event_perc'][event['start']:event['end'] ]
        event['event_total'] = np.sum(event['event_rainfall'][event['start']:event['end'] ]) 
        event['event_avg'] = np.mean(event['event_rainfall'][event['start']:event['end'] ]) 
        
    ts_1 = pulse_split[0]
    ts_1['end'] = pd.Timestamp('2017-2-4-00',tz='UTC')
    ts_1['synoptic_event'] = 1
    ts_1['pulse_event'] = 1
    adjust_events(ts_1) 
    ts1_pd = ts_1['event_rainfall'].index[np.argsort(ts_1['event_rainfall'])][::-1][0:4].sort_values()
    
    
    ts_2 = pulse_split[1].copy()
    ts_2['start'] = pd.Timestamp('2017-2-6-00',tz='UTC')
    ts_2['end'] = pd.Timestamp('2017-2-6-08',tz='UTC')
    ts_2['synoptic_event'] = 2
    ts_2['pulse_event'] = 2
    adjust_events(ts_2) 
    ts2_pd = ts_2['event_rainfall'].index[np.argsort(ts_2['event_rainfall'])][::-1][0:4].sort_values()
    
    ts_3 = pulse_split[1].copy()
    ts_3['start']  = pd.Timestamp('2017-2-6-22',tz='UTC')
    ts_3['synoptic_event'] = 2
    ts_3['pulse_event'] = 3
    adjust_events(ts_3) 
    ts3_pd = ts_3['event_rainfall'].index[np.argsort(ts_3['event_rainfall'])][::-1][0:4].sort_values()
    
    ts_4 = pulse_split[2].copy()
    ts_4['start'] = pd.Timestamp('2017-2-8-07',tz='UTC')
    ts_4['end'] = pd.Timestamp('2017-2-8-22',tz='UTC')
    ts_4['synoptic_event'] = 3
    ts_4['pulse_event'] = 4
    adjust_events(ts_4) 
    ts4_pd = ts_4['event_rainfall'].index[np.argsort(ts_4['event_rainfall'])][::-1][0:4].sort_values()
    
    ts_5 = pulse_split[2].copy()
    ts_5['start'] = pd.Timestamp('2017-2-9-08',tz='UTC')
    ts_5['synoptic_event'] = 3
    ts_5['pulse_event'] = 5
    adjust_events(ts_5) 
    ts5_pd = ts_5['event_rainfall'].index[np.argsort(ts_5['event_rainfall'])][::-1][0:4].sort_values()
    
    ts_6 = pulse_split[3].copy()
    ts_6['synoptic_event'] = 4
    ts_6['pulse_event'] = 6
    adjust_events(ts_6) 
    ts6_pd = ts_6['event_rainfall'].index[np.argsort(ts_6['event_rainfall'])][::-1][0:4].sort_values()
    
    ts_7 = pulse_split[4].copy()
    ts_7['start'] = pd.Timestamp('2017-2-16-06',tz='UTC')
    ts_7['synoptic_event'] = 5
    ts_7['pulse_event'] = 7
    adjust_events(ts_7) 
    ts7_pd = ts_7['event_rainfall'].index[np.argsort(ts_7['event_rainfall'])][::-1][0:4].sort_values()
    
    ts_8 = pulse_split[5].copy()
    ts_8['start'] = pd.Timestamp('2017-2-17-05',tz='UTC')
    ts_8['end'] = pd.Timestamp('2017-2-18-18',tz='UTC')
    ts_8['synoptic_event'] = 6
    ts_8['pulse_event'] = 8
    adjust_events(ts_8) 
    ts8_pd = ts_8['event_rainfall'].index[np.argsort(ts_8['event_rainfall'])][::-1][0:4].sort_values()
    
    ts_9 = pulse_split[5].copy()
    ts_9['start'] = pd.Timestamp('2017-2-19-06',tz='UTC')
    ts_9['synoptic_event'] = 7
    ts_9['pulse_event'] = 9
    adjust_events(ts_9) 
    ts9_pd = ts_9['event_rainfall'].index[np.argsort(ts_9['event_rainfall'])][::-1][0:4].sort_values()
    
    
    paper_dates = [ts1_pd,ts2_pd,ts3_pd,ts4_pd,ts5_pd,ts6_pd,ts7_pd,ts8_pd,ts9_pd]
    

    ts_total = [ts_1]#[ts_1,ts_2,ts_3,ts_4,ts_5,ts_7,ts_8,ts_9]
    event_classification = []
    
    for i in ts_total:
        event_classification.append([i['synoptic_event'],i['pulse_event'],str(i['start'].month) +"-"+ str(i['start'].day) +"-"+ str(i['start'].hour),str(i['end'].month) +"-"+ str(i['end'].day) +"-"+ str(i['end'].hour),np.round(i['event_total'],1),np.round(i['event_avg'],1),np.round(np.max(i['event_rainfall']),1),i['end']-i['start']])
    ec = pd.DataFrame(event_classification)
    ec.columns = ['Synoptic Event #','Pulse Event #','Event start date (UTC)','Event end date (UTC)','Total precipiation (mm)','Average precipiation (mm)',"Maxium Precipitation (mm)",'Day Difference']
    ec.to_csv('Event_Classification.csv',index=False)
    
    ts_ams = [ts_1]
    #plot all pulse events
    fig = plt.figure(figsize=(20, 20))
    multimodal= 1
    pp_date = 0
    for i in ts_total:
        
        start_date = i['start']
        end_date = i['end']
        print(start_date)
        print(end_date)
        outdirck = 'G:\\NCFR Thesis\\NCFR_Thesis\\combined_'+station_name + '_'+str(start_date.year)+str(start_date.month)+str(start_date.day)+ str(end_date.hour)+'_'+ str(end_date.year)+ str(end_date.month)+ str(end_date.day)+ str(end_date.hour) +'\\'
        if os.path.exists(outdirck):
            shutil.rmtree(outdirck)
        os.mkdir(outdirck)
        
        start_date1 = i['start']
        end_date1 = i['end']
        station_data = station_data_hourly
        station_data_temp = station_data_hourly
        station_data_wind = station_data_hourly_t
        station_data_direction = station_data_hourly_d
        ts_selected = i
        station_name= station_name
        radar_name = 'KBBX'
        slat = station_lat
        slon = station_lon
        
        radar.pull_radar(i['start'], i['end'],station_data_hourly,station_data_hourly_t,station_data_hourly_w,station_data_hourly_d,i,station_name,'KBBX',station_lat,station_lon)
       
        cb.combine_figures(paper_dates[pp_date],outdirck)
        
        date_filter = station_data_hourly[i['start']:i['end']]
        date_filter_cumsum = date_filter.cumsum()
        ax = fig.add_subplot(4,2,multimodal)
        ax.bar(date_filter.index,date_filter,width=0.01)
        ax.set_title('Event #' + str(multimodal)  , fontsize=30)
        plt.xticks(fontsize=20,rotation=30)
        plt.yticks(fontsize=20)
        plt.ylim([0,4.5])
        ax.grid()
        date_form = DateFormatter("%d:%H")
        ax.xaxis.set_major_formatter(date_form)
        
        #fig.savefig(station_name + "_"+args['vars'] +str(multimodal)+ '.jpeg')  
        print(multimodal)
        multimodal = multimodal + 1
        pp_date = pp_date + 1
    fig.tight_layout(rect=[0.05, 0.05, 1, 0.93])   
    fig.supylabel('Precipiation (mm)', fontsize=35)
    fig.supxlabel('Day:Hour', fontsize=35)
    fig.suptitle("February 2017 Rainfall Pulses at Oroville Municipal Airport",fontsize=35)
    fig.savefig(station_name + "_Pulses"+ '.jpeg')
    plt.close('all')
    
    
   
    
    
    
    # fig1,ax1 = plt.subplots(figsize=(20, 12.5))
    
    cwe_series = {'CDEC_MFF_2002_2023_V2':'CDEC_MFF_FBS_2840','CDEC_NFF_2002_2023':'CDEC_NFF_BRS_3560','CDEC_UYB_2002_2023':'CDEC_UYB_PKC_3714'}
    cwe_titles = ['Forbestown (2840 ft)','Brush Creek (3560 ft)','Pike County (3714 ft)']
    colors = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']
    title_i = 1
    fig2,ax2 = plt.subplots(figsize=(40, 20))
    plot_text = ['(a)','(b)','(c)']
    for i in cwe_series.keys():
        pt = plot_text[title_i-1]
        additional_time_series = pd.read_csv('G:\\NCFR Thesis\\NCFR_Thesis\\'+i+'.csv')
        site = cwe_series[i].split('_')[2]
        #fig2,ax2 = plt.subplots(figsize=(20, 12.5))
        #ax2 = fig2.add_subplot(1,1,1)
        date_filter = additional_time_series.loc[additional_time_series['YYYYMMDDHH'].apply(str).str.contains('201702'),:]#station_data_hourly.loc[(station_data_hourly.index.year == 2017) & (station_data_hourly.index.month == 2)]
        dt = pd.to_datetime(date_filter['YYYYMMDDHH'].apply(str),format='%Y%m%d%H')
        ax2.bar(dt[date_filter[cwe_series[i]]>=0],date_filter[cwe_series[i]][date_filter[cwe_series[i]]>=0],width=0.15,color=colors[title_i-1],label=cwe_titles[title_i-1])
        #fig2.suptitle(cwe_titles[title_i-1], fontsize=30)
        #ax2.set_title(str(cwe_titles[title_i-1])  , fontsize=40)
        #date_form = DateFormatter("%d-%Y-%H")
        #fig2.savefig(site.upper()+'.jpeg')
        title_i = title_i + 1
        
        #plt.close('all')
    fig2,ax2 = plt.subplots(figsize=(40, 20))
    date_filter = station_data_hourly.loc[(station_data_hourly.index.year == 2017) & (station_data_hourly.index.month == 2)]
    ax2.bar(date_filter.index,date_filter,color=colors[3],width=0.15,label="Oroville Municipal Airport (194 ft)")
    fig2.suptitle("Hourly Precipitation Record near \n Oroville Dam in February 2017", fontsize=50,fontweight='bold')
    plt.ylabel('Precipiation (mm)', fontsize=70)
    plt.xlabel('Day of Month', fontsize=70)
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=60)
    ax2.grid()
    plt.legend(fontsize=39,loc='upper right')
    #ax2.set_title('Oroville Municipal Airport' , fontsize=40)
    date_form = DateFormatter("%d")
    ax2.xaxis.set_major_formatter(date_form)
    ax2 = fig2.gca()
    ax2.set_ylim([0, None])
    #fig2.tight_layout(rect=[0.05, 0.05, 1, 0.93]) 
    fig2.savefig('Precip_Time_Series_Oroville_alone'+'.jpeg')
    plt.close('all')
    
    fig3,ax3 = plt.subplots(figsize=(40, 20))
    #date_filter = station_data_hourly.loc[(station_data_hourly.index.year == 2017) & (station_data_hourly.index.month == 2)]
    ax3.bar(ts_9['event_rainfall'].index[21:],ts_9['event_rainfall'][21:],color=colors[0],width=0.01,label="Pulse Example")
    fig3.suptitle("Precipitation Pulse Inspection Example", fontsize=50,fontweight='bold')
    plt.ylabel('Precipiation (mm)', fontsize=70)
    plt.xlabel('Day:Hour', fontsize=70)
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=60)
    ax3.grid()
    plt.legend(fontsize=39,loc='upper right')
    #ax2.set_title('Oroville Municipal Airport' , fontsize=40)
    date_form = DateFormatter("%d:%H")
    ax3.xaxis.set_major_formatter(date_form)
    ax3 = fig3.gca()
    ax3.set_ylim([0, None])
    #fig2.tight_layout(rect=[0.05, 0.05, 1, 0.93]) 
    fig3.savefig('Precip_Pulse_Example'+'.jpeg')

        
        