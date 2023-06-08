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

import MRMS_data_pull as MRMS

# stage four quanitiative precipitation (hourly radar product)
# find the ARs that are colocated
# Specify request parameters (as strings)
token = 'aa31874e86fb42d9b2ea6b293f1bb004' # use your own token

def count_and_sum_events(time_series):
    event_start = []
    event_end = []
    event_acum = []
    for i in np.arange(0,len(time_series)-1):
        if time_series[i] != 0 and np.all(time_series[i-24:i] == 0):
            event_start.append(time_series.index[i])
            start_i = i
        if time_series[i] != 0 and np.all(time_series[i+1:i+25] == 0):
            event_end.append(time_series.index[i])
            end_i = i
            event_total = np.sum(time_series[start_i:end_i+1])
            event_acum.append(np.sum(time_series[start_i:end_i+1]))
            
    return(np.array([event_start,event_end,event_acum]))

def find_nearest_value(target, values):
    nearest = None
    min_difference = float('inf')

    for value in values:
        difference = abs(value - target)
        if difference < min_difference:
            min_difference = difference
            nearest = value

    return nearest

#open AR time series
with open('D:\\PSU Thesis\\data\\AR_California_Landfall.pickle', 'rb') as data:
    AR_TS = pickle.load(data)
    
def myround(x, prec=2, base=.5):
  return round(base * round(float(x)/base),prec)

#def mergedata(AR_TS,event):
    #merge AR with event by date and AR landfall lat

def data_availability_check(station_name,token,var):
    if responseDict_por['SUMMARY']['RESPONSE_MESSAGE'] == 'OK':
        return(True)
    else:
        return(False)

AR_point = []

for i in range(0,len(AR_TS)):
    latlon_time = [AR_TS[i][0],AR_TS[i][1]-360,AR_TS[i][2].astype(object).year,AR_TS[i][2].astype(object).month,AR_TS[i][2].astype(object).day,AR_TS[i][2].astype(object).hour]
    AR_point.append(latlon_time)
AR_pd = pd.DataFrame(AR_point,columns=["latitude","longitude","year","month","day","hour"])

#define boundaries for reserviors
AR_pd["Yuba-Feather"] = (AR_pd.latitude >= 37) & (AR_pd.latitude <= 40)
AR_pd["Santa Ana"] = (AR_pd.latitude >= 32) & (AR_pd.latitude <= 37)
AR_pd["Russian River"] = (AR_pd.latitude >= 38) & (AR_pd.latitude <= 40)

    
    
#['ksfo','klax','krdd','ksmo','kcic','kpdx']'ksfo','klax','krdd','ksmo','kcic',#

#'kuki','ksts','kapc','k069','kmvy','kove','kgoo','ko05','ktrk','kblu','kcno','kral','kl35','kajo'


#russian river
#kuki = Ukiah
#ksts = santa rosa
#k069 = petaluma
#kapc = Napa County

#yuba feather
#kmvy = Yuba city
#kove = oroville
#kgoo = nevada county
#ko05 = rogers field airport (lake alamour): no precip data
#ktrk = lake tahoe
#kblu = emigrant pass (Blue Canyon)

#santa ana
#kcno = Chino
#klgb = Long Beach
#kral = Riverside
#kl35 = Big Bear
#kajo = Corona

station_name_list = ['kcic','kuki','ksts','kapc','k069','kmvy','kove','kgoo','ko05','ktrk','kblu','kcno','kl35','kajo','kral']#["kmvy","kove"]#,"kcno",'klax','krdd','ksmo','kcic','kpdx','ksfo','klax','krdd','ksmo','kcic','C3BCC','C3BVS','C3CAT','C3DLA','C3DRW','C3FRC','C3GPO','C3HDC','C3HRD','C3NBB','C3NCM','C3POR','C3PVN','C3SKI','C3SKY','C3SOD','C3WDG','C3WPO']
#station_boundary = ['yfrr','yfrr','yfrr','yfrr']

station_window = {'kcic':'yfrr','kuki':'yfrr','ksts':'yfrr','kapc':'yfrr','k069':'yfrr','kmvy':'yfrr','kove':'yfrr','kgoo':'yfrr','ko05':'yfrr','ktrk':'yfrr','kblu':'yfrr','kcno':'sa','kl35':'sa','kajo':'sa','kral':'sa'}

AR_Catalog = pd.read_excel('D:\\PSU Thesis\\data\\ARcatalog_NCEP_NEW_1948-2018_Comprehensive_FINAL_29JAN18.xlsx',"AR_Events")
# check for highest resolution of precip data 


var = 'precip_accum_one_hour'
unit = 'precip|mm'
plot_full_record = False
output_5_perc_data = 'D:\\PSU Thesis\\data\\precip_threshold_dates\\'

#define cumsum dates
start_date = dtetme.datetime(2019,1,12,9)
end_date = dtetme.datetime(2019,1,13,1)



for i in station_name_list:
    print(i)
    bound = station_window[i]
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
        
        #find AR lat
        near_AR_lat = find_nearest_value(station_lat,AR_pd['latitude'])
        near_AR_lat_p1 = near_AR_lat+0.5
        near_AR_lat_m1 = near_AR_lat-0.5
        
        AR_intersect = AR_pd.loc[(AR_pd['latitude']==near_AR_lat) | (AR_pd['latitude']==near_AR_lat_p1) | (AR_pd['latitude']==near_AR_lat_m1),:]
        AR_intersect['date'] = pd.to_datetime(AR_intersect[['year','month','day']])
        
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
        
        s_e_sum = count_and_sum_events(station_data_hourly)
        
        sum_events = s_e_sum[2]
        
        sum_events_uq_thres = np.quantile(sum_events,0.75)
        
        thres_boolean = sum_events >= sum_events_uq_thres
        
        uq_start = s_e_sum[0][thres_boolean]
        uq_end = s_e_sum[1][thres_boolean]
        sum_events_uq = sum_events[thres_boolean]
        
        time_spans = []
        for i in np.arange(0,len(sum_events_uq)-1):
            time_spans.append({'start':uq_start[i],'end':uq_end[i],'event_total':sum_events_uq[i]})
        
        event_id = []
        for r in time_spans:
            event_percentage = 100 * station_data_hourly[r['start']:r['end']]/r['event_total']
            dates_of_event = pd.DataFrame(np.unique(pd.to_datetime(station_data_hourly[r['start']:r['end']].index.date)))
            AR_DF = np.unique(pd.to_datetime(AR_intersect['date']))
            event_id.append({'start':r['start'],'end':r['end'],'event_rainfall':station_data_hourly[r['start']:r['end']],'event_perc':event_percentage,'associated_AR':any(any(row) for row in dates_of_event.isin(AR_DF).values)})
            
        #limit to MRMS data    
       # events_with_AR = [item for item in event_id if item['associated_AR'] ==True and item['start'].year >= 2016]
        events_with_AR = [item for item in event_id if item['associated_AR'] ==True and (item['start'].year >= 2016) and (any(item['event_perc'] >= 20)) ]
            #station_data_hourly[r['start']:r['end']]
            
            
        
        

            
        

        
        #non-zero rainfall reports
        station_data_hourly_nz = station_data_hourly.loc[station_data_hourly > 0]
        station_data_hourly_nnull = station_data_hourly.loc[station_data_hourly.notnull()]
        
        #compute statistics
        percentile_75 = station_data_hourly.quantile(0.75)
        percentile_95 = station_data_hourly.quantile(0.95)
        percentile_99 = station_data_hourly.quantile(0.99)
        
        percentile_75_nz = station_data_hourly_nz.quantile(0.75)
        percentile_95_nz = station_data_hourly_nz.quantile(0.95)
        percentile_99_nz = station_data_hourly_nz.quantile(0.99)
        
        
        threshold_exceedance_points_95_nz = station_data_hourly.loc[station_data_hourly >= percentile_95_nz]
        
        threshold_exceedance_points_75_nz = station_data_hourly.loc[station_data_hourly >= percentile_75_nz]
        
        thres_to_plot = threshold_exceedance_points_95_nz.groupby([threshold_exceedance_points_95_nz.index.year, 
                                                                 threshold_exceedance_points_95_nz.index.month,
                                                                 threshold_exceedance_points_95_nz.index.day]).count()
        thres_to_plot_75 = threshold_exceedance_points_75_nz.groupby([threshold_exceedance_points_75_nz.index.year, 
                                                                 threshold_exceedance_points_75_nz.index.month,
                                                                 threshold_exceedance_points_75_nz.index.day]).count()
        thres_to_plot_75[thres_to_plot_75>1]
        
        
        thres_to_plot_75[thres_to_plot_75>1].to_csv('precip_nonzero_hourly.csv')
        threshold_exceedance_points_75_nz.to_csv(station_name + ' 75th_percentile.csv')
        #plotting for slide
        
        fig,ax = plt.subplots(figsize=(40, 20))
        thres_to_plot.plot(kind="bar")
        ax.set_xticks(np.arange(1,len(thres_to_plot),1))
        ax.set_xticks(np.arange(1,len(thres_to_plot),5))
        fig.suptitle('# of Rainfall Events per day at '+station_name.upper() +' with Hourly Precipitation Rate above 95th Percentile', fontsize=60)
        plt.ylabel('Precipiation (mm)', fontsize=50)
        plt.xlabel('Date', fontsize=50)
        plt.xticks(fontsize=40,rotation=45)
        plt.yticks(fontsize=40)
        date_form = DateFormatter("%m-%d-%Y")
        ax.grid()
        ax.xaxis.set_major_formatter(date_form)
        fig.savefig(station_name + "_"+args['vars'] + '_hist.png')
        
        fig1,ax1 = plt.subplots(figsize=(40, 25))
        ax1.plot(station_data_hourly.index.date,station_data_hourly)
        #ax1.set_xticks(np.arange(1,len(thres_to_plot),1))
        #ax1.set_xticks(np.arange(1,len(thres_to_plot),5))
        fig1.suptitle('Hourly Precipiation at Riverside, CA ('+station_name.upper()+ ')', fontsize=60)
        plt.ylabel('Precipiation (mm)', fontsize=50)
        plt.xlabel('Date', fontsize=50)
        plt.xticks(fontsize=30,rotation=45)
        plt.yticks(fontsize=30)
        ax1.grid()
        date_form = DateFormatter("%m-%d-%Y")
        ax1.xaxis.set_major_formatter(date_form)
        fig1.savefig(station_name + "_"+args['vars'] + '_line.png')
        
        
        
        threshold_exceedance_points_95_nz.to_csv(output_5_perc_data+'5_perc_rainfalls_'+station_name+'.csv')
        
        threshold_exceedance_points_95 = station_data_hourly.loc[station_data_hourly >= percentile_95]
        threshold_exceedance_points_99 = station_data_hourly.loc[station_data_hourly >= percentile_99]
        
        
        threshold_grouping_95 = threshold_exceedance_points_95.groupby([threshold_exceedance_points_95.index.year, 
                                                                 threshold_exceedance_points_95.index.month,
                                                                 threshold_exceedance_points_95.index.day]).count()
        
        threshold_grouping_99 = threshold_exceedance_points_99.groupby([threshold_exceedance_points_99.index.year, 
                                                                 threshold_exceedance_points_99.index.month,
                                                                 threshold_exceedance_points_99.index.day]).count()
        
        cum_sum_precip = station_data_hourly.groupby([station_data_hourly.index.year,
                                                 station_data_hourly.index.month,
                                                 station_data_hourly.index.day]).cumsum()
        
        
        cum_sum_precip_t = cum_sum_precip.reset_index()
        
        cols= ["index"]
      #  cum_sum_precip_t['Date'] = datetime.strptime(cum_sum_precip_t["Date"],"%Y-%m-%dT%H:%M:%S.%fZ")
        cum_sum_precip_t['Date'] = cum_sum_precip_t[cols].apply(lambda x: '-'.join(x.values.astype(str)), axis="columns")
        cum_sum_precip_t['Date'] = pd.to_datetime(cum_sum_precip_t["Date"])
        
        cum_precip = station_data_hourly.groupby([station_data_hourly.index.year,
                                                 station_data_hourly.index.month,
                                                 station_data_hourly.index.day]).sum()
        
        
        t=cum_sum_precip_t

        t=t.reset_index()


        ar_cols = ['year','month','day']
        
        AR_pd['Date'] = AR_pd[ar_cols].apply(lambda x: '-'.join(x.values.astype(str)), axis="columns")
       # AR_pd['Date'] = pd.to_datetime(AR_pd["Date"])
        t['mergedate'] = t.Date.dt.year.astype(str) + "-" + t.Date.dt.month.astype(str) + "-" + t.Date.dt.day.astype(str)

        t2 = t.merge(AR_pd, left_on='mergedate', right_on='Date',how="left")
        YRR = t2.loc[(t2["Yuba-Feather"]==True)]
        SA = t2.loc[(t2["Santa Ana"]==True)]
        
        merged_cumsum = t2[["index", 0, "Yuba-Feather","Santa Ana"]].drop_duplicates()
        merged_cumsum.columns = ["Date","precip_mm","Yuba_Russian","Santa_Ana"]
      
        
     #isolate day, plot cumulation
        if(plot_full_record):
            date_filter_cumsum = merged_cumsum
        else: 
            date_filter_cumsum = station_data_hourly[start_date:end_date].cumsum()
            
        
        fig3,ax3 = plt.subplots(figsize=(100, 60))
        ax3.grid()
        
        
        cumsum_date = date_filter_cumsum.index
        cumsum_precip_value  = date_filter_cumsum
  
        
        #t_Yuba_Feather_RR.plot(ax=ax2,x="Date",y=0,kind="bar",color="red",label="Yuba Feather or Russian River AR")
        #t_Santa_Ana.plot(ax=ax2,x="Date",y=0,kind="bar",color="blue",label="Santa Ana AR")
        ax3.plot(cumsum_date, date_filter_cumsum,color="red")
        #ax3.plot(Yuba_Russian.Date,Yuba_Russian.precip_mm,color="green",label="Yuba Feather/Russian River")
        #ax3.plot(Santa_Ana.Date,Santa_Ana.precip_mm,color="blue",label="Santa Ana")
        #t2.plot(ax=ax2,x="Date",y=0,kind="bar",color="green",label="all")
        ax3.set_ylabel('daily cumlative precip (mm)',fontsize=70)
        ax3.set_title('Precip Daily Cumulative Accumulation for ' + station_name,fontsize=70)
        ax3.tick_params(axis='both', which='major', labelsize=40)
        ax3.legend(loc="upper left",fontsize=70)
        #ax3.set_xticks(np.arange(1,len(cumsum_date),1))
        #ax3.set_xticks(np.arange(1,len(cumsum_date),10))
        #fig3.savefig(station_name + "_"+args['vars'] + '_cumsum.png')
        #  #cum_sum_precip_t.to_csv('test.csv')
        # # station_data_out.to_csv('test2.csv')
        print(station_name)
        for i in events_with_AR:
            MRMS.map_event(i['start'].to_pydatetime(), i['end'].to_pydatetime(), station_lon, station_lat, bound,station_name)
            
    else:
        print("No data available for " + station_name)
        
        