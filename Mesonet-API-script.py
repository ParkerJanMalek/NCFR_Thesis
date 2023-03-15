import urllib
import json
import pandas as pd
import matplotlib.pyplot as plt
import datetime
from datetime import datetime
import scipy.io
import numpy as np
import pickle
from matplotlib.dates import DayLocator, DateFormatter

# stage four quanitiative precipitation (hourly radar product)
# find the ARs that are colocated
# Specify request parameters (as strings)
token = 'aa31874e86fb42d9b2ea6b293f1bb004' # use your own token


#open AR time series
with open('D:\\PSU Thesis\\data\\AR_California_Landfall.pickle', 'rb') as data:
    AR_TS = pickle.load(data)

AR_point = []

for i in range(0,len(AR_TS)):
   latlon_time = [AR_TS[i][0],AR_TS[i][1]-360,AR_TS[i][2].astype(object).year,AR_TS[i][2].astype(object).month,AR_TS[i][2].astype(object).day,AR_TS[i][2].astype(object).hour]
   AR_point.append(latlon_time)
AR_pd = pd.DataFrame(AR_point,columns=["latitude","longitude","year","month","day","hour"])

#define boundaries for reserviors
AR_pd["Yuba-Feather"] = (AR_pd.latitude >= 37) & (AR_pd.latitude <= 40)
AR_pd["Santa Ana"] = (AR_pd.latitude >= 32) & (AR_pd.latitude <= 37)

def data_availability_check(station_name,token,var):
    if responseDict_por['SUMMARY']['RESPONSE_MESSAGE'] == 'OK':
        return(True)
    else:
        return(False)
    
    
#['ksfo','klax','krdd','ksmo','kcic','kpdx']'ksfo','klax','krdd','ksmo','kcic',#

#russian river
#kuki = Ukiah

#yuba feather
#kmvy = Yuba city
#kove = oroville

#santa ana
#kcno = chino
station_name_list = ["klax"]#["kmvy","kove"]#,"kcno",'klax','krdd','ksmo','kcic','kpdx','ksfo','klax','krdd','ksmo','kcic','C3BCC','C3BVS','C3CAT','C3DLA','C3DRW','C3FRC','C3GPO','C3HDC','C3HRD','C3NBB','C3NCM','C3POR','C3PVN','C3SKI','C3SKY','C3SOD','C3WDG','C3WPO']


AR_Catalog = pd.read_excel('D:\\PSU Thesis\\data\\ARcatalog_NCEP_NEW_1948-2018_Comprehensive_FINAL_29JAN18.xlsx',"AR_Events")
# check for highest resolution of precip data 


var = 'precip_accum_one_hour'
unit = 'precip|mm'
plot_full_record = True
output_5_perc_data = 'D:\\PSU Thesis\\data\\precip_threshold_dates\\'

for i in station_name_list:

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
        
        por_start_fmt = datetime.strftime(datetime.strptime(por_start[0:10],"%Y-%m-%d"),"%Y%m%d%H%S")
        por_end_fmt = datetime.strftime(datetime.strptime(por_end[0:10],"%Y-%m-%d"),"%Y%m%d%H%S")
        
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
        
        #non-zero rainfall reports
        station_data_hourly_nz = station_data_hourly.loc[station_data_hourly > 0]
        station_data_hourly_nnull = station_data_hourly.loc[station_data_hourly.notnull()]
        
        #compute statistics
        
        percentile_95 = station_data_hourly.quantile(0.95)
        percentile_99 = station_data_hourly.quantile(0.99)
        
        percentile_95_nz = station_data_hourly_nz.quantile(0.95)
        percentile_99_nz = station_data_hourly_nz.quantile(0.99)
        
        
        threshold_exceedance_points_95_nz = station_data_hourly.loc[station_data_hourly >= percentile_95_nz]
        
        thres_to_plot = threshold_exceedance_points_95_nz.groupby([threshold_exceedance_points_95_nz.index.year, 
                                                                 threshold_exceedance_points_95_nz.index.month,
                                                                 threshold_exceedance_points_95_nz.index.day]).count()
        
        
        
        #plotting for slide
        
        fig,ax = plt.subplots(figsize=(40, 20))
        thres_to_plot.plot(kind="bar")
        ax.set_xticks(np.arange(1,len(thres_to_plot),1))
        ax.set_xticks(np.arange(1,len(thres_to_plot),5))
        fig.suptitle('# of Rainfall Events per day at KLAX with Hourly Precipitation Rate above 95th Percentile', fontsize=60)
        plt.ylabel('Precipiation (mm)', fontsize=50)
        plt.xlabel('Date', fontsize=50)
        plt.xticks(fontsize=20,rotation=45)
        plt.yticks(fontsize=20)
        date_form = DateFormatter("%m-%d-%Y")
        ax.xaxis.set_major_formatter(date_form)
        fig.savefig(station_name + "_"+args['vars'] + '_hist.png')
        
        fig1,ax1 = plt.subplots(figsize=(40, 25))
        ax1.plot(station_data_hourly_nz.index.date,station_data_hourly_nz)
        #ax1.set_xticks(np.arange(1,len(thres_to_plot),1))
        #ax1.set_xticks(np.arange(1,len(thres_to_plot),5))
        fig1.suptitle('Hourly Precipiation at KLAX', fontsize=60)
        plt.ylabel('Precipiation (mm)', fontsize=50)
        plt.xlabel('Date', fontsize=50)
        plt.xticks(fontsize=30,rotation=45)
        plt.yticks(fontsize=30)
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
            #define min and max date for range
            minDate = pd.to_datetime("2001-01-10")
            maxDate = pd.to_datetime("2001-01-11")
            date_filter_cumsum = merged_cumsum.loc[((merged_cumsum["Date"].dt.month >= minDate.month) & (merged_cumsum["Date"].dt.month <= maxDate.month))
                                                       & ((merged_cumsum["Date"].dt.day >= minDate.day) & (merged_cumsum["Date"].dt.day <= maxDate.day))
                                                       & ((merged_cumsum["Date"].dt.year >= minDate.year) & (merged_cumsum["Date"].dt.year <= maxDate.year)) ]
            
        
        fig3,ax3 = plt.subplots(figsize=(100, 60))
        ax3.grid()
        
        
        cumsum_date = date_filter_cumsum["Date"]
        cumsum_precip_value  = date_filter_cumsum["precip_mm"]
        
        Yuba_Russian = date_filter_cumsum.loc[date_filter_cumsum.Yuba_Russian == True]
        Santa_Ana = date_filter_cumsum.loc[date_filter_cumsum.Santa_Ana == True]
        
        #t_Yuba_Feather_RR.plot(ax=ax2,x="Date",y=0,kind="bar",color="red",label="Yuba Feather or Russian River AR")
        #t_Santa_Ana.plot(ax=ax2,x="Date",y=0,kind="bar",color="blue",label="Santa Ana AR")
        ax3.plot(cumsum_date, cumsum_precip_value,color="red",label="unaffiliated")
        #ax3.plot(Yuba_Russian.Date,Yuba_Russian.precip_mm,color="green",label="Yuba Feather/Russian River")
        #ax3.plot(Santa_Ana.Date,Santa_Ana.precip_mm,color="blue",label="Santa Ana")
        #t2.plot(ax=ax2,x="Date",y=0,kind="bar",color="green",label="all")
        ax3.set_ylabel('daily cumlative precip (mm)',fontsize=70)
        ax3.set_title('Precip Daily Cumulative Accumulation for ' + station_name,fontsize=70)
        ax3.tick_params(axis='both', which='major', labelsize=40)
        ax3.legend(loc="upper left",fontsize=70)
        #ax3.set_xticks(np.arange(1,len(cumsum_date),1))
        #ax3.set_xticks(np.arange(1,len(cumsum_date),10))
        fig3.savefig(station_name + "_"+args['vars'] + '_cumsum.png')
        #  #cum_sum_precip_t.to_csv('test.csv')
        # # station_data_out.to_csv('test2.csv')
        print(station_name)
    else:
        print("No data available for " + station_name)
        
        