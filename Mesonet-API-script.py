import urllib
import json
import pandas as pd
import matplotlib.pyplot as plt
import datetime
from datetime import datetime
import scipy.io
import numpy as np


# f

# checks for missingness in kstatus for each timestep
# 
# for i in range(0,len(kstatus[0,:])):
#     print(str(np.all(np.isnan(np.array(kstatus[:,i])))) + " " + str(i) + " " + str(float(kstatus.lat[i])))

# for i in range(0,len(landfall[0,:,0])):
#     for j in range(0,len(landfall[0,0,:])):
#         print(str(np.all(np.isnan(np.array(landfall[:,i,j])))) + " " 
#               + str(i) + "," + str(j) + " " + str(float(landfall.lat[i]))
#               + "," +str(float(landfall.lon[j])))


#AR_termination_CA = kstatus[:,(kstatus.lat<= 60) & (kstatus.lat >= 0)]
#latlon[:,(latlon.lat<= 40) & (latlon.lat >= 37.5),(latlon.lon<= 180+125) & (latlon.lon >= 180+110)]
AR_termination_CA = kstatus.isel((kstatus.lat<= 60) & (kstatus.lat >= 0))



AR_Catalog_MERRA = pd.DataFrame(fds.time,fds.lat,fds.lon)
# Specify request parameters (as strings)
token = 'aa31874e86fb42d9b2ea6b293f1bb004' # use your own token

#test
# Create API query string


def data_availability_check(station_name,token,var):
    if responseDict_por['SUMMARY']['RESPONSE_MESSAGE'] == 'OK':
        return(True)
    else:
        return(False)
#['ksfo','klax','krdd','ksmo','kcic','kpdx']'ksfo','klax','krdd','ksmo','kcic',#
station_name_list = ['klax']#,'krdd','ksmo','kcic','kpdx','ksfo','klax','krdd','ksmo','kcic']#,'C3BCC','C3BVS','C3CAT','C3DLA','C3DRW','C3FRC','C3GPO','C3HDC','C3HRD','C3NBB','C3NCM','C3POR','C3PVN','C3SKI','C3SKY','C3SOD','C3WDG','C3WPO']


AR_Catalog = pd.read_excel('D:\\PSU Thesis\\data\\ARcatalog_NCEP_NEW_1948-2018_Comprehensive_FINAL_29JAN18.xlsx',"AR_Events")
# check for highest resolution of precip data 


ar_DMY = pd.DataFrame([AR_Catalog["Day"],AR_Catalog["Month"],AR_Catalog["Year"],
                        AR_Catalog["Russian River"],AR_Catalog["Yuba-Feather"],
                        AR_Catalog["Santa Ana"],AR_Catalog["California_AR"]]).T


cols=["Year","Month","Day"]

ar_DMY['Date'] = ar_DMY[cols].apply(lambda x: '-'.join(x.values.astype(str)), axis="columns")

ar_CA = ar_DMY.loc[ar_DMY["California_AR"] == "yes"]

ar_Dedupe = ar_CA.drop_duplicates(subset=['Year','Month','Day','California_AR'])


var = 'precip_accum_one_hour'
unit = 'precip|mm'



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
        
        # Retain only the hourly observations
        station_data_hourly = station_data_out.resample('60min').mean()
        #station_data_hourly = station_data_out.groupby(pd.Grouper(freq='60Min', offset=0, label='right')).first()
        
        # Plotting code
        # fig,ax = plt.subplots(figsize=(20, 20))
        # plt.plot(station_data_hourly,linewidth=3)
        # ax.set_ylabel(args['units'],fontsize=40)
        # ax.tick_params(axis='both', which='major', labelsize=40)
        # plt.xticks(rotation=90)
        # ax.set_title(args['stids'] + ' '+ args['units'],fontsize=40)
        # Clean up the plot a bit
        from matplotlib.dates import DayLocator, DateFormatter
        #ax.xaxis.set_major_locator(DayLocator())
        #ax.xaxis.set_major_formatter(DateFormatter('%d-%b'))
        # ax.grid()
        
        # fig.savefig(args['vars'] + '_'+station_name+'.png')
        
        
        #compute 95 percentile threshold
        
        percentile_95 = station_data_hourly.quantile(0.95)
        percentile_99 = station_data_hourly.quantile(0.99)
        
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
        
        
        t=cum_precip

        t=t.reset_index()


        cols=["level_0","level_1","level_2"]
        
        t['Date'] = t[cols].apply(lambda x: '-'.join(x.values.astype(str)), axis="columns")


        t2 = t.merge(ar_Dedupe, left_on='Date', right_on='Date',how="left")
        
        #t2['Date']=pd.datetime(t2['Date'])
        date_begin = t["Date"][0]
        date_end = t["Date"][len(t["Date"])-1]
        t_Yuba_Feather_RR = t2["Yuba-Feather"] == "yes"
        t_Santa_Ana = t2["Santa Ana"] == "yes"
       
        # fig1,ax1 = plt.subplots(figsize=(400, 50))
        # threshold_grouping_95.plot(kind="bar")
        # ax1.set_ylabel('counts',fontsize=70)
        # ax1.set_title('Precip occurrances for ' + station_name + ' where hourly rainfall rates exceed ' + str(percentile_95) + ' mm',fontsize=70)
        # ax1.tick_params(axis='both', which='major', labelsize=40)
        # fig1.savefig(args['vars'] + '_exceedance_'+station_name+'.png')
        
        fig2,ax2 = plt.subplots(figsize=(400, 50))
        #t_Yuba_Feather_RR.plot(ax=ax2,x="Date",y=0,kind="bar",color="red",label="Yuba Feather or Russian River AR")
       # t_Santa_Ana.plot(ax=ax2,x="Date",y=0,kind="bar",color="blue",label="Santa Ana AR")
        ax2.bar(t2["Date"],t2[0],color="red",label="Unassociated")
        ax2.bar(t2["Date"][t_Yuba_Feather_RR],t2[0][t_Yuba_Feather_RR],color="green",label="Yuba Feather/Russian River")
        ax2.bar(t2["Date"][t_Santa_Ana],t2[0][t_Santa_Ana],color="blue",label="Santa Ana")
        #t2.plot(ax=ax2,x="Date",y=0,kind="bar",color="green",label="all")
        ax2.set_ylabel('daily precip accumulation (mm)',fontsize=70)
        ax2.set_title('Precip test',fontsize=70)
        ax2.tick_params(axis='both', which='major', labelsize=40)
        ax2.legend(loc="upper left",fontsize=70)
        ax2.set_xticks(np.arange(1,len(t["Date"]),200))
        fig2.savefig(args['vars'] + 'test'+station_name+'TEST2.png')
        
        t2[t_Yuba_Feather_RR].to_csv(station_name + '_YubaFeather.csv')
        t2[t_Santa_Ana].to_csv(station_name + '_SantaAna.csv')
        t2.to_csv(station_name + '.csv')
        
     #isolate day, plot cumulation
     
        Date = pd.to_datetime("2010-12-19")
        fig3,ax3 = plt.subplots(figsize=(60, 60))
        cumsum_date = cum_sum_precip_t["Date"].loc[(cum_sum_precip_t["Date"].dt.month == Date.month) & (cum_sum_precip_t["Date"].dt.day == Date.day) & (cum_sum_precip_t["Date"].dt.year == Date.year) ]
        cumsum_precip_value  = cum_sum_precip_t[0].loc[(cum_sum_precip_t["Date"].dt.month == Date.month) & (cum_sum_precip_t["Date"].dt.day == Date.day) & (cum_sum_precip_t["Date"].dt.year == Date.year) ]
        #t_Yuba_Feather_RR.plot(ax=ax2,x="Date",y=0,kind="bar",color="red",label="Yuba Feather or Russian River AR")
       # t_Santa_Ana.plot(ax=ax2,x="Date",y=0,kind="bar",color="blue",label="Santa Ana AR")
        ax3.plot(cumsum_date,cumsum_precip_value,color="red")
        #t2.plot(ax=ax2,x="Date",y=0,kind="bar",color="green",label="all")
        ax3.set_ylabel('daily precip accumulation (mm)',fontsize=70)
        ax3.set_title('Precip test',fontsize=70)
        ax3.tick_params(axis='both', which='major', labelsize=40)
       # ax2.legend(loc="upper left",fontsize=70)
        #ax3.set_xticks(np.arange(1,len(cumsum_date),1))
        fig3.savefig(station_name + "_"+args['vars'] + 'cumsumt2010.png')
        #cum_sum_precip_t.to_csv('test.csv')
       # station_data_out.to_csv('test2.csv')
        print(station_name)
    else:
        print("No data available for " + station_name)
        
        