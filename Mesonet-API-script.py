import urllib
import json
import pandas as pd
import matplotlib.pyplot as plt
import datetime
from datetime import datetime

# Specify request parameters (as strings)
token = 'aa31874e86fb42d9b2ea6b293f1bb004' # use your own token

#test
# Create API query string

def data_availability_check(station_name,token,var):
    if responseDict_por['SUMMARY']['RESPONSE_MESSAGE'] == 'OK':
        return(True)
    else:
        return(False)

station_name_list = 'C3WDG'#['ksfo','klax','krdd']#['C3BCC','C3BVS','C3CAT','C3DLA','C3DRW','C3FRC','C3GPO','C3HDC','C3HRD','C3NBB','C3NCM','C3POR','C3PVN','C3SKI','C3SKY','C3SOD','C3WDG','C3WPO'],

# check for highest resolution of precip data 

var = 'precip_accum_one_hour'


for i in station_name_list:
    station_name = i
    args_por = {
        'obtimezone':'UTC',
        'vars':var,
        'stids': station_name,
        'units':'precip|mm',
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
            'units':'precip|mm',
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
        #station_data_hourly = station_data_out.where(station_data_out.index.minute == 53).dropna()
        
        # Plotting code
        fig,ax = plt.subplots(figsize=(20, 20))
        plt.plot(station_data_out,linewidth=3)
        ax.set_ylabel(args['units'],fontsize=40)
        ax.tick_params(axis='both', which='major', labelsize=40)
        ax.set_title(args['stids'] + ' Hourly Precip (mm)',fontsize=40)
        # Clean up the plot a bit
        from matplotlib.dates import DayLocator, DateFormatter
        #ax.xaxis.set_major_locator(DayLocator())
        #ax.xaxis.set_major_formatter(DateFormatter('%d-%b'))
        ax.grid()
        
        fig.savefig('temp2_'+station_name+'.png')
        
        #compute 95 percentile threshold
        
        percentile_95 = station_data_out.quantile(0.95)
        
        threshold_exceedance_points = station_data_out.loc[station_data_out >= percentile_95]
        
        
        threshold_grouping= threshold_exceedance_points.groupby([threshold_exceedance_points.index.year, 
                                                                 threshold_exceedance_points.index.month,
                                                                 threshold_exceedance_points.index.day]).count()
       
        fig1,ax1 = plt.subplots(figsize=(300, 50))
        threshold_grouping.plot(kind="bar")
        ax1.set_ylabel('counts',fontsize=70)
        ax1.set_title('Precip occurrances for ' + station_name + ' where hourly rainfall rates exceed ' + str(percentile_95) + ' mm',fontsize=70)
        ax1.tick_params(axis='both', which='major', labelsize=40)
        fig1.savefig('freq_over_exceedance_'+station_name+'.png')
else:
    print("No data available")