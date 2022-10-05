import urllib
import json
import pandas as pd
import matplotlib.pyplot as plt

# Specify request parameters (as strings)
token = 'aa31874e86fb42d9b2ea6b293f1bb004' # use your own token

#test
# Create API query string
args = {
    'start':'202106010000',
    'end':'202106070000',
    'obtimezone':'UTC',
    'vars':'precip_accum_one_hour',
    'stids': 'ksea',
        #['C3BCC','C3BVS','C3CAT','C3DLA','C3DRW','C3FRC','C3GPO','C3HDC','C3HRD','C3NBB','C3NCM','C3POR','C3PVN','C3SKI','C3SKY','C3SOD','C3WDG','C3WPO'],
    'units':'precip|in',
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
airT = responseDict['STATION'][0]['OBSERVATIONS'][args['vars']+'_set_1']
ksea = pd.Series(airT,index=pd.to_datetime(dateTime))

# Retain only the hourly observations
ksea = ksea.where(ksea.index.minute == 53).dropna()

# Plotting code
fig,ax = plt.subplots()
ksea.plot(ax=ax)
plt.title(args['stids'] + ' Hourly Precip (mm)')

# Clean up the plot a bit
from matplotlib.dates import DayLocator, DateFormatter
ax.xaxis.set_major_locator(DayLocator())
ax.xaxis.set_major_formatter(DateFormatter('%d-%b'))
ax.set_ylabel(args['units'])
ax.grid()

plt.savefig('temp2.png')