# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 18:10:10 2023

Pulls data based from MRMS database based on start and end dates

@author: parke
"""
import datetime
import wget

start_date = datetime.datetime(2020,1,1,1)
end_date = datetime.datetime(2020,1,31,24)


from datetime import date, timedelta

def daterange(start_date, end_date):
    for n in range(int((end_date - start_date).days)):
        yield start_date + timedelta(n)
        
for single_date in daterange(start_date, end_date):
    url = 'https://mtarchive.geol.iastate.edu/'+str(single_date.year)+'/'+str(single_date.month)+'/'+str(single_date.day)+'/mrms/ncep/'
    print(url)
        
