# -*- coding: utf-8 -*-
"""
Created on Sun Aug  4 11:34:49 2024

@author: MalekP
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import pickle 
import numpy as np
import os
from os import listdir
from os.path import isfile, join
import re


def combine_figures(paper_date,outdir):
        
    #figure consolidation
    
    #pull figure data associated with first instance of hourly data for the hours that exhibit the top 4 rainfall rate figures
    mypath = outdir
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    year = [x.zfill(2) for x in list(map(str, list(paper_date.year)))]
    month = [x.zfill(2) for x in list(map(str, list(paper_date.month)))]
    day = [x.zfill(2) for x in list(map(str, list(paper_date.day)))]
    hour =[x.zfill(2) for x in list(map(str, list(paper_date.hour)))]
    
    dates_to_combine = [i[0]+i[1]+i[2]+i[3] for i in zip(year,month,day, hour)]
    
    date_files = [re.sub(r'\D', '', i)[0:10] for i in onlyfiles]
    
    
    index_selection = [([idx for idx, val in enumerate(date_files) if val == sub] if sub in date_files else [None])
      for sub in dates_to_combine]
    
    index_first = [item[0] for item in index_selection]
    
    
    files_to_combine  = [onlyfiles[i] for i in index_first] 
    
    paths_to_combine = [mypath + '\\' + s for s in files_to_combine]
    
    figures = {}
    for file in paths_to_combine:
        print(file)
        file_to_open = open(file,'rb')
    
        # dump information to that file
        figures[file]=pickle.load(file_to_open)
        print('Loaded')
        # close the file
        #file_to_open.close()
    
    # #figure 2
    # file = open('G:\\NCFR Thesis\\NCFR_Thesis\\combined_kove_2017220_2017240\\KDAX20170203_071101_V06_paper.pickle','rb')
    
    # # dump information to that file
    # tt2=pickle.load(file)
    
    # # close the file
    # file.close()
    
    # #figure 3
    # file = open('G:\\NCFR Thesis\\NCFR_Thesis\\combined_kove_2017220_2017240\\KDAX20170203_120525_V06_paper.pickle','rb')
    
    # # dump information to that file
    # tt3=pickle.load(file)
    
    # # close the file
    # file.close()
    
    # #figure 4
    # file = open('G:\\NCFR Thesis\\NCFR_Thesis\\combined_kove_2017220_2017240\\KDAX20170203_172257_V06_paper.pickle','rb')
    
    # # dump information to that file
    # tt4=pickle.load(file)
    
    # # close the file
    # file.close()
    file_to_open.close()
    
    plt.close('all')
    
    
    backend = mpl.get_backend()
    mpl.use('agg')
    

    
    
    # sb1 = list(figures.values())[0].colorbar
    # sb1.remove()
    c1 = list(figures.values())[0]
    # sb2 = list(figures.values())[1].colorbar()
    c2 = list(figures.values())[1]
    #sb3 = list(figures.values())[2].colorbar()
    c3 = list(figures.values())[2]
    c4 = list(figures.values())[3]
    for i in c1.axes:
        if i.get_label() == "<colorbar>" or i.get_ylabel() == "NWSReflectivity" :
            i.remove()
            
    for i in c2.axes:
        if i.get_label() == "<colorbar>" or i.get_ylabel() == "NWSReflectivity" :
            i.remove()
            
    for i in c3.axes:
        if i.get_label() == "<colorbar>" or i.get_ylabel() == "NWSReflectivity" :
            i.remove()
    
    
    c1.canvas.draw()
    c2.canvas.draw()
    c3.canvas.draw()
    c4.canvas.draw()
    dpi=300
    a1 = np.array(c1.canvas.buffer_rgba())
    a2 = np.array(c2.canvas.buffer_rgba())
    a3 = np.array(c3.canvas.buffer_rgba())
    a4 = np.array(c4.canvas.buffer_rgba())
    a = np.hstack((a1,a2,a3,a4))
    
    mpl.use(backend)
    fig,ax = plt.subplots(figsize=(40, 35))
    fig.subplots_adjust(0, 0, 1, 1)
    
    font = {'family' : 'normal',
        'size'   : 22}

    plt.rc('font', **font)
    
    ax.set_axis_off()
    ax.matshow(a)
    
    fig.savefig(outdir+"combined_figure_for_paper.jpeg",dpi=300)