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
from PIL import Image


def combine_figures(paper_date,outdir):
        
    #figure consolidation
    
    #pull figure data associated with first instance of hourly data for the hours that exhibit the top 4 rainfall rate figures
    mypath = outdir+'paper\\'
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
    
    # figures = {}
    # for file in paths_to_combine:

        
        
        
    #     print(file)
    #     file_to_open = open(file,'rb')
    
    #     # dump information to that file
    #     figures[file]=pickle.load(file_to_open)
    #     print('Loaded')
    #     # close the file
    #     file_to_open.close()
    
    # # #figure 2
    # # file = open('G:\\NCFR Thesis\\NCFR_Thesis\\combined_kove_2017220_2017240\\KDAX20170203_071101_V06_paper.pickle','rb')
    
    # # # dump information to that file
    # # tt2=pickle.load(file)
    
    # # # close the file
    # # file.close()
    
    # # #figure 3
    # # file = open('G:\\NCFR Thesis\\NCFR_Thesis\\combined_kove_2017220_2017240\\KDAX20170203_120525_V06_paper.pickle','rb')
    
    # # # dump information to that file
    # # tt3=pickle.load(file)
    
    # # # close the file
    # # file.close()
    
    # # #figure 4
    # # file = open('G:\\NCFR Thesis\\NCFR_Thesis\\combined_kove_2017220_2017240\\KDAX20170203_172257_V06_paper.pickle','rb')
    
    # # # dump information to that file
    # # tt4=pickle.load(file)
    
    # # # close the file
    # # file.close()
    # #file_to_open.close()
    
    # plt.close('all')
    
    
    # backend = mpl.get_backend()
    # mpl.use('agg')
    

    # # for fig in figures.values():
    # #     if fig.get_ylabel() == "kg $\\mathregular{m^{-1}}$ $\\mathregular{s^{-1}}$" or fig.get_ylabel() == "NWSReflectivity" or  fig.get_ylabel() == "hPa" or fig.get_ylabel() == "Degrees/hr" or fig.get_ylabel() == "K":
    # #         fig.remove()
    # #     fig.set_figwidth(120)  # Increase width of each figure
    # #     fig.set_figheight(120) 
    
    # # sb1 = list(figures.values())[0].colorbar
    # # sb1.remove()
    # c1 = list(figures.values())[0]
    # c2 = list(figures.values())[1]
    # #sb3 = list(figures.values())[2].colorbar()
    # c3 = list(figures.values())[2]
    # c4 = list(figures.values())[3]
    # c1.tight_layout(pad=0)
    # c2.tight_layout(pad=0)
    # c3.tight_layout(pad=0)
    # c4.tight_layout(pad=0)
    # for i in c1.axes:
    #     #i.set_anchor('E')
    #     if i.get_ylabel() == "kg $\\mathregular{m^{-1}}$ $\\mathregular{s^{-1}}$" or i.get_ylabel() == "NWSReflectivity" or  i.get_ylabel() == "hPa" or i.get_ylabel() == "Degrees/hr" or i.get_ylabel() == "K":
    #         i.remove()
            
    # for i in c2.axes:
    #     #i.set_anchor('E')
    #     if i.get_ylabel() == "kg $\\mathregular{m^{-1}}$ $\\mathregular{s^{-1}}$" or i.get_ylabel() == "NWSReflectivity" or  i.get_ylabel() == "hPa" or i.get_ylabel() == "Degrees/hr" or i.get_ylabel() == "K":
    #         i.remove()
            
    # for i in c3.axes:
    #     #i.set_anchor('E')
    #     if i.get_ylabel() == "kg $\\mathregular{m^{-1}}$ $\\mathregular{s^{-1}}$" or i.get_ylabel() == "NWSReflectivity" or  i.get_ylabel() == "hPa" or i.get_ylabel() == "Degrees/hr" or i.get_ylabel() == "K":
    #         i.remove()
    # # for i in c4.axes:
    # #     #i.set_anchor('E')
    # #     if i.get_label() == "<colorbar>" or i.get_ylabel() == "NWSReflectivity" :
    # #         i.remove()
    
    
   
    
    # # figure_arrays = []
    # # for fig in figures.values():
    # #     if i.get_ylabel() == "kg $\\mathregular{m^{-1}}$ $\\mathregular{s^{-1}}$" or i.get_ylabel() == "NWSReflectivity" or  i.get_ylabel() == "hPa" or i.get_ylabel() == "Degrees/hr" or i.get_ylabel() == "K":
    # #         i.remove()
    # #     fig.set_figwidth(120)  # Increase width of each figure
    # #     fig.set_figheight(120)  # Increase height of each figure
    # #     fig.tight_layout(pad=0)  # Minimize internal padding
    # #     fig.canvas.draw()  # Draw the figure to update the canvas
    # #     figure_arrays.append(np.array(fig.canvas.buffer_rgba()))  # Convert canvas to RGBA array
    
    # # # Combine the arrays horizontally
    # # combined_array = np.hstack(figure_arrays)
    
    # # # Calculate the aspect ratio
    # # aspect_ratio = combined_array.shape[0] / combined_array.shape[1]
    
    # # # Create a new figure to display the combined image
    # # fig, ax = plt.subplots(figsize=(60, 160))  # Adjust figsize to fit combined image
    # # fig.subplots_adjust(left=0.05, bottom=0.07, right=0.95, top=0.95, wspace=0, hspace=0)
    # # ax.set_axis_off()  # Turn off the axis
    # # ax.imshow(combined_array)  # Display the combined image and stretch it to fill
    
    # # fig.savefig(outdir+"combined_figure_for_paper.jpeg",bbox_inches='tight',dpi=1300, pad_inches=0)
    
    # # f, axarr = plt.subplots(1,4)
    # # axarr[0].imshow(c1.canvas.buffer_rgba())
    # # axarr[1].imshow(c2.canvas.buffer_rgba())
    # # axarr[2].imshow(c3.canvas.buffer_rgba())
    # # axarr[3].imshow(c4.canvas.buffer_rgba())
    
    
    # dpi=1300
    # c1.canvas.draw()
    
    # a1 = np.array(c1.canvas.buffer_rgba())
    
    # c2.canvas.draw()
    # a2 = np.array(c2.canvas.buffer_rgba())
    
    # c3.canvas.draw()
    # a3 = np.array(c3.canvas.buffer_rgba())
    
    # c4.canvas.draw()
    # a4 = np.array(c4.canvas.buffer_rgba())
    # a = np.hstack((a1,a2,a3,a4))
    
    # mpl.use(backend)
    # fig,ax = plt.subplots(figsize=(60, 60),frameon=False)
    # ax.matshow(a)
    # ax.axis('off')
    # fig.tight_layout(pad=0)
    # fig.subplots_adjust(left=0,
    #                 bottom=0,
    #                 right=1,
    #                 top=1)
    # ax.margins(0,0)
    
    # ax.set_axis_off()
    # fig.subplots_adjust(bottom=0, top=1, left=0, right=1, hspace=0, wspace=0)
    # font = {'family' : 'normal',
    #     'size'   : 22}

    # plt.rc('font', **font)
    
    #ax.set_axis_off()
    
    # Open the images
    if len(paths_to_combine) > 1:
        img1 = Image.open(paths_to_combine[0])
        img2 = Image.open(paths_to_combine[1])
        img3 = Image.open(paths_to_combine[2])
        img4 = Image.open(paths_to_combine[3])
        
        # Find the total width and height (assuming they all have the same height)
        total_width = img1.width + img2.width + img3.width + img4.width
        max_height = max(img1.height, img2.height, img3.height, img4.height)
    
        # Create a new blank image with the combined size
        combined_img = Image.new('RGB', (total_width, max_height))
    
        # Paste the images side by side
        combined_img.paste(img1, (0, 0))
        combined_img.paste(img2, (img1.width, 0))
        combined_img.paste(img3, (img1.width + img2.width, 0))
        combined_img.paste(img4, (img1.width + img2.width + img3.width, 0))
    
        # Show or save the combined image
        #combined_img.show()
        combined_img.save(outdir+'paper\\combined_figures.png')
    
   
    
   # fig.savefig(outdir+"combined_figure_for_papert.jpg",bbox_inches='tight',dpi=300, pad_inches=0)