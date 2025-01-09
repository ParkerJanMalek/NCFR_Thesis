# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 19:22:43 2024

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
import imageio


png_dir = "G://NCFR Thesis//NCFR_Thesis//combined_kove_2017268_2017268//"

pic_dir = "G://NCFR Thesis//NCFR_Thesis//"


pic_dir=['G://NCFR Thesis//NCFR_Thesis//combined_kove_201721613_201721613',
'G://NCFR Thesis//NCFR_Thesis//combined_kove_201721718_201721818',
'G://NCFR Thesis//NCFR_Thesis//combined_kove_20172194_20172214',
'G://NCFR Thesis//NCFR_Thesis//combined_kove_2017220_2017240',
'G://NCFR Thesis//NCFR_Thesis//combined_kove_20172615_20172715',
'G://NCFR Thesis//NCFR_Thesis//combined_kove_2017268_2017268',
'G://NCFR Thesis//NCFR_Thesis//combined_kove_20172822_20172822',
'G://NCFR Thesis//NCFR_Thesis//combined_kove_2017290_20172100']

#pic_dir = "G://NCFR Thesis//NCFR_Thesis//combined_kove_2017268_2017268"
for i in pic_dir:
    print(i)
    images = []
    for file_name in sorted(os.listdir(i)):
        if file_name.endswith('.png'):
            print(file_name)
            file_path = os.path.join(i, file_name)
            images.append(imageio.imread(file_path))
    
    # Make it pause at the end so that the viewers can ponder
    # for _ in range(5):
    #     images.append(imageio.imread(os.path.join(pic_dir,i)))
      
    imageio.mimsave(i+"//"+'combined_kove.gif', images, loop=0)
#images=[]