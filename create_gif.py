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
images = []
for i in sorted(os.listdir(pic_dir)):
    if i.startswith('combined_kove'):
        print(i)
        for file_name in sorted(os.listdir(pic_dir+i)):
            if file_name.endswith('.png'):
                file_path = os.path.join(pic_dir,i, file_name)
                images.append(imageio.imread(file_path))

        # Make it pause at the end so that the viewers can ponder
        # for _ in range(5):
        #     images.append(imageio.imread(os.path.join(pic_dir,i)))

        imageio.mimsave(pic_dir+i+"//"+i+'_presentation.gif', images)