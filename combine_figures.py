# -*- coding: utf-8 -*-
"""
Created on Sun Aug  4 11:34:49 2024

@author: MalekP
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import pickle 
import numpy as np

#figure 1
file = open('G:\\NCFR Thesis\\NCFR_Thesis\\combined_kove_2017220_2017240\\KDAX20170202_120155_V06_paper.pickle','rb')

# dump information to that file
tt=pickle.load(file)

# close the file
file.close()

#figure 2
file = open('G:\\NCFR Thesis\\NCFR_Thesis\\combined_kove_2017220_2017240\\KDAX20170203_071101_V06_paper.pickle','rb')

# dump information to that file
tt2=pickle.load(file)

# close the file
file.close()

#figure 3
file = open('G:\\NCFR Thesis\\NCFR_Thesis\\combined_kove_2017220_2017240\\KDAX20170203_120525_V06_paper.pickle','rb')

# dump information to that file
tt3=pickle.load(file)

# close the file
file.close()

#figure 4
file = open('G:\\NCFR Thesis\\NCFR_Thesis\\combined_kove_2017220_2017240\\KDAX20170203_172257_V06_paper.pickle','rb')

# dump information to that file
tt4=pickle.load(file)

# close the file
file.close()


plt.close('all')


backend = mpl.get_backend()
mpl.use('agg')

c1 = tt.canvas
c2 = tt2.canvas
c3 = tt3.canvas
c4 = tt4.canvas

c1.draw()
c2.draw()
c3.draw()
c4.draw()
dpi=300
a1 = np.array(c1.buffer_rgba())
a2 = np.array(c2.buffer_rgba())
a3 = np.array(c3.buffer_rgba())
a4 = np.array(c4.buffer_rgba())
a = np.hstack((a1,a2,a3,a4))

mpl.use(backend)
fig,ax = plt.subplots(figsize=(40, 35))
fig.subplots_adjust(0, 0, 1, 1)
ax.set_axis_off()
ax.matshow(a)

fig.savefig("test1.jpeg",dpi=300)