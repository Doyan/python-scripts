#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 16:28:03 2019

Makes a certain bar chart for  a paper comparing different methods of 
estimating lateral dispersion in a fluidised bed.

@author: gabgus
"""

import matplotlib.pyplot as plt
import numpy as np

# -----------------------------------------------------------------------------
exp_corr = 1.5625
# case 1, 6, 8, 9 10 
D_fits = [0.0016, 0.0098, 0.0226, 0.00840,0.01495]
k_fits = [0.00075, 0.0023, 0.00655, 0.00585, 0.0079]
q_signal = [0.00243839, 0.006174958, 0.019210092, 0.009364, 0.016194 ]
exp = [0.0025,0.0025,0.0145,0.0025,0.0145]
vdotmethod=[0,0,0,0.0058,0]
ticklabels = ['nozzles', 'pp', 'sette', 'pp large', 'sette large']

n_groups = 5

# rearrange because i'm lazy
for array in [D_fits,k_fits,q_signal,exp,vdotmethod, ticklabels]:
    array[2], array[3] = array[3], array[2]


# -----------------------------------------------------------------------------
    
fig, ax = plt.subplots()

index = np.arange(n_groups)
bar_width = 0.17

opacity = 0.4
error_config = {'ecolor': '0.3'}

index[-2:] = index[-2:] + 1
ax.zorder = 0
ax.grid(linestyle = '--', alpha=0.5,zorder=0,axis='y')


rects1 = ax.bar(index+ bar_width*1, D_fits, bar_width,
                alpha=opacity, color='b',
                label='M2',zorder=40)


rects2 = ax.bar(index + bar_width*2, q_signal, bar_width,
                alpha=opacity, color='g',
                label='M3',zorder=41)


rects3 = ax.bar(index + bar_width*3, k_fits, bar_width,
                alpha=opacity, color='r',
                label='M1',zorder=42)

rects4 = ax.bar(index + bar_width*4, vdotmethod, bar_width,
                alpha=opacity, color='y',
                label='M4',zorder=43)

#rects5 = ax.bar(index + bar_width*0, exp, bar_width,
#                alpha=opacity, color='k',
#                label='experiment',zorder=44)


# -----------------------------------------------------------------------------


ax.axvline(index[-2]-bar_width*4, linestyle='--',color='gray',alpha=0.3,linewidth=1 )

#ax.set_xlabel('Case')
ax.set_ylabel('$D_L$ $[m^2/s]$')
ax.set_title('Calculated dispersion coefficients for the free bed')
ax.set_xticks(index + bar_width*4 / 2)
ax.set_xticklabels(ticklabels)
ax.legend(loc='upper left')


fig.tight_layout()
plt.show()
fig.savefig('paper/bar_chart.pdf')





