#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 14:23:33 2018

@author: gabgus
"""

import numpy as np

import matplotlib.pyplot as plt

#----------------------------------------------------------------------------------

# high mesh
TG_avg_h = [765.1691502546852, 773.7080757144821]
TC_max_h = [872.8852783203125, 882.412255859375]
TG_min_h = [746.1984497070313, 754.5776611328125]
Cells_h =  4.821680e+05

# mid mesh

TG_avg_m = [763.3680105652624, 771.8211752373446]
TC_max_m = [868.464501953125, 877.8099609375]
TG_min_m = [745.0793701171875, 753.389794921875]
Cells_m = 1.200910e+05

# low mesh

TG_avg_l = [7.616326e+02, 7.684362e+02]
TC_max_l = [8.675207e+02, 8.774951e+02]
TG_min_l = [7.430160e+02, 7.492616e+02]
Cells_l = 6.317400e+04


x = [Cells_l, Cells_m, Cells_h]


TG_avg = []
dt = []
for i in range(len(TG_avg_h)):
    TG_avg.append([TG_avg_l[i], TG_avg_m[i], TG_avg_h[i] ])
    dt.append([TC_max_l[i] - TG_min_l[i], TC_max_m[i]- TG_min_m[i], TC_max_h[i] - TG_min_h[i]])

plt.figure(1,figsize=(8,5))
plt.plot(x,TG_avg[0],label='case 0')
plt.plot(x,TG_avg[1], label = 'case 1')
plt.plot([Cells_m,Cells_m],[750,780],'k--', label = 'Used mesh size')
plt.ylim(761,774)

plt.xlabel('Cells')
plt.ylabel('$T_{G,avg} \;\; [\degree C]$')
plt.title('Mesh comparison')
plt.legend(loc='best')