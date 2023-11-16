#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 11:18:16 2019
paper D-fit example
@author: gabgus
"""

import numpy as np
import matplotlib.pyplot as plt
import comparator as C

# ---------------------------------------------------------------------------

ybar, gx, t = C.get1ddata(10,'uds')

kbar = C.fourier1d(gx,t,0.014)


L = (gx[-1] + gx[1]) - 2*gx[0]

sNo = 0


fig = plt.figure(1,(8,6))

axsim = fig.add_subplot(212)
axana = fig.add_subplot(211,sharex=axsim)

sampleNos = [0,2,10,50,100,150]
for sNo in sampleNos:
    axsim.plot(gx,ybar[sNo],'.-',label=f't={t[sNo]}s')
    axana.plot(gx,kbar[sNo],'.-',label=f't={t[sNo]}s')


for ax in [axana,axsim]:
    ax.grid()
    ax.set_ylabel('species fraction')
    ax.plot([L/2,L/2],[0,1],'k--',alpha=0.3)
    ax.plot([0,L],[0.5,0.5],'k--',alpha=0.3)
    ax.set_xlim(gx[0],gx[-1])
    ax.legend(loc='lower right',ncol=2)

axsim.set_xlabel('x-coordinate [m]')

plt.setp(axana.get_xticklabels(),visible=False)
fig.suptitle('Analytic vs Simulated concentration profiles')

fig.text(0.2,0.84,'Analytic')
fig.text(0.2,0.42,'Simulated')

fig.savefig('paper/Dfit_example.pdf',bbox_inches='tight')





