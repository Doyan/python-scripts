# -*- coding: utf-8 -*-
"""
Created on Sun Jun  3 23:06:52 2018

@author: Gabriel
"""

import numpy as np, matplotlib.pyplot as plt
import pandas as pd
# -----------------------------------------------------------------------------
path = './datafiles/'

fnam ='Qbed.csv'


data = pd.read_csv(path+fnam,sep=',',skiprows=1, names=['time','q'] )

#%%

N = data.shape[0]

qframe = data.q.iloc[0:N]
tframe = data.time.iloc[0:N]

t = np.array(tframe)
q = np.array(qframe)

qsum =qframe.cumsum()
qavg= qframe.expanding().mean()

qrev = qframe.sort_index(ascending=False, axis=0)
trev = tframe.sort_index(ascending=False, axis=0)

qavgr = qrev.expanding().mean()

avgrdiff = qavgr.diff() 

plt.figure(0, figsize=(4,4))
plt.plot(t,qframe)
plt.xlabel('Time [s]')
plt.ylabel('Net energy flow [W]')
plt.title('Raw monitor signal')
plt.xlim(7.9,24.2)
plt.grid()
plt.savefig('raw.pdf', bbox_inches='tight')


plt.figure(1,figsize=(4,4))
plt.plot(t,qavg,label='Expanding average')
plt.plot(trev,qavgr,label='Reverse expanding average')
plt.grid()
plt.legend(loc='best')
plt.xlabel('Time [s]')
plt.ylabel('Average net energy flow [W]')
plt.title('Expanding time average')
plt.xlim(7.9,24.2)
plt.savefig('cumavg.pdf', bbox_inches='tight')


plt.figure(2)
plt.plot(t,qavg,label='Expanding average')
plt.plot(trev,qavgr,label='Reverse expanding average')
plt.plot([7,25],[555000,555000],'--',color=(0.5,0.5,0.5),label='Expected lower bound')
plt.plot([7,25],[2777000,2777000],'-.',color=(0.5,0.5,0.5),label='Expected higher bound')

plt.xlabel('Time [s]')
plt.ylabel('Average net energy flow [W]')
plt.title('Expanding time average')

plt.xlim(8.01,24.2)
plt.ylim(-0.8e7,0.8e7)
plt.grid()
plt.legend(loc='best')

plt.savefig('cumavgzoom.pdf', bbox_inches='tight')


plt.figure(3)
plt.plot(t,qsum)
plt.grid()
plt.xlabel('Time [s]')
plt.ylabel('Cumulative energy flow [W]')
plt.title('Net accumulation')

print(qavg.iloc[-1] * 1.4 / (0.3*1.2*(60)))


