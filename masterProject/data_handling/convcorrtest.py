#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 15:30:07 2019

@author: gabgus
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 14:14:13 2019

@author: gabgus
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
import time
# -----------------------------------------------------------------------------
path = '/scratch/gabgus/fluent/henriks/noWall/planes/'


## interpolate the velocity of interest onto the created mesh
meshname='half'
elevation='high'
with h5py.File(path + 'plane_data.h5','r') as f:    
   maxN = f.attrs['last_no'][-1]
   xi = f['mesh'][meshname]['xi'][:]
   zi = f['mesh'][meshname]['zi'][:]
   vi = f[elevation][meshname]['xv'][0:maxN]
   t = f['time'][0:maxN]
#%%

def RL(kdot,v_inst,vbar,vdot):
    top=0
    bot=0
    ksum=0
    
    if type(kdot) == 'int':
        kmax=kdot
    else:
        kmax=kdot[-1]
        
    for k in range(int(len(v_inst) - kmax)):
        top = top + (v_inst[k] - vbar)*(vinst[k+kdot] - vbar)
        bot  = bot+ (v_inst[k] - vbar)**2
        ksum +=1
    
    top=top/ksum
    bot=bot/ksum
    
    return top/bot


def origwrap(x,lags):
    vbar = vinst.mean(axis=0)

    vdot = vinst-vbar
    return RL(lags,x,vbar,vdot)


def autocorr1(x,lags):
    '''np.corrcoef, partial'''

    corr=[1. if l==0 else np.corrcoef(x[l:],x[:-l])[0][1] for l in lags]
    return np.array(corr)

def autocorr2(x,lags):
    '''manually compute, non partial'''

    mean=np.mean(x)
    var=np.var(x)
    xp=x-mean
    corr=[1. if l==0 else np.sum(xp[l:]*xp[:-l])/len(x)/var for l in lags]

    return np.array(corr)

def autocorr3(x,lags):
    '''fft, pad 0s, non partial'''

    n=len(x)
    # pad 0s to 2n-1
    ext_size=2*n-1
    # nearest power of 2
    fsize=2**np.ceil(np.log2(ext_size)).astype('int')

    xp=x-np.mean(x)
    var=np.var(x)

    # do fft and ifft
    cf=np.fft.fft(xp,fsize)
    sf=cf.conjugate()*cf
    corr=np.fft.ifft(sf).real
    corr=corr/var/n

    return corr[:len(lags)]

def autocorr5(x,lags):
    '''np.correlate, non partial'''
    mean=x.mean()
    var=np.var(x)
    xp=x-mean
    corr=np.correlate(xp,xp,'full')[len(x)-1:]/var/len(x)

    return corr[:len(lags)]








n=9
frac=0.4
vinst = vi[::n,7,25]

N=int(maxN*frac)#10050
karr=np.arange(int(N/n))

fig,ax=plt.subplots()

runtimes = []
for funcii, labelii in zip([origwrap, autocorr2, autocorr3,
    autocorr5], ['original function', 'manual, non-partial',
        'fft, pad 0s, non-partial',
        'np.correlate, non-partial']):
    start = time.time()
    cii=funcii(vinst,karr)
    runtimes.append(time.time() - start) 
    #print(labelii)
    #print(cii)
    ax.plot(karr*n*5e-5,cii,label=labelii + ' {:.4f} s'.format(runtimes[-1]) )

plt.plot([0,4],[0,0],'--k')
ax.set_xlabel('lag')
ax.set_ylabel('correlation coefficient')
ax.legend()
plt.show()






