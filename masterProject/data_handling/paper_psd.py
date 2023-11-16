#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 15:28:28 2019

@author: gabgus
"""

import h5py
import numpy as np, matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
path = '/scratch/gabgus/fluent/henriks/noWall/planes/'


direction= 'x'
x=40
z=20

with h5py.File(path + 'plane_data.h5','r') as f:
    new_last = f.attrs['last_no'][-1]
    N=new_last
    data = f['high/half/{}v'.format(direction)][0:N,x,z] #- f['high/half/xv'][0:N,x,z].mean(axis=0)
    t = f['time'][0:N]
#%%

data = data

plt.figure()
plt.plot(t,data)
plt.xlabel('Time  [s]')
plt.ylabel('x-velocity [m/s]')
plt.title('x-velocity over time at a single point.')
plt.savefig('paper/characterisation/xvel.pdf',bbox_inches ='tight')

n= int(N/2)
ps=np.abs(np.fft.fft(data))**2
freqs=np.fft.fftfreq(data.size,5e-5)
idx=np.argsort(freqs)

plt.figure()
plt.semilogy(freqs[0:n],ps[0:n])
plt.xlim(-.1,100)
plt.title('PSD for the single point x-velocity')
plt.xlabel('Frequency  [Hz]')
plt.ylabel('Power density $[(m/s)^2 / Hz]$')
plt.ylim(0.1,1e9)
plt.savefig('paper/characterisation/psd.pdf',bbox_inches ='tight')


Etot=np.trapz(ps[0:n],freqs[0:n])
array=[]
for i in range(n):
    array.append(np.trapz(ps[0:i],freqs[0:i])/Etot)

for zoom in [(20,0.85,'zoom_out'),(2,0.01,'zoom_in')]:
    plt.figure()
    plt.plot(freqs[0:n],array,label = 'PSD integral')  
    plt.plot([freqs[0],freqs[n]],[1,1],label = '100%')
    plt.plot([freqs[0],freqs[n]],[0.99,0.99],'k--',alpha=0.2,label = '99%')
    plt.xlim(0,zoom[0])
    plt.ylim(zoom[1],1.001)
    plt.grid()
    plt.title('Integrated power spectral density for case 9')
    plt.xlabel('Frequency  [Hz]')
    plt.ylabel('Fraction of total energy')
    plt.legend(loc='lower right')
    plt.savefig('paper/characterisation/psd_int_{}.pdf'.format(zoom[2]),bbox_inches ='tight')
#%%

fig = plt.figure(figsize=(10,6)) # create a new figure with a default 111 subplot


ax = fig.add_subplot(2,1,1) 
axhigh = fig.add_subplot(2,2,4)
axlow = fig.add_subplot(2,2,3)


ax.plot(freqs[0:n],array,label = 'PSD integral')  
#ax.plot([freqs[0],freqs[n]],[1,1],label = '100%')
ax.plot([freqs[0],freqs[n]],[0.99,0.99],'k--',alpha=0.3,label = '99%')
ax.plot([0.25,0.25],[0,1.5],'r--',alpha=0.4,label = '0.25 Hz')

ax.set_xlim(0,20)   
ax.set_ylim(0,1.05)    
ax.legend(loc='lower right')
ax.set_xlabel('Frequency [Hz]')    
ax.set_ylabel('Fraction of total energy')
ax.grid(linestyle='--',alpha=0.5)
ax.set_title('Integrated and normalised PSD for case 4')



axhigh.plot(freqs[0:n],array,label = 'PSD integral')  
#axhigh.plot([freqs[0],freqs[n]],[1,1],label = '100%')
axhigh.plot([freqs[0],freqs[n]],[0.99,0.99],'k--',alpha=0.3,label = '99%')
axhigh.set_xlim(0,20)
axhigh.set_ylim(0.8,1.001)    
axhigh.legend()
axhigh.grid(linestyle='--',alpha=0.5)    
axhigh.set_xlabel('Frequency [Hz]')    
axhigh.set_title('Upper portion')
#axhigh.set_ylabel('Fraction of total energy')


axlow.plot(freqs[0:n],array,label = 'PSD integral')  
axlow.plot([0.25,0.25],[0,1.0],'r--',alpha=0.4,label = '0.25 Hz')
axlow.set_xlim(0,2)
axlow.set_ylim(0.1,0.85)    
axlow.legend()
axlow.grid(linestyle='--',alpha=0.5)
axlow.set_xlabel('Frequency [Hz]')    
axlow.set_title('Lower portion')
axlow.set_ylabel('Fraction of total energy')

    
fig.subplots_adjust(hspace=0.4,wspace=0.3)  
    
    
    