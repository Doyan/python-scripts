import comparator as C
import numpy as np
import matplotlib.pyplot as plt

mtemp = C.mtemp

imagefolder = 'paper/profile_symmetries/'

knumbers = [4000]*6  + [8000]*2 + [25000] + [10000, 10000]
dnumbers = [0.0025]*6 + [0.008]*2 + [0.025] + [0.009,0.015]

caseNo=1
scalar = 'uds'

frac=1



for caseNo in [1,2,3,4,6,7,8,9,10]:
    for scalar  in ['temp','uds']:
        for frac in [1.0,0.4,0.1]:
            
            knumber = knumbers[caseNo]
            dnumber = dnumbers[caseNo]
            x1 = C.Klimits['x1'][caseNo]
            x0 = C.Klimits['x0'][caseNo]
            
            ybar, gx,t  =  C.get1ddata(caseNo,scalar=scalar)
            
            if scalar == 'temp':
                xslice = slice(x0,x1)
                kbar = C.heateq1D(gx[xslice],t,knumber, vof_s=C.vof_s[caseNo],mtemp=mtemp[caseNo],Tcorr=0)    
                
            
            
            if scalar == 'uds':
                xslice = slice(0,500)
                kbar = C.fourier1d(gx,t,dnumber)
            
            fig, ax = plt.subplots(2,1,sharex=True)
            
            
            
            
            for sNo in range(int(len(ybar)*frac)):
                ax[0].plot(gx[xslice],kbar[sNo])
                ax[1].plot(gx[xslice],ybar[sNo][xslice])
                
                
            
            if scalar == 'temp':
                ax[0].text(x1-0.4,mtemp[caseNo]+25,'Analytic')
                ax[1].text(x1-0.4,mtemp[caseNo]+25,'Simulated')
                for axes in ax:
                    axes.grid()
                    axes.plot([gx[x0-2],gx[x1+2]],[mtemp[caseNo]]*2,'k--',alpha=0.5)
                    axes.plot([(gx[x1]-gx[x0])/2+gx[x0],(gx[x1]-gx[x0])/2+gx[x0]],[mtemp[caseNo]+25,mtemp[caseNo]-25],'k--',alpha=0.5)
                    axes.plot([gx[x0],gx[x1]],[mtemp[caseNo]+25,mtemp[caseNo]-25],'k--',alpha=0.3)
                    axes.set_ylabel('Temperature [K]')
                    axes.set_xlim(gx[x0-1],gx[x1])
                    
            if scalar == 'uds':
                for axes in ax:
                    axes.grid()
                    axes.plot([gx[0],gx[-1]],[0.5]*2,'k--',alpha=0.5)
                    axes.plot([(gx[-1]-gx[0])/2+gx[0],(gx[-1]-gx[0])/2+gx[0]],[0,1],'k--',alpha=0.5)
                    axes.set_ylabel('Species fraction')
            
            
            ax[1].set_xlabel('x-distance [m]')
            fig.suptitle('Comparison all profiles, analytic vs sim')
            fig.savefig(imagefolder + 'c{}_{}_frac{}.pdf'.format(caseNo,scalar,frac))

