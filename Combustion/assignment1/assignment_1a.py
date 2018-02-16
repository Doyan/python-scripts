# -*- coding: utf-8 -*-
"""
Created on Fri Sep 08 13:32:07 2017

@author: Gabriel
"""

import cantera as ct, numpy as np

Rm=ct.gas_constant/1000
Tref=298
Nu=np.array([1, 3.5, 13.16, 0, 0])
Nb=np.array([0, 0, 13.16, 2, 3])

Href=np.array([-84667, 0, 0, -393546, -241845])
Cp=np.array([0, 0, 35.099, 58.836, 48.035])

Cv= Cp-Rm

Tb= Tref + (sum(Nu*Href)- sum(Nb*Href))/sum(Nb*Cv)

print Tb

