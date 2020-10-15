# -*- coding: utf-8 -*-
"""
Created on Sun May  5 14:57:02 2019
DND dice
@author: Gabriel Gustafsson
"""

import numpy as np



# 5
def d(N,n=1):
    
    ans = []
    for n in range(n):
        ans.append(np.random.randint(1,N+1))
    return np.array(ans)

def d4():
    
    return np.random.randint(1,4)

def d6():
    return np.random.randint(1,6)

def d10():
    return np.random.randint(1,10)

def d20() :
    return np.random.randint(1,20)







"""

1
------
1
--
10
--
13
--
12
--






"""


