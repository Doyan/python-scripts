# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 14:00:10 2018

@author: Gabriel
"""

import re
import numpy as np

file=open("file.txt",'r') 

floatpattern = re.compile("\d\.\d*")

i=0

span={}

span[i]=[]

row=file.readlines()

for line in row:      
    if floatpattern.match(line): 
        data=line.split()
        span[i].append(data)
        
    if not floatpattern.match(line):
        if not len(span[i]) == 0:
            i=i+1
            span[i]=[]

file.close()