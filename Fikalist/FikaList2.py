# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 16:21:07 2018

Script to generate a "fairly" randomized list to tell who'se time it is 
to provide fika every day. 


@author: Gabriel
"""

from openpyxl import load_workbook

# ----------------------------------------------------------------------------

# Path to fikalist file
fikalistPath = './FikaListDataBase.xlsx' 

# Days for person to stay off the active drawing pool after having fika once
sTime = 10


# Read fikalist file, 
# assumes 1 header row and columns ordered: Name, fikaScore, kitchenScore 
def readDataBase(fikalistPath):  
    # Loading Workbook and getting active sheet
    wb = load_workbook(fikalistPath)

    sheet = wb.active

    cells = sheet[sheet.dimensions]
    
    # Parse workbook and sort columns into lists
    names = []
    fikaScore = []
    kitchenScore = []
    for row in cells:
        names.append(row[0].value)
        fikaScore.append(row[1].value)
        kitchenScore.append(row[2].value)
    
    # Remove header lines from lists
    namestring = names.pop(0)
    fikastring = fikaScore.pop(0)
    kitchenstring = kitchenScore.pop(0)
    
    # print pop capture to help identify if columns have been messed with
    print('sorted "' + namestring + '" into "names"')
    print('sorted "' + fikastring + '" into "fikaScore"')
    print('sorted "' + kitchenstring + '" into "kitchenScore"\n')

    return names, fikaScore, kitchenScore


def makeFikaSequence(nameList,scoreList,days):
    
#
#names, fikaScore, kitchenScore = readDataBase(fikalistPath)
#
#for i in range(len(names)):
#    print('{0:8}\t\t{1:8}\t{2:8}'.format(names[i],fikaScore[i],kitchenScore[i]))





