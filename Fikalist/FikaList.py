

import xlsxwriter
import xlrd
import re
import random

def ReadDataBase(fileName, dataBaseSheetName, rowMin, rowMax):
    workbook = xlrd.open_workbook(fileName)
    sheet = workbook.sheet_by_name(dataBaseSheetName)
    names = []
    fikaCount = []
    kitchenCount = []
    for i in range(rowMin-1,rowMax):
        names.append(sheet.cell(i, 0).value)
        fikaCount.append(sheet.cell(i, 1).value)
        kitchenCount.append(sheet.cell(i, 2).value)

    return [names, fikaCount, kitchenCount]

def WriteDataBase(fileName, dataBaseSheetName, names, fikaCount, kitchenCount):
    workbook = xlsxwriter.Workbook(fileName)
    worksheet = workbook.add_worksheet(dataBaseSheetName)

    # Formatting of the excel sheet
    bold = workbook.add_format({'bold': True})
    worksheet.set_column(0, 0, 23)
    worksheet.set_column(1, 1, 10)
    worksheet.set_column(2, 2, 12)

    # Write header of database
    worksheet.write(0, 0, 'Name', bold)
    worksheet.write(0, 1, 'Fika count', bold)
    worksheet.write(0, 2, 'Kitchen count', bold)
    row = 1
    for i in range(len(names)):
        worksheet.write(row, 0, names[i])
        worksheet.write(row, 1, fikaCount[i])
        worksheet.write(row, 2, kitchenCount[i])
        row += 1

    workbook.close()

    return

def AdvanceDay(day, month, year):
    daysPerMont = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    day += 1
    if month == 2 and year%4 == 0:
        daysPerMont[1] = 29

    if(day>daysPerMont[month-1]):
        day = 1
        month += 1
        if month > 12:
            month = 1
            year +=1

    return [day, month, year]

def CreateKitchenList(kitchenListName, startingDate, names, kitchenCount, numberOfWeeks):
    months=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    dateMatch = re.search(r'([\d]+)/([\d]+)/([\d]+)', startingDate)
    if dateMatch:
        startingMonth = int(dateMatch.group(1))
        startingDay = int(dateMatch.group(2))
        startingYear = int(dateMatch.group(3))

    month = startingMonth
    day = startingDay
    year = startingYear

    listCompleted = False

    selectedPeople = []

    timesInThisList = []
    for i in range(len(names)):
        timesInThisList.append(0)

    while not listCompleted:
        # Search for the person with max fikaCount
        maxFikaCount = max(kitchenCount)
        minFikaCount = min(kitchenCount)

        if maxFikaCount == minFikaCount:
            maxFikaCount += 1

        # Get the persons with less fika than the maximum
        lessFikaPeopleIndex = []
        for i in range(len(kitchenCount)):
            if kitchenCount[i] < maxFikaCount:
                lessFikaPeopleIndex.append(i)

        while len(lessFikaPeopleIndex)>0:
            print('lessFikaPeopleIndex length = %f' % len(lessFikaPeopleIndex))
            selectedIndex = random.randint(0,len(lessFikaPeopleIndex)-1)
            selected = lessFikaPeopleIndex[selectedIndex]
            if timesInThisList[selected] == 0:
                timesInThisList[selected] = 1
                if all(b == 1 for b in timesInThisList):
                    for i in range(len(names)):
                        timesInThisList[i]=0
                selectedPeople.append(selected)
                print('Selected %s' % names[selected])
                kitchenCount[selected] += 1
                del lessFikaPeopleIndex[selectedIndex]
            if len(selectedPeople)==2*numberOfWeeks:
                listCompleted = True
                break

    # Create the excel workbook and sheet
    workbook = xlsxwriter.Workbook(kitchenListName)
    worksheet = workbook.add_worksheet()
    
    # Formatting of the excel sheet
    bold = workbook.add_format({'bold': True})
    worksheet.set_column(0, 0, 12)
    worksheet.set_column(1, 1, 21)
    grey = workbook.add_format()
    grey.set_bg_color('#DCDCDC')
    forNames = workbook.add_format()
    forNames.set_bg_color('#DCDCDC')
    forNames.set_border(2)

    firstDayGray = workbook.add_format()
    firstDayGray.set_bg_color('#DCDCDC')
    firstDayGray.set_top(2)
    firstDayGray.set_left(2)
    firstDayGray.set_right(2)

    secondDayGray = workbook.add_format()
    secondDayGray.set_bg_color('#DCDCDC')
    secondDayGray.set_bottom(2)
    secondDayGray.set_left(2)
    secondDayGray.set_right(2)
    secondDayGray.set_align('right')

    firstDayWhite = workbook.add_format()
    firstDayWhite.set_top(2)
    firstDayWhite.set_left(2)
    firstDayWhite.set_right(2)

    secondDayWhite = workbook.add_format()
    secondDayWhite.set_bottom(2)
    secondDayWhite.set_left(2)
    secondDayWhite.set_right(2)
    secondDayWhite.set_align('right')

    # write the header
    worksheet.write(0, 0, 'Date', bold)
    worksheet.write(0, 1, 'Name', bold)

    row = 0

    GrayDate = False

    for k in range(int(len(selectedPeople)/2)):
        # Write both names
        worksheet.write(2*row, 1, names[selectedPeople[2*k]], forNames)
        worksheet.write(2*row+1, 1, names[selectedPeople[2*k+1]], forNames)

        # Write the first date
        if day<10:
            dayStr = '0%d' % (day)
        else:
            dayStr = str(day)
        if GrayDate:
            worksheet.write(2*row, 0, '%s-%s-%d' % (dayStr, months[month - 1], year), firstDayGray)
        else:
            worksheet.write(2*row, 0, '%s-%s-%d' % (dayStr, months[month - 1], year), firstDayWhite)

        # Advance to the second date
        for _ in range(4):
            [day, month, year] = AdvanceDay(day, month, year)

        # Write the first date
        if day < 10:
            dayStr = '0%d' % (day)
        else:
            dayStr = str(day)
        if GrayDate:
            worksheet.write(2 * row+1, 0, '%s-%s-%d' % (dayStr, months[month - 1], year), secondDayGray)
        else:
            worksheet.write(2 * row+1, 0, '%s-%s-%d' % (dayStr, months[month - 1], year), secondDayWhite)

        # Advance to the second date
        for _ in range(3):
            [day, month, year] = AdvanceDay(day, month, year)

        if GrayDate:
            GrayDate = False
        else:
            GrayDate = True

        row += 1

    workbook.close()

    return kitchenCount

def CreateFikaList(fikaListName, startingDate, names, fikaCount, numberOfWeeks):
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    dateMatch = re.search(r'([\d]+)/([\d]+)/([\d]+)', startingDate)
    if dateMatch:
        startingMonth = int(dateMatch.group(1))
        startingDay = int(dateMatch.group(2))
        startingYear = int(dateMatch.group(3))

    month = startingMonth
    day = startingDay
    year = startingYear

    listCompleted = False

    selectedPeople = []

    timesInThisList = []
    for i in range(len(names)):
        timesInThisList.append(0)

    while not listCompleted:
        # Search for the person with max fikaCount
        maxFikaCount = max(fikaCount)
        minFikaCount = min(fikaCount)

        if maxFikaCount == minFikaCount:
            maxFikaCount += 1

        # Get the persons with less fika than the maximum
        lessFikaPeopleIndex = []
        for i in range(len(fikaCount)):
            if fikaCount[i] < maxFikaCount:
                lessFikaPeopleIndex.append(i)

        while len(lessFikaPeopleIndex) > 0:
            selectedIndex = random.randint(0, len(lessFikaPeopleIndex) - 1)
            selected = lessFikaPeopleIndex[selectedIndex]
            if timesInThisList[selected] == 0:
                timesInThisList[selected] = 1
                if all(b == 1 for b in timesInThisList):
                    for i in range(len(names)):
                        timesInThisList[i]=0
                selectedPeople.append(selected)
                print('Selected %s' % names[selected])
                fikaCount[selected] += 1
                del lessFikaPeopleIndex[selectedIndex]
            if len(selectedPeople) == 5 * numberOfWeeks:
                listCompleted = True
                break

    # Create the excel workbook and sheet
    workbook = xlsxwriter.Workbook(fikaListName)
    worksheet = workbook.add_worksheet()

    # Formatting of the excel sheet
    bold = workbook.add_format({'bold': True})
    worksheet.set_column(0, 0, 12)
    worksheet.set_column(1, 1, 21)
    grey = workbook.add_format()
    grey.set_bg_color('#DCDCDC')
    greyAndDoubleLine = workbook.add_format()
    greyAndDoubleLine.set_bg_color('#DCDCDC')
    greyAndDoubleLine.set_top(6)
    greyAndDoubleLine.set_left(2)
    greyAndDoubleLine.set_right(2)
    greyAndDoubleLine.set_bottom(2)
    greyAndSingleLine = workbook.add_format()
    greyAndSingleLine.set_bg_color('#DCDCDC')
    greyAndSingleLine.set_border(2)

    row = 0

    topRow = True

    counter = 0
    for i in selectedPeople:
        if day < 10:
            dayStr = '0%d' % (day)
        else:
            dayStr = str(day)
        if topRow:
            worksheet.write(row, 0, '%s-%s-%d' % (dayStr, months[month - 1], year), greyAndDoubleLine)
            worksheet.write(row, 1, names[i], greyAndDoubleLine)
        else:
            worksheet.write(row, 0, '%s-%s-%d' % (dayStr, months[month - 1], year), greyAndSingleLine)
            worksheet.write(row, 1, names[i], greyAndSingleLine)

        [day, month, year] = AdvanceDay(day, month, year)

        topRow = False
        if counter == 4:
            topRow = True
            counter = -1
            row += 1
            for _ in range(2):
                [day, month, year] = AdvanceDay(day, month, year)
        counter += 1
        row += 1

    workbook.close()

    return fikaCount


def main():

    databaseName = 'FikaListDataBase.xlsx'
    dataBaseSheetName = 'Names & data'
    fikaListName = 'FikaList.xlsx'
    kitchenListName = 'KitchenList.xlsx'
    rowMin = 2
    rowMax = 33 
    startingDate = '10/29/2018'  # mm/dd/yyyy
    createFikaList = True
    createKitchenList = True
    numberOfWeeks = 3


    [names, fikaCount, kitchenCount] = ReadDataBase(databaseName, dataBaseSheetName, rowMin, rowMax)

    if createFikaList:
        print('creating fikalist')
        fikaCount = CreateFikaList(fikaListName, startingDate, names, fikaCount, numberOfWeeks)

    if createKitchenList:
        print('creating kitchenlist')
        kitchenCount = CreateKitchenList(kitchenListName, startingDate, names, kitchenCount, numberOfWeeks)

    if createFikaList or createKitchenList:
        print('creating database')
        WriteDataBase(databaseName, dataBaseSheetName, names, fikaCount, kitchenCount)


    return

# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
