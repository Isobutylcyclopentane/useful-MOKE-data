#MOKE Azimuthal angel analysis

from MOKE_data import *
import math
import statistics
import os
import glob
import ntpath
from matplotlib import pyplot as plt
import numpy as np

class Experiment_Result(object):
    def __init__(self):
        self.storage = dict()

    def append(self, other):
        if other.aziAngle > 360:
            other.aziAngle -= 360
        if other.aziAngle in self.storage:
            self.storage[other.aziAngle].update(other)
        else:
            self.storage[other.aziAngle] = other
        return

    def get_all_azimuthal_Hc (self):
        dataX = list()
        dataY = list()
        all_angle = sorted(self.storage)
        for azi in all_angle:
            x = (azi / 180) * math.pi
            dataX.append(x)
            y = self.storage[azi].Hc
            dataY.append(y)
        return dataX, dataY

    def get_all_skewness (self):
        dataX = list()
        dataY = list()
        all_angle = sorted(self.storage)
        for azi in all_angle:
            x = (azi / 180) * math.pi
            dataX.append(x)
            y = self.storage[azi].skewness
            dataY.append(y)
        dataX = np.array(dataX)
        dataY = np.array(dataY)
        return dataX, dataY

    def get_amount(self):
        return len(self.storage)

    def check_completeness(self):
        flag = True
        maxAngle = max(self.storage)
        for a in range(maxAngle, -1, -15):
            if a not in self.storage:
                flag = False
                break
        return flag

    def rough_fill(self):
        maxAngle = max(self.storage)
        for a in range(maxAngle, -1, -15):
            if a not in self.storage:
                self.storage[a] = MOKE_result(a, (0,0), (0,0), (0,0))
        return

    def __str__(self):
        result = ""
        order = sorted(self.storage)
        for key in order:
            result += str(self.storage[key])
            result += '\n'
        return result

class MOKE_result(object):
    def __init__(self, aziAngle, Hc, Mr, Ms):
        self.aziAngle = aziAngle
        self.Hc = Hc[0]
        self.HcErr = Hc[1]
        self.Mr = Mr[0]
        self.MrErr = Mr[1]
        self.Ms = Ms[0]
        self.MsErr = Ms[1]
        self.skewness, self.skewErr = self.get_skewness()

    def __str__(self):
        result = "{:0>3d}".format(self.aziAngle) + \
                 '|' +  str(list((self.Hc, self.Mr, self.Ms)))
        return result

    def get_skewness_and_error(self):

        if self.Ms == 0:
            return (0,0)

        skewness = self.Mr / self.Ms

        skewnessMax = (self.Mr + self.MrErr) / (self.Ms - self.MsErr)
        skewnessMin = (self.Mr - self.MrErr) / (self.Ms + self.MrErr)
        error = (skewnessMax - skewnessMin) / 2

        return skewness, error

    def get_skewness(self):
        if self.Ms == 0:
            return 0, 0
        skewness = self.Mr / self.Ms
        return skewness, 0


    def update(self, other):

        if self.aziAngle != other.aziAngle:
            print("***Warning, updating data measured \
                    at different azimuthal angle")
            print("self.azimuthal", self.aziAngle, \
                  "other.azimuthal", other.aziAngle)

        self.Hc = (self.Hc + other.Hc) / 2
        self.HcErr = (self.HcErr + other.HcErr) / 2
        self.Mr = (self.Mr + other.Mr) / 2
        self.MrErr = (self.MrErr + other.MrErr) / 2
        self.Ms = (self.Ms + other.Ms) / 2
        self.MsErr = (self.MsErr + other.MsErr) / 2
        self.skewness, self.skewErr = self.get_skewness()
        return

#***************************************************************************
#helper functions
#***************************************************************************

def read_file(filename):
    resultData = []
    file = open(filename, encoding="UTF-8-sig")

    for line in file:

        if '\t' in line:
            tempLst = line.split('\t')
            tempVal = float(tempLst[1]), float(tempLst[2])
            resultData.append(tempVal)

    return resultData


def rough_centering(dataLst):
    result = []

    numOfAna = math.ceil(len(dataLst)/10)

    dataSortX = sorted(dataLst)
    dataSortY = sorted(dataLst, key=lambda x: x[1])

    dataXSort = [x[0] for x in dataSortX]
    dataYSort = [y[1] for y in dataSortY]

    allXMedium = statistics.median(dataXSort)
    allYMedium = statistics.median(dataYSort)

    xMinLst = dataSortX[:numOfAna]
    xMin = statistics.mean([x[0] for x in xMinLst])

    xMaxLst = dataSortX[-numOfAna:]
    xMax = statistics.mean([x[0] for x in xMaxLst])

    yMinLst = dataSortY[:numOfAna]
    yMin = statistics.mean([y[1] for y in yMinLst])

    yMaxLst = dataSortY[-numOfAna:]
    yMax = statistics.mean([y[1] for y in yMaxLst])

    adjX = (xMax + xMin) / 2
    adjY = (yMax + yMin) / 2

    for point in dataLst:
        newPoint = (point[0] - allXMedium),(point[1] - adjY)
        result.append(newPoint)
    return result

def separate_quardum(dataLst):
    result = dict()
    result[1] = []
    result[2] = []
    result[3] = []
    result[4] = []

    for point in dataLst:
        x,y = point
        if x >= 0 and y >= 0:
            result[1].append((x, y))
        elif x < 0 and y >= 0:
            result[2].append((x, y))
        elif x < 0 and y < 0:
            result[3].append((x, y))
        elif x >=0 and y < 0:
            result[4].append((x,y))

    return result


def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

def get_azimuthal_angle(filename):
    temp = filename.split('_')[0].split(' ')[-1]
    result = int(temp)
    return result

def find_proper_Hc(leftHc, rightHc):
    uncent = (rightHc + leftHc) / 2
    Hc = rightHc - uncent
    return Hc

def find_proper_Mr (upperMr, lowerMr):
    uncent = (upperMr + lowerMr) / 2
    Mr = upperMr - uncent
    return Mr


def adjust_data(data, adjX=0, adjY=0):
    result = []
    for point in data:
        x,y = point
        x -= adjX
        y -= adjY
        result.append((x, y))
    return result


def plot_and_save (filename, location, expData):
    path = "C:\\Users\\J. Cheng\\Dropbox\\AFM_image\\1024_ADV_MOKE\\10052018"
    newFileName = filename.replace(".xlsx", ".jpg")
    newPathName = path + '\\'+location+'\\' + "CENTERED " + newFileName
    # newFile = open(newPathName, "w+")
    dataX = np.array([x[0] for x in expData])
    dataY = np.array([x[1] for x in expData])
    plt.plot(dataX, dataY)
    plt.grid(True)
    plt.savefig(newPathName)
    plt.close()
    return

def delete_all_files ():
    pathBase = "C:\\Users\\J. Cheng\\Dropbox\\AFM_image\\1024_ADV_MOKE\\10052018"
    fileDelete = open("DELETE.txt")
    fileStr = fileDelete.read()
    fileStr = fileStr.split('\n')
    for filename in fileStr:
        filename = filename.strip()
        pathname = os.path.join(pathBase, filename)
        filelist = [f for f in os.listdir(pathname)]
        for f in filelist:
            temp = os.path.join(pathname, f)
            os.remove(temp)
    fileDelete.close()
    return



if __name__ == '__main__':

    # delete_all_files()

    result = Experiment_Result()
    path = "C:\\Users\\J. Cheng\\Dropbox\\AFM_image\\1111_MULTILOOP\\DATA"

    goodDataCount = 0
    badDataCount = 0

    for pathname in glob.glob(os.path.join(path, '*.xlsx')):
        filename = path_leaf(pathname)
        azimuthal = get_azimuthal_angle(filename)


        print("\nDEBUG filename", filename)

        expData = read_file(pathname)
        expData = rough_centering(expData)
        #expData = separate_quardum(expData)
        leftHc, rightHc = find_raw_Hc(expData)

        if ((leftHc == False) and (rightHc == True)):
            plot_and_save(filename, "INCOMPLETE", expData)
            badDataCount += 1
            continue

        if ((leftHc == False) and (rightHc == False)):
            plot_and_save(filename, "MULTILOOP", expData)
            badDataCount += 1
            continue

        uncentAdjX = (rightHc + leftHc) / 2
        expData = adjust_data(expData, adjX=uncentAdjX)

        Hc = find_proper_Hc(leftHc, rightHc), 0

        upperMr, lowerMr = find_raw_Mr(expData)
        Mr = find_proper_Mr(upperMr, lowerMr), 0

        plot_and_save(filename, "ANALYZED", expData)

        expDataDict = separate_quardum(expData)

        '''
        flag = False
        for quardum in expDataDict:
            if len(expDataDict[quardum]) == 0:
                plot_and_save(filename, "MISSPOINT", expData)
                flag = True
        if flag:
            badDataCount += 1
            continue
        '''

        '''
        try:

            Ms = find_Ms(expDataDict)
            fileResult = MOKE_result(azimuthal, Hc, Mr, Ms)
            result.append(fileResult)
            goodDataCount += 1
        except:
                #for debugging the auto-centering
                badDataCount += 1
                plot_and_save(filename, "CENTERED", expData)

                #newFile.close()
        '''
        Ms = find_Ms(expDataDict)
        fileResult = MOKE_result(azimuthal, Hc, Mr, Ms)
        result.append(fileResult)
        goodDataCount += 1

    print("Data accepted and analyzed:", goodDataCount)
    print("Data rejected:", badDataCount)
    print("Data storage amount:", result.get_amount())
    print("Data completeness:", result.check_completeness())
    print("all analyzed data before: \n"+str(result))
    result.rough_fill()
    print("all analyzed data after: \n"+str(result))

    skewX, skewY = result.get_all_skewness()
    plt.polar(skewX, skewY)
    plt.title("squareness")
    plt.show(block = True)
    HcX, HcY = result.get_all_azimuthal_Hc()
    plt.polar(HcX, HcY)
    plt.title("Hc")
    plt.show(block=True)



    '''
        #for debugging the auto-centering
        
        newFileName = filename.replace(".xlsx", ".txt")
        newPathName = path + '\\CENTERED\\' + "CENTERED " + filename
        newFile = open(newPathName, "w+")
        for point in expData:
            line =  '\t' + str(point[0]) + '\t' + str(point[1]) + "\t\n"
            newFile.write(line)
        newFile.close()
    '''

    '''
        try:
            expParameters = MOKE_analysis(expData)
            expResult = MOKE_result(azimuthal, expParameters["Hc_result"], \
                                expParameters["Mr_result"], expParameters["Ms_result"])
            goodDataCount += 1
        except:
            badDataCount += 1
    '''



