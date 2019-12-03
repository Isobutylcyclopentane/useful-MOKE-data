#MOKE_data analysis
import math
import csv
from scipy import optimize
import statistics
import warnings

def read_MOKE_data(fileName):

    result = dict()
    result[1] = []
    result[2] = []
    result[3] = []
    result[4] = []

    with open(fileName, 'r',encoding="UTF-8-sig") as csvfile:
        reader = csv.reader(csvfile)
        for line in reader:
            x, y = line
            x, y = float(x), float(y)
            if x>=0 and y>=0:
                result[1].append((x,y))
            elif x<=0 and y>=0:
                result[2].append((x,y))
            elif x<=0 and y<=0:
                result[3].append((x,y))
            elif x>=0 and y<=0:
                result[4].append((x,y))

    for index in result:
        result[index].sort()

    return result

def find_Mr(data):
    # get the data for residual magnetization and uncertainty
    # 1in: dict; 2out: float, float

    upperMr = find_intersection(data[2],data[1],0)
    lowerMr = find_intersection(data[2],data[1],0)

    # minimizing experimental error by centering the data
    # and find uncertainty by the uncentered value
    uncentered = abs(upperMr) - abs(lowerMr)
    Mr = abs(upperMr - uncentered)

    uncertainty = abs(uncentered) / 2

    return Mr,uncertainty

def find_raw_Hc(data):
    #dataUpper = [p for p in data and p[1] >= 0]
    #dataLower = [p for p in data and p[1] < 0]
    pointAnalysis = get_points(data, 1)
    if len(pointAnalysis) == 1:
        print("***ERROR: incomplete data set,",\
              "the current experimental data aborted")
        return False, True
    elif len(pointAnalysis) > 2:
        print("***WARNING: multiple loop in one data set,",\
              "currently program cannot analyze this type of data")
        warnings.warn("***WARNING: multiple loops in one data set")
        return False, False


    #leftHc = find_intersection(data[2],data[3], 1)
    #rightHc = find_intersection(data[1],data[4], 1)
    leftHc = find_point_intersection(pointAnalysis[0], 0)
    rightHc = find_point_intersection(pointAnalysis[1], 0)

    return leftHc, rightHc

def find_raw_Mr (data):
    pointAnalysis = get_points(data, 0)
    if len(pointAnalysis) == 1:
        print("***ERROR: incomplete data set,", \
              "the currect experimental data aborted")
        return False, True
    elif len(pointAnalysis) > 2:
        print("***WARNING: multiple loop in one data set", \
              "currently program cannot analyze this type of data")

        return False, False

    upperMr = find_point_intersection(pointAnalysis[0], 1)
    lowerMr = find_point_intersection(pointAnalysis[1], 1)

    return upperMr, lowerMr

def get_points(data, mode):
    result = []
    numOfData = math.ceil(len(data) / 20)
    numGet = math.ceil(numOfData / 2)
    for index in range(len(data)):
        if index == 0:
            continue
        currentSign = int(abs(data[index][mode]) / data[index][mode])
        previousSign = int(abs(data[index-1][mode]) / data[index-1][mode])
        if currentSign != previousSign:
            temp = data[index-numGet:index+numGet]
            result.append(temp)
    return result

def find_point_intersection(points, mode):
    dataX = [x[0] for x in points]
    dataY = [x[1] for x in points]
    popt, pcov = list(optimize.curve_fit(linear, dataX, dataY))
    a, b = popt
    if mode == 0:
        return (-b/a)
    else:
        return b

#old find_intersection
'''
def find_intersection(quadrum1, quadrum2, mode):
    quadrum1.sort(key= lambda x: abs(x[mode]))
    quadrum2.sort(key= lambda x: abs(x[mode]))

    numPtsAna = math.ceil(min(len(quadrum1),len(quadrum2)) / 5)

    point_data = quadrum1[0:numPtsAna] + quadrum2[0:numPtsAna]
    point_data.sort()

    dataX = list(map(lambda x: x[0],point_data))
    dataY = list(map(lambda x: x[1],point_data))

    popt,pcov = list(optimize.curve_fit(linear, dataX, dataY))
    a,b = popt
    if mode == 1:
        return (b/a)
    elif mode == 0:
        return b
'''

def find_Ms (data):

    Quad1Result = find_saturation(data[1])
    Quad3Result = find_saturation(data[3])

    uncent = Quad1Result[0] + Quad3Result[0]

    Ms = Quad1Result[0] - uncent
    uncertantity = Quad1Result[1] + abs(uncent)
    result = Ms, uncertantity

    return result


def find_saturation(quadrum):

    temp_y = []
    temp = []
    for i in range(6):
        tmp_val = quadrum[-(i+1)]
        temp_y.append(tmp_val[1])
        temp.append(tmp_val)

    for index in range((len(quadrum)-7),-1,-1):
        tmp_val = quadrum[index]

        stdv = statistics.stdev(temp_y)
        avg = statistics.mean(temp_y)
        upper_bound = avg + stdv
        lower_bound = avg - stdv

        if not (lower_bound < tmp_val[1] < upper_bound):
            break
        temp.append(tmp_val)
        temp_y.append(tmp_val[1])
    avg = statistics.mean(temp_y)
    stdv = statistics.stdev(temp_y)
    result = (avg,stdv)
    return result


def linear(x,a,b):
    return a * x + b


def MOKE_analysis(data):
    # complete analysis with a set of MOKE data
    # include find coercivity(Hc), residual magnetization(Mr), and saturation(Ms)
    # 1in: dict; 1out: dict

    result = dict()

    #result["Hc_result"] = find_Hc(data)
    result["Mr_result"] = find_Mr(data)
    result["Ms_result"] = find_Ms(data)

    return result


def moke_main(fileName = None):
    if fileName == None:
        fileName = input("please enter the name of the data file (.csv required) => ")
    fileName += ".csv"
    data = read_MOKE_data(fileName)

    result = MOKE_analysis(data)
    Hc_result = result["Hc_result"]
    Mr_result = result["Mr_result"]
    Ms_result = result["Ms_result"]

    print("The coercivity strength (Hc) of the material is %.5f , with an uncertainty of %.5f" % Hc_result)
    print("The residual magnetization (Mr) of the material is %.5f, with an uncertainty of %.5f" % Mr_result)
    print("The magnetic saturation (Ms) of the material is %.5f, with an uncertainty of %.5f" % Ms_result)

    command = input("\nmore accurate data? (y/n) => ")
    if command.strip().lower() == 'y':
        print()
        print("Hc:",Hc_result[0],'±',Hc_result[1])
        print("Mr:",Mr_result[0],'±',Mr_result[1])
        print("Ms:",Ms_result[0],'±',Ms_result[1])
    elif command.strip().lower() =='n':
        return
    else:
        print("Invalid command")

if __name__ == "__main__":
    #fileName = input("please type in the name of the data file => ")
    fileName = "MOKE_CuCo"
    moke_main(fileName)


