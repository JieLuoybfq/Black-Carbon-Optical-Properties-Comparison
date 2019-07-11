import DBMethods as DB
from ConfigParserM import logging
import GeneralFunctions as GF
from decimal import Decimal
import numpy as np
from matplotlib import rcParams
from matplotlib import pyplot as plt
from math import pi

####### Plotting Parameters

rcParams['mathtext.fontset'] = 'stix'
rcParams['font.family'] = 'STIXGeneral'


##############################


def toSaveHistogram(Folder, Name, Array, Figure_DPI=400):
    try:
        n, bins, patches = plt.hist(Array, 50, density=True, facecolor='b', alpha=0.75)
        plt.ylabel('Probability')
        plt.title('Histogram of ' + str(Name))
        plt.grid(True)
        Address = GF.getAddressTo(Folder, None, Name, "jpg")
        plt.savefig(Address, format='jpg', dpi=Figure_DPI, bbox_inches='tight')
        plt.clf()
        plt.close()

    except Exception as e:
        logging.exception(e)
        raise


def isNewInputDB(INFO, TableName, InputDictionary):
    try:

        # DB.createTable(INFO=INFO, TableName=TableName, arrHeaderNamesInput=InputHeader, arrHeaderNamesOutput=OutputHeader)
        hash_Database = DB.getHashofRow(INFO=INFO, TableName=TableName, ValuesDictionary=InputDictionary)
        ID_Database, rowCount_Database = DB.checkHashValue(INFO=INFO, TableName=TableName, Hash=hash_Database)

        return ID_Database, rowCount_Database

    except Exception as e:
        logging.exception(e)
        raise


def uniqueList(list1):
    try:
        unique_list = []
        for x in list1:
            if x not in unique_list:
                unique_list.append(x)
        return unique_list

    except Exception as e:
        logging.exception(e)
        raise


def uniqueEntry(Array):
    try:
        col = len(Array[0])
        ranges = []
        for k in range(col):
            ranges.append(uniqueList(GF.selectColumnsList(ColumnIndex=[k], List=Array)))

        return ranges
    except Exception as e:
        logging.exception(e)
        raise


def findNearest(Array, input):
    try:
        inp = Decimal(input)
        A = min(Array, key=lambda x: abs(x - inp))
        return A

    except Exception as e:
        logging.exception(e)
        raise


def findIndexNearest(Array, input):
    try:
        inp = Decimal(input)
        A = min(range(len(Array)), key=lambda i: abs(Array[i] - inp))
        return A

    except Exception as e:
        logging.exception(e)
        raise


def checkIndex(tolerance, MainIndex, max):
    try:
        up = MainIndex + tolerance
        down = MainIndex - tolerance
        if down < 0:
            down = 0
        if up > max:
            up = max
        return [up, down]


    except Exception as e:
        logging.exception(e)
        raise


def createRandomNormalArr(Center, Width, Number):
    try:
        A = np.random.normal(Center, Width, int(Number))
        return np.around(A, 3)

    except Exception as e:
        logging.exception(e)
        raise


def calcMonomerParameter(dpArray, WaveLengthArray):
    try:
        rows = len(dpArray)
        A = []
        for i in range(rows):
            A.append(pi * dpArray[i] / WaveLengthArray[i])

        return A

    except Exception as e:
        logging.exception(e)
        raise


def createRandomNormalArrINT(Center, Width, Number):
    try:
        A = np.random.normal(Center, Width, int(Number))
        return np.around(A, 0)

    except Exception as e:
        logging.exception(e)
        raise


def findCommonIndex(*args):
    try:
        A = len(args)
        C = args[0]
        for i in range(1, A):
            C = set(args[i]) & set(C)
        return list(C)

    except Exception as e:
        logging.exception(e)
        raise


def checkRDGDBforIndexes(*args):
    try:
        A = len(args)
        C = args[0]
        for i in range(1, A):
            C = set(args[i]) & set(C)
        return list(C)

    except Exception as e:
        logging.exception(e)
        raise


def checkTMatrixDBforIndexes(INFO, TableName, **kwargs):
    try:
        rows = len(kwargs['Df'])
        A = kwargs['Df'][1]
        # for i in range(rows):
        # dict1={
        #    'Df':kwargs['Df'][i], 'kf':, 'R_RI':, 'I_RI':, 'WaveL':, 'dp':, 'Np':, 'Version':
        # }
        A = 3




    except Exception as e:
        logging.exception(e)
        raise


def getPossibleArray(Array, Indexes):
    try:
        B = []
        for i in Indexes:
            B.append(Decimal(Array[i]))
        return B

    except Exception as e:
        logging.exception(e)
        raise


def getRandomFromArr(Array, Number):
    try:
        # uniform choice
        A = np.random.choice(Array, size=int(Number), replace=False)
        return A

    except Exception as e:
        logging.exception(e)
        raise


def getGoodIndexes(Array, Bound):
    try:
        min = Bound[0]
        max = Bound[1]
        Index = []
        rows = len(Array)
        for i in range(rows):
            if (Array[i] >= min) and (Array[i] <= max):
                Index.append(i)
        return Index

    except Exception as e:
        logging.exception(e)
        raise


def getToleratedArray(Array, Input, Tolerance, uniques):
    try:
        B = []
        index = []
        rows = len(Array)
        col = len(Array[0])
        nearest = []
        nearest_Index = []
        accepted_Index = []

        for i in range(col):
            nearest.append(findNearest(uniques[i], Input[i]))
            nearest_Index.append(findIndexNearest(uniques[i], Input[i]))

        for i in range(col):
            accepted_Index.append(checkIndex(Tolerance[i], nearest_Index[i], len(uniques[i]) - 1))

        for i in range(rows):
            k = 0
            for j in range(col):
                if (Array[i][j] >= uniques[j][accepted_Index[j][1]]) and (Array[i][j] <= uniques[j][accepted_Index[j][0]]):
                    k += 1
            if k == col:
                B.append(Array[i][:])
                index.append(i)

        return B, index
        '''
        for i in range(rows):
            k = 0
            for j in range(col):
                if (Array[i][j] >= Input[j] - Tolerance[j]) and (Array[i][j] <= Input[j] + Tolerance[j]):
                    k += 1
            if k == col:
                B.append(Array[i][:])
                index.append(i)

        return B, index
        '''
    except Exception as e:
        logging.exception(e)
        raise
