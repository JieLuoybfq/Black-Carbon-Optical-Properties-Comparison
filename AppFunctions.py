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
        Arr = []
        for i in range(len(Array)):
            Arr.append(float(Array[i]))
        n, bins, patches = plt.hist(Arr, 50, density=True, facecolor='b', alpha=0.75)
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


def isNewArrayDB(INFO, TableName, Header, Array):
    try:
        hash_array = DB.getHashofArray(INFO=INFO, TableName=TableName, Header=Header, Array=Array)
        Hashes = DB.checkHashArray(INFO=INFO, TableName=TableName, HashArray=hash_array)
        return Hashes

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


def createRandomNormalArr(Center, Width, Number, digit=3):
    try:
        A = np.random.normal(Center, Width, int(Number))
        length = len(A)
        B = []
        for i in range(length):
            B.append(round(Decimal(A[i]), digit))
        return B

    except Exception as e:
        logging.exception(e)
        raise


def calcMonomerParameter(dpArray, WaveLengthArray):
    try:
        rows = len(dpArray)
        A = []
        for i in range(rows):
            A.append(Decimal(pi) * dpArray[i] / WaveLengthArray[i])

        return A

    except Exception as e:
        logging.exception(e)
        raise


def createRandomNormalArrINT(Center, Width, Number):
    try:
        A = np.random.normal(Center, Width, int(Number))
        length = len(A)
        B = []
        for i in range(length):
            B.append(round(Decimal(A[i]), 0))
        return B

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


def TMatrixOutputDictoArray(Dictionary):
    try:
        A = []
        A.append(Dictionary['TMatrix_ABS_CRS'])
        A.append(Dictionary['TMatrix_SCA_CRS'])
        return A

    except Exception as e:
        logging.exception(e)
        raise


def checkMethodDBforTMatrixIndexes(INFO, TableName, Header, Array):
    try:

        Found = isNewArrayDB(INFO=INFO, TableName=TableName, Header=Header, Array=Array)
        Newinput = []
        OldOutput = []
        OldInput = []
        Indexes = []
        for i in range(len(Found)):
            if Found[i] == -1:
                Newinput.append(Array[i][:])
                Indexes.append(-1)
            else:
                OldOutput.append(TMatrixOutputDictoArray(DB.getOutputRowByHash(INFO=INFO, TableName=TableName, Hash=Found[i])))
                OldInput.append(Array[i][:])
                Indexes.append(1)

        return OldInput, OldOutput, Newinput, Indexes

    except Exception as e:
        logging.exception(e)
        raise


def RDGOutputDictoArray(Dictionary):
    try:
        A = []
        A.append(Dictionary['RDG_ABS_CRS'])
        A.append(Dictionary['RDG_SCA_CRS'])
        return A

    except Exception as e:
        logging.exception(e)
        raise


def checkMethodDBforRDGIndexes(INFO, TableName, Header, Array):
    try:

        Found = isNewArrayDB(INFO=INFO, TableName=TableName, Header=Header, Array=Array)
        Newinput = []
        OldOutput = []
        OldInput = []
        Indexes = []
        for i in range(len(Found)):
            if Found[i] == -1:
                Newinput.append(Array[i][:])
                Indexes.append(-1)
            else:
                OldOutput.append(RDGOutputDictoArray(DB.getOutputRowByHash(INFO=INFO, TableName=TableName, Hash=Found[i])))
                OldInput.append(Array[i][:])
                Indexes.append(1)

        return OldInput, OldOutput, Newinput, Indexes

    except Exception as e:
        logging.exception(e)
        raise


def joinArray(*args):
    try:
        rows = len(args[0])
        Arr = []
        for i in range(rows):
            A = []
            for j in range(len(args)):
                A.append(args[j][i])
            Arr.append(A)
        return Arr

    except Exception as e:
        logging.exception(e)
        raise


def createConstantArray(Number, Howmany):
    try:
        B = []
        for i in range(Howmany):
            B.append(round(Decimal(Number), 1))
        return B

    except Exception as e:
        logging.exception(e)
        raise


def createIndexArray(start, len):
    try:
        B = []
        for i in range(start, len):
            B.append(int(i))
        return B

    except Exception as e:
        logging.exception(e)
        raise


def getPossibleArray(Array, Indexes):
    try:
        B = []
        for i in Indexes:
            B.append(Array[i])
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


def addMlinkToArray(Array, LastID):
    try:
        arrLen = len(Array)
        arrCol = len(Array[0])
        MainArr = []
        Mlink = list(range(LastID, arrLen + LastID))
        for i in range(arrLen):
            A = []
            for j in range(arrCol):
                A.append(Array[i][j])
            A.append(Mlink[i])
            MainArr.append(A)
        return MainArr

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
    except Exception as e:
        logging.exception(e)
        raise
