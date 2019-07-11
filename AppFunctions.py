import DBMethods as DB
from ConfigParserM import logging
import GeneralFunctions as GF
from decimal import Decimal


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
        b = 3

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
