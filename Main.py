# Programmed by Keyhan Babaee Under Prof. Steven Rogak supervision, https://github.com/KeyhanB
# Using Fortran Code from: C. Liu, X. Xu, Y. Yin, M. Schnaiter, Y. L. Yung, 2018: Black carbon aggregate: A database for optical properties
# Version 0.1
# July 2019
import ConfigParserM as CP
import DBMethods as DB
from ConfigParserM import logging
import GeneralFunctions as GF
import AppFunctions as FN
from scipy.interpolate import griddata

if __name__ == "__main__":
    CP.readLogConfig()
    DB_Info = CP.readConfigToDict(SectionName="DatabaseInfo")
    FF_Info = CP.readConfigToDict(SectionName="FilesFoldersInfo")
    AGG_Info = CP.readConfigToDict(SectionName="AggregateDetails", ConvertParseTo='float', hasComment=True)

    # DB.createDB(INFO=DB_Info)
    logging.info("Application Started!")
    tableName = 'Raw_V1'
    columnName, dataFull = DB.readAllRowsfromTable(INFO=DB_Info, TableName=tableName)
    columnName = GF.selectColumnsList(ColumnIndex=[1, 2, 3, 4, 5], List=columnName, Dimension=1)
    EXT_Coeff_Full = GF.selectColumnsList(ColumnIndex=[6], List=dataFull)
    SCT_Coeff_Full = GF.selectColumnsList(ColumnIndex=[7], List=dataFull)
    ABS_Coeff_Full = GF.selectColumnsList(ColumnIndex=[8], List=dataFull)
    inputData_Full = GF.selectColumnsList(ColumnIndex=[1, 2, 3, 4, 5], List=dataFull)

    uniqueColumn = FN.uniqueEntry(inputData_Full)
    In = [2.3, 1.4, 0.6, 0.2, 52]
    Tolerance = [1, 2, 2, 2, 2]

    X, ind = FN.getToleratedArray(inputData_Full, In, Tolerance, uniqueColumn)
    out = GF.selectColumnsList(ind, EXT_Coeff_Full, Dimension=1)
    A = griddata(X, out, In)

    In = [2.5, 1.3, 0.7, 0.2, 2100]
    # Tolerance = [2, 2, 2, 2, 7]
    X, ind = FN.getToleratedArray(inputData_Full, In, Tolerance, uniqueColumn)
    out = GF.selectColumnsList(ind, EXT_Coeff_Full, Dimension=1)
    A1 = griddata(X, out, In, rescale=True)

    In = [2.5, 1.3, 0.7, 0.2, 2100]
    # Tolerance = [2, 2, 2, 2, 7]
    X, ind = FN.getToleratedArray(inputData_Full, In, Tolerance, uniqueColumn)
    out = GF.selectColumnsList(ind, ABS_Coeff_Full, Dimension=1)
    A2 = griddata(X, out, In, rescale=True)
    b = 3
