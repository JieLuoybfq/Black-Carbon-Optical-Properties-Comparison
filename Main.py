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
    #######################
    DB.showAllTablesInDBSummary(DB_Info)
    appDirectory = GF.getRootDirectory()
    #######################
    saveHistogram = True
    ################################################################################################################
    RDGTableName = "RDG_V1"
    RDGInputHeaders = ['Df', 'kf', 'R_RI', 'I_RI', 'WaveL', 'dp', 'Np', 'sigma', 'Version']
    RDGOutputHeaders = ['RDG_ABS_CRS', 'RDG_SCA_CRS']
    DB.createTable(INFO=DB_Info, TableName=RDGTableName, arrHeaderNamesInput=RDGInputHeaders, arrHeaderNamesOutput=RDGOutputHeaders)
    ################################################################################################################
    TMatrixTableName = "TMatrix_V1"
    TMatrixInputHeaders = ['Df', 'kf', 'R_RI', 'I_RI', 'WaveL', 'dp', 'Np', 'Version']
    TMatrixOutputHeaders = ['RDG_ABS_CRS', 'RDG_SCA_CRS']
    DB.createTable(INFO=DB_Info, TableName=TMatrixTableName, arrHeaderNamesInput=TMatrixInputHeaders, arrHeaderNamesOutput=TMatrixOutputHeaders)
    ################################################################################################################
    ErrorTableName = "Error_V1"
    ErrorInputHeaders = ['Df', 'kf', 'R_RI', 'I_RI', 'WaveL', 'dp', 'Np', 'sigma', 'Version']
    ErrorOutputHeaders = ['RDG_ABS_CRS', 'RDG_SCA_CRS']
    DB.createTable(INFO=DB_Info, TableName=ErrorTableName, arrHeaderNamesInput=ErrorInputHeaders, arrHeaderNamesOutput=ErrorOutputHeaders)
    ################################################################################################################
    # DB.createDB(INFO=DB_Info)
    DB.showAllTablesInDBSummary(DB_Info)
    # DB.dropTableSet(DB_Info,RDGTableName)
    # DB.dropTableSet(DB_Info, TMatrixTableName)
    # DB.dropTableSet(DB_Info, ErrorTableName)
    # DB.reinitializeDB(DB_Info)
    ################################################################################################################
    RDGInputTableHeaderDB, RDGInputDataFullDB = DB.readAllRowsfromTable(INFO=DB_Info, TableName=RDGTableName)
    TMatrixInputTableHeaderDB, TMatrixInputDataFullDB = DB.readAllRowsfromTable(INFO=DB_Info, TableName=TMatrixTableName)
    ErrorInputTableHeaderDB, ErrorInputDataFullDB = DB.readAllRowsfromTable(INFO=DB_Info, TableName=ErrorTableName)
    ################################################################################################################
    RDGOutputTableHeaderDB, RDGOutputDataFullDB = DB.readAllRowsfromTable(INFO=DB_Info, TableName=RDGTableName + "_Out")
    TMatrixOutputTableHeaderDB, TMatrixOutputDataFullDB = DB.readAllRowsfromTable(INFO=DB_Info, TableName=TMatrixTableName + "_Out")
    ErrorOutputTableHeaderDB, ErrorOutputDataFullDB = DB.readAllRowsfromTable(INFO=DB_Info, TableName=ErrorTableName + "_Out")
    ################################################################################################################
    arrAgg_Fractal_Dimension = FN.createRandomNormalArr(Center=AGG_Info['AGG_FRACTAL_DIMENSION_CENTER'], Width=AGG_Info['AGG_FRACTAL_DIMENSION_STANDARD_DEVIATION'], Number=AGG_Info['MONTECARLO_ARRAY_SIZE'])
    arrAgg_Fractal_Prefactor = FN.createRandomNormalArr(Center=AGG_Info['AGG_FRACTAL_PREFACTOR_CENTER'], Width=AGG_Info['AGG_FRACTAL_PREFACTOR_STANDARD_DEVIATION'], Number=AGG_Info['MONTECARLO_ARRAY_SIZE'])
    arrAgg_RI_Real = FN.createRandomNormalArr(Center=AGG_Info['AGG_RI_REAL_CENTER'], Width=AGG_Info['AGG_RI_REAL_STANDARD_DEVIATION'], Number=AGG_Info['MONTECARLO_ARRAY_SIZE'])
    arrAgg_RI_Imag = FN.createRandomNormalArr(Center=AGG_Info['AGG_RI_IMAG_CENTER'], Width=AGG_Info['AGG_RI_IMAG_STANDARD_DEVIATION'], Number=AGG_Info['MONTECARLO_ARRAY_SIZE'])
    arrAgg_WLength = FN.createRandomNormalArr(Center=AGG_Info['AGG_WLENGTH_CENTER'], Width=AGG_Info['AGG_WLENGTH_STANDARD_DEVIATION'], Number=AGG_Info['MONTECARLO_ARRAY_SIZE'])
    arrAgg_Primary_Diameter = FN.createRandomNormalArr(Center=AGG_Info['AGG_PRIMARY_DIAMETER_CENTER'], Width=AGG_Info['AGG_PRIMARY_DIAMETER_STANDARD_DEVIATION'], Number=AGG_Info['MONTECARLO_ARRAY_SIZE'])
    arrAgg_Polydispersity_Sigma_Within = FN.createRandomNormalArr(Center=AGG_Info['AGG_POLYDISPERSITY_SIGMA_WITHIN_CENTER'], Width=AGG_Info['AGG_POLYDISPERSITY_SIGMA_WITHIN_STANDARD_DEVIATION'], Number=AGG_Info['MONTECARLO_ARRAY_SIZE'])
    arrAgg_Polydispersity_Sigma_Each_Mobility = FN.createRandomNormalArr(Center=AGG_Info['AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER'], Width=AGG_Info['AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_STANDARD_DEVIATION'], Number=AGG_Info['MONTECARLO_ARRAY_SIZE'])
    arrAgg_Prefactor_Projected_Area = FN.createRandomNormalArr(Center=AGG_Info['AGG_PREFACTOR_PROJECTED_AREA_COEFFICIENT_CENTER'], Width=AGG_Info['AGG_PREFACTOR_PROJECTED_AREA_COEFFICIENT_STANDARD_DEVIATION'], Number=AGG_Info['MONTECARLO_ARRAY_SIZE'])
    arrAgg_Exponent_Projected_Area = FN.createRandomNormalArr(Center=AGG_Info['AGG_EXPONENT_PROJECTED_AREA_COEFFICIENT_CENTER'], Width=AGG_Info['AGG_EXPONENT_PROJECTED_AREA_COEFFICIENT_STANDARD_DEVIATION'], Number=AGG_Info['MONTECARLO_ARRAY_SIZE'])
    arrAgg_Primary_Number = FN.createRandomNormalArrINT(Center=AGG_Info['AGG_PRIMARY_NUMBER_CENTER'], Width=AGG_Info['AGG_PRIMARY_NUMBER_STANDARD_DEVIATION'], Number=AGG_Info['MONTECARLO_ARRAY_SIZE'])

    if saveHistogram:
        histogram_Folder = GF.getAddressTo(appDirectory, FF_Info['FOLDER_NAME_GRAPH'] + "/Normal Histograms")
        FN.toSaveHistogram(Folder=histogram_Folder, Name="Agg_Fractal_Dimension", Array=arrAgg_Fractal_Dimension)
        FN.toSaveHistogram(Folder=histogram_Folder, Name="Agg_Fractal_Prefactor", Array=arrAgg_Fractal_Prefactor)
        FN.toSaveHistogram(Folder=histogram_Folder, Name="Agg_RI_Real", Array=arrAgg_RI_Real)
        FN.toSaveHistogram(Folder=histogram_Folder, Name="Agg_RI_Imag", Array=arrAgg_RI_Imag)
        FN.toSaveHistogram(Folder=histogram_Folder, Name="Agg_WLength", Array=arrAgg_WLength)
        FN.toSaveHistogram(Folder=histogram_Folder, Name="Agg_Primary_Diameter", Array=arrAgg_Primary_Diameter)
        FN.toSaveHistogram(Folder=histogram_Folder, Name="Agg_Polydispersity_Sigma_Within", Array=arrAgg_Polydispersity_Sigma_Within)
        FN.toSaveHistogram(Folder=histogram_Folder, Name="Agg_Polydispersity_Sigma_Each_Mobility", Array=arrAgg_Polydispersity_Sigma_Each_Mobility)
        FN.toSaveHistogram(Folder=histogram_Folder, Name="Agg_Prefactor_Projected_Area", Array=arrAgg_Prefactor_Projected_Area)
        FN.toSaveHistogram(Folder=histogram_Folder, Name="Agg_Exponent_Projected_Area", Array=arrAgg_Exponent_Projected_Area)
        FN.toSaveHistogram(Folder=histogram_Folder, Name="Agg_Primary_Number", Array=arrAgg_Primary_Number)

    arrAgg_Fractal_Dimension_Random = FN.getRandomFromArr(Array=arrAgg_Fractal_Dimension, Number=AGG_Info['MONTECARLO_RANDOM_SIZE'])
    arrAgg_Fractal_Prefactor_Random = FN.getRandomFromArr(Array=arrAgg_Fractal_Prefactor, Number=AGG_Info['MONTECARLO_RANDOM_SIZE'])
    arrAgg_RI_Real_Random = FN.getRandomFromArr(Array=arrAgg_RI_Real, Number=AGG_Info['MONTECARLO_RANDOM_SIZE'])
    arrAgg_RI_Imag_Random = FN.getRandomFromArr(Array=arrAgg_RI_Imag, Number=AGG_Info['MONTECARLO_RANDOM_SIZE'])
    arrAgg_WLength_Random = FN.getRandomFromArr(Array=arrAgg_WLength, Number=AGG_Info['MONTECARLO_RANDOM_SIZE'])
    arrAgg_Primary_Diameter_Random = FN.getRandomFromArr(Array=arrAgg_Primary_Diameter, Number=AGG_Info['MONTECARLO_RANDOM_SIZE'])
    arrAgg_Polydispersity_Sigma_Within_Random = FN.getRandomFromArr(Array=arrAgg_Polydispersity_Sigma_Within, Number=AGG_Info['MONTECARLO_RANDOM_SIZE'])
    arrAgg_Polydispersity_Sigma_Each_Mobility_Random = FN.getRandomFromArr(Array=arrAgg_Polydispersity_Sigma_Each_Mobility, Number=AGG_Info['MONTECARLO_RANDOM_SIZE'])
    arrAgg_Prefactor_Projected_Area_Random = FN.getRandomFromArr(Array=arrAgg_Prefactor_Projected_Area, Number=AGG_Info['MONTECARLO_RANDOM_SIZE'])
    arrAgg_Exponent_Projected_Area_Random = FN.getRandomFromArr(Array=arrAgg_Exponent_Projected_Area, Number=AGG_Info['MONTECARLO_RANDOM_SIZE'])
    arrAgg_Primary_Number_Random = FN.getRandomFromArr(Array=arrAgg_Primary_Number, Number=AGG_Info['MONTECARLO_RANDOM_SIZE'])

    ################################################################################################################
    # Check Boundary for Inputs
    arrAgg_Monomer_Parameter_Random = FN.calcMonomerParameter(dpArray=arrAgg_Primary_Diameter_Random, WaveLengthArray=arrAgg_WLength_Random)

    Df_Bound = [1.801, 2.799]
    RI_Real_Bound = [1.06, 1.99]
    RI_Imag_Bound = [0.01, 0.99]
    Np_Bound = [2, 2999]
    MonomerParameter_Bound = [0.06, 0.49]

    Df_Index = FN.getGoodIndexes(Array=arrAgg_Fractal_Dimension_Random, Bound=Df_Bound)
    RI_Real_Index = FN.getGoodIndexes(Array=arrAgg_RI_Real_Random, Bound=RI_Real_Bound)
    RI_Imag_Index = FN.getGoodIndexes(Array=arrAgg_RI_Imag_Random, Bound=RI_Imag_Bound)
    Np_Index = FN.getGoodIndexes(Array=arrAgg_Primary_Number_Random, Bound=Np_Bound)
    MonomerParameter_Index = FN.getGoodIndexes(Array=arrAgg_Monomer_Parameter_Random, Bound=MonomerParameter_Bound)

    possible_Indexes = FN.findCommonIndex(Df_Index, RI_Real_Index, RI_Imag_Index, Np_Index, MonomerParameter_Index)

    Df_Possible = FN.getPossibleArray(Array=arrAgg_Fractal_Dimension_Random, Indexes=possible_Indexes)
    RI_Real_Possible = FN.getPossibleArray(Array=arrAgg_RI_Real_Random, Indexes=possible_Indexes)
    RI_Imag_Possible = FN.getPossibleArray(Array=arrAgg_RI_Imag_Random, Indexes=possible_Indexes)
    Np_Possible = FN.getPossibleArray(Array=arrAgg_Primary_Number_Random, Indexes=possible_Indexes)
    Primary_Diameter_Possible = FN.getPossibleArray(Array=arrAgg_Primary_Diameter_Random, Indexes=possible_Indexes)
    WLength_Possible = FN.getPossibleArray(Array=arrAgg_WLength_Random, Indexes=possible_Indexes)
    # to do add more
    ################################################################################################################
    # Check with Databases
    storedOutputTMatrixIndexes, plannedTMatrixInputIndexes = FN.checkTMatrixDBforIndexes(INFO=DB_Info, TableName=TMatrixTableName, Df=Df_Possible, RI_R=RI_Real_Possible, RI_I=RI_Imag_Possible, Np=Np_Possible, dp=Primary_Diameter_Possible, Wlength=WLength_Possible)
    ################################################################################################################
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
