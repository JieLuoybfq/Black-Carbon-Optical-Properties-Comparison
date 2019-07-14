# Programmed by Keyhan Babaee Under Prof. Steven Rogak supervision, https://github.com/KeyhanB
# Using Fortran Code from: C. Liu, X. Xu, Y. Yin, M. Schnaiter, Y. L. Yung, 2018: Black carbon aggregate: A database for optical properties
Version = 0.1
# July 2019
import ConfigParserM as CP
import DBMethods as DB
from ConfigParserM import logging
import GeneralFunctions as GF
import AppFunctions as FN
import TMatrix as TM
from concurrent.futures import ThreadPoolExecutor
from functools import partial

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
    RDGInputHeaders = ['Df', 'kf', 'R_RI', 'I_RI', 'WaveL', 'dp', 'Np', 'Version']
    RDGOutputHeaders = ['RDG_ABS_CRS', 'RDG_SCA_CRS']
    DB.createTable(INFO=DB_Info, TableName=RDGTableName, arrHeaderNamesInput=RDGInputHeaders, arrHeaderNamesOutput=RDGOutputHeaders)
    ################################################################################################################
    TMatrixTableName = "TMatrix_V1"
    TMatrixInputHeaders = ['Df', 'kf', 'R_RI', 'I_RI', 'WaveL', 'dp', 'Np', 'Version']
    TMatrixOutputHeaders = ['TMatrix_ABS_CRS', 'TMatrix_SCA_CRS']
    DB.createTable(INFO=DB_Info, TableName=TMatrixTableName, arrHeaderNamesInput=TMatrixInputHeaders, arrHeaderNamesOutput=TMatrixOutputHeaders)
    ################################################################################################################
    ErrorTableName = "Error_V1"
    ErrorInputHeaders = ['Df', 'kf', 'R_RI', 'I_RI', 'WaveL', 'dp', 'Np', 'sigma', 'Version']
    ErrorOutputHeaders = ['ERR_RDG_M_TMatrix_ABS_CRS', 'ERR_RDG_M_TMatrix_SCA_CRS']
    DB.createTable(INFO=DB_Info, TableName=ErrorTableName, arrHeaderNamesInput=ErrorInputHeaders, arrHeaderNamesOutput=ErrorOutputHeaders)
    ################################################################################################################
    # DB.createDB(INFO=DB_Info)
    DB.showAllTablesInDBSummary(DB_Info)
    # DB.dropTableSet(DB_Info,RDGTableName)
    # DB.dropTableSet(DB_Info, TMatrixTableName)
    # DB.dropTableSet(DB_Info, ErrorTableName)
    # DB.reinitializeDB(DB_Info)
    ################################################################################################################
    '''
    TMatrixInputHeaders = ['Df', 'kf', 'R_RI', 'I_RI', 'WaveL', 'dp', 'Np', 'Version']
    TMatrixOutputHeaders = ['TMatrix_ABS_CRS', 'TMatrix_SCA_CRS','MLink']
    Array1=[[2.3,1.2,1.6,0.6,860,33,50,0.1],[2.4,1.2,1.6,0.6,860,33,50,0.1],[2.6,1.2,1.6,0.6,860,33,50,0.1]]
    Array2=[[0.056,0.666,1],[0.053,0.666,2],[0.051,0.666,3]]
    DB.insertArrayIntoTable(INFO=DB_Info,TableName=TMatrixTableName,NameArray=TMatrixInputHeaders,Array=Array1)
    DB.insertArrayIntoTable(INFO=DB_Info, TableName=TMatrixTableName+"_out", NameArray=TMatrixOutputHeaders, Array=Array2)
    DB.showAllTablesInDBSummary(DB_Info)
    '''
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
    Kf_Bound = [1.19, 1.21]
    RI_Real_Bound = [1.06, 1.99]
    RI_Imag_Bound = [0.01, 0.99]
    Np_Bound = [2, 2999]
    MonomerParameter_Bound = [0.06, 0.49]

    Df_Index = FN.getGoodIndexes(Array=arrAgg_Fractal_Dimension_Random, Bound=Df_Bound)
    Kf_Index = FN.getGoodIndexes(Array=arrAgg_Fractal_Prefactor_Random, Bound=Kf_Bound)
    RI_Real_Index = FN.getGoodIndexes(Array=arrAgg_RI_Real_Random, Bound=RI_Real_Bound)
    RI_Imag_Index = FN.getGoodIndexes(Array=arrAgg_RI_Imag_Random, Bound=RI_Imag_Bound)
    Np_Index = FN.getGoodIndexes(Array=arrAgg_Primary_Number_Random, Bound=Np_Bound)
    MonomerParameter_Index = FN.getGoodIndexes(Array=arrAgg_Monomer_Parameter_Random, Bound=MonomerParameter_Bound)

    possible_Indexes = FN.findCommonIndex(Df_Index, Kf_Index, RI_Real_Index, RI_Imag_Index, Np_Index, MonomerParameter_Index)

    Df_Possible = FN.getPossibleArray(Array=arrAgg_Fractal_Dimension_Random, Indexes=possible_Indexes)
    Kf_Possible = FN.getPossibleArray(Array=arrAgg_Fractal_Prefactor_Random, Indexes=possible_Indexes)
    RI_Real_Possible = FN.getPossibleArray(Array=arrAgg_RI_Real_Random, Indexes=possible_Indexes)
    RI_Imag_Possible = FN.getPossibleArray(Array=arrAgg_RI_Imag_Random, Indexes=possible_Indexes)
    Np_Possible = FN.getPossibleArray(Array=arrAgg_Primary_Number_Random, Indexes=possible_Indexes)
    Primary_Diameter_Possible = FN.getPossibleArray(Array=arrAgg_Primary_Diameter_Random, Indexes=possible_Indexes)
    WLength_Possible = FN.getPossibleArray(Array=arrAgg_WLength_Random, Indexes=possible_Indexes)
    # to do add more
    ################################################################################################################
    # Check with Databases
    Version_Array = FN.createConstantArray(Number=Version, Howmany=len(possible_Indexes))
    # ['Df', 'kf', 'R_RI', 'I_RI', 'WaveL', 'dp', 'Np', 'Version']
    TMatrix_Input_Array = FN.joinArray(Df_Possible, Kf_Possible, RI_Real_Possible, RI_Imag_Possible, WLength_Possible, Primary_Diameter_Possible, Np_Possible, Version_Array)
    # ['Df', 'kf', 'R_RI', 'I_RI', 'WaveL', 'dp', 'Np', 'Version']
    RDG_Input_Array = FN.joinArray(Df_Possible, Kf_Possible, RI_Real_Possible, RI_Imag_Possible, WLength_Possible, Primary_Diameter_Possible, Np_Possible, Version_Array)

    storedInput_TMatrix, storedOutput_TMatrix, plannedInput_TMatrix = FN.checkMethodDBforIndexes(INFO=DB_Info, TableName=TMatrixTableName, Header=TMatrixInputHeaders, Array=TMatrix_Input_Array)
    storedInput_RDG, storedOutput_RDG, plannedInput_RDG = FN.checkMethodDBforIndexes(INFO=DB_Info, TableName=RDGTableName, Header=RDGInputHeaders, Array=RDG_Input_Array)
    ################################################################################################################

    TMatrixDBTableName = 'Raw_V1'

    TMatrixDBColumnName, TMatrixDBDataFull = DB.readAllRowsfromTable(INFO=DB_Info, TableName=TMatrixDBTableName)
    TMatrixDBColumnName = GF.selectColumnsList(ColumnIndex=[1, 2, 3, 4, 5], List=TMatrixDBColumnName, Dimension=1)

    TMatrixDBSCT_Coeff_Full = GF.selectColumnsList(ColumnIndex=[7], List=TMatrixDBDataFull)
    TMatrixDBABS_Coeff_Full = GF.selectColumnsList(ColumnIndex=[8], List=TMatrixDBDataFull)
    TMatrixDBInputData_Full = GF.selectColumnsList(ColumnIndex=[1, 2, 3, 4, 5], List=TMatrixDBDataFull)

    TMatrixDBUniqueColumn = FN.uniqueEntry(TMatrixDBInputData_Full)
    ################################################################################################################
    Division = 4
    plannedInput_TMatrix_Chopped = GF.divideArray(NumberofDivisions=Division, List=plannedInput_TMatrix)
    func = partial(TM.TmatrixInterpolator, TMatrixDBInputData_Full, TMatrixDBUniqueColumn, TMatrixDBABS_Coeff_Full, TMatrixDBSCT_Coeff_Full)
    arrInTMatrix = []
    arrOutTMatrix = []
    with ThreadPoolExecutor() as executor:
        for inputs, outputs in zip(plannedInput_TMatrix_Chopped, executor.map(func, plannedInput_TMatrix_Chopped)):
            arrInTMatrix.append(inputs)
            arrOutTMatrix.append(outputs)
    ####################################
    plannedInput_TMatrix_Out = []
    plannedOutput_TMatrix_Out = []
    for i in range(Division):
        for j in range(len(arrInTMatrix[i])):
            plannedInput_TMatrix_Out.append(arrInTMatrix[i][j])

        for j in range(len(arrOutTMatrix[i][0])):
            A = []
            A.append(arrOutTMatrix[i][0][j])
            A.append(arrOutTMatrix[i][1][j])
            plannedOutput_TMatrix_Out.append(A)
    ################################################################################################################
    A = 3
    ################################################################################################################
    # TMatrix_ABS_CS, TMatrix_SCA_CS = TM.TmatrixInterpolator(FullMainDB=TMatrixDBInputData_Full, MainDBUniques=TMatrixDBUniqueColumn, ABS_MainDB=TMatrixDBABS_Coeff_Full, SCA_MainDB=TMatrixDBSCT_Coeff_Full, TargetArray=plannedInput_TMatrix)
    ################################################################################################################
    '''
    TMatrixInputHeaders = ['Df', 'kf', 'R_RI', 'I_RI', 'WaveL', 'dp', 'Np', 'Version']
    TMatrixOutputHeaders = ['TMatrix_ABS_CRS', 'TMatrix_SCA_CRS', 'MLink']
    C = []
    counter = 1
    for i in range(len(plannedInput_TMatrix)):
        C.append([0.056 + counter * 0.001, 0.666, counter])
        counter += 1
    DB.insertArrayIntoTable(INFO=DB_Info, TableName=TMatrixTableName, NameArray=TMatrixInputHeaders, Array=plannedInput_TMatrix)
    DB.insertArrayIntoTable(INFO=DB_Info, TableName=TMatrixTableName + "_out", NameArray=TMatrixOutputHeaders, Array=C)
    DB.showAllTablesInDBSummary(DB_Info)

    TMatrix_Input_Array.append([2.3, 1.2, 1.6, 0.6, 860, 33, 50, 0.1])
    TMatrix_Input_Array.append([2.4, 1.2, 1.6, 0.6, 860, 33, 50, 0.1])
    TMatrix_Input_Array.append([2.6, 1.2, 1.6, 0.6, 860, 33, 50, 0.1])
    storedOutput_TMatrix1, plannedInput_TMatrix1 = FN.checkMethodDBforIndexes(INFO=DB_Info, TableName=TMatrixTableName, Header=TMatrixInputHeaders, Array=TMatrix_Input_Array)
    '''
    ################################################################################################################
