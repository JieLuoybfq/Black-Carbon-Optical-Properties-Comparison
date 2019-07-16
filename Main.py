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
    appDirectory = GF.getRootDirectory()
    #######################
    saveHistogram = True
    ################################################################################################################
    RDG_Table_Name = "RDG_V1"
    RDG_Table_Input_Headers = ['Df', 'kf', 'R_RI', 'I_RI', 'WaveL', 'dp', 'Np', 'Version']
    RDG_Table_Output_Headers = ['RDG_ABS_CRS', 'RDG_SCA_CRS']
    DB.createTable(INFO=DB_Info, TableName=RDG_Table_Name, arrHeaderNamesInput=RDG_Table_Input_Headers, arrHeaderNamesOutput=RDG_Table_Output_Headers)
    ################################################################################################################
    TMatrix_Table_Name = "TMatrix_V1"
    TMatrix_Table_Input_Headers = ['Df', 'kf', 'R_RI', 'I_RI', 'WaveL', 'dp', 'Np', 'Version']
    TMatrix_Table_Output_Headers = ['TMatrix_ABS_CRS', 'TMatrix_SCA_CRS']
    DB.createTable(INFO=DB_Info, TableName=TMatrix_Table_Name, arrHeaderNamesInput=TMatrix_Table_Input_Headers, arrHeaderNamesOutput=TMatrix_Table_Output_Headers)
    ################################################################################################################
    Error_Table_Name = "Error_V1"
    Error_Table_Input_Headers = ['Df', 'kf', 'R_RI', 'I_RI', 'WaveL', 'dp', 'Np', 'sigma', 'Version']
    Error_Output_Headers = ['ERR_RDG_M_TMatrix_ABS_CRS', 'ERR_RDG_M_TMatrix_SCA_CRS']
    DB.createTable(INFO=DB_Info, TableName=Error_Table_Name, arrHeaderNamesInput=Error_Table_Input_Headers, arrHeaderNamesOutput=Error_Output_Headers)
    ################################################################################################################
    # DB.createDB(INFO=DB_Info)
    DB.showAllTablesInDBSummary(DB_Info)
    # DB.dropTableSet(DB_Info,RDG_Table_Name)
    # DB.dropTableSet(DB_Info, TMatrix_Table_Name)
    # DB.dropTableSet(DB_Info, Error_Table_Name)
    # DB.reinitializeDB(DB_Info)
    # DB.dumpDB(INFO=DB_Info, FileAddress=GF.getAddressTo(appDirectory, FF_Info['FOLDER_NAME_DATABASE'], FileName=GF.getDateandTimeUTC(), Extension="sql.gz"))
    # DB.loadDB(INFO=DB_Info, FileAddress=GF.findLatestFile(GF.getFilesNameAddressinFolder(GF.getAddressTo(appDirectory, FF_Info['FOLDER_NAME_DATABASE']), Extension="sql.gz")))
    # CSV_Address = GF.getAddressTo(Main=appDirectory, FolderName=FF_Info['FOLDER_NAME_CSV'], FileName=GF.getDateandTimeUTC(), Extension="csv")
    # DB.dumpTableCSV(INFO=DB_Info, TableName=TMatrix_Table_Name, Address=CSV_Address)
    ################################################################################################################
    '''
    TMatrix_Table_Input_Headers1 = ['Df', 'kf', 'R_RI', 'I_RI', 'WaveL', 'dp', 'Np', 'Version']
    TMatrix_Table_Output_Headers1 = ['TMatrix_ABS_CRS', 'TMatrix_SCA_CRS', 'MLink']
    Array1 = [[2.3, 1.2, 1.6, 0.6, 860, 33, 50, 0.1], [2.4, 1.2, 1.6, 0.6, 860, 33, 50, 0.1], [2.6, 1.2, 1.6, 0.6, 860, 33, 50, 0.1]]
    Array2 = [[0.056, 0.666, 20], [0.053, 0.666, 21], [0.051, 0.666, 22]]
    DB.insertArrayIntoTable(INFO=DB_Info, TableName=TMatrix_Table_Name, NameArray=TMatrix_Table_Input_Headers1, Array=Array1)
    DB.insertArrayIntoTable(INFO=DB_Info, TableName=TMatrix_Table_Name + "_out", NameArray=TMatrix_Table_Output_Headers1, Array=Array2)
    DB.showAllTablesInDBSummary(DB_Info)
    '''
    ################################################################################################################
    arrAgg_Fractal_Dimension = FN.createRandomNormalArr(Center=AGG_Info['AGG_FRACTAL_DIMENSION_CENTER'], Width=AGG_Info['AGG_FRACTAL_DIMENSION_STANDARD_DEVIATION'], Number=AGG_Info['MONTECARLO_ARRAY_SIZE'], digit=2)
    arrAgg_Fractal_Prefactor = FN.createRandomNormalArr(Center=AGG_Info['AGG_FRACTAL_PREFACTOR_CENTER'], Width=AGG_Info['AGG_FRACTAL_PREFACTOR_STANDARD_DEVIATION'], Number=AGG_Info['MONTECARLO_ARRAY_SIZE'], digit=2)
    arrAgg_RI_Real = FN.createRandomNormalArr(Center=AGG_Info['AGG_RI_REAL_CENTER'], Width=AGG_Info['AGG_RI_REAL_STANDARD_DEVIATION'], Number=AGG_Info['MONTECARLO_ARRAY_SIZE'], digit=2)
    arrAgg_RI_Imag = FN.createRandomNormalArr(Center=AGG_Info['AGG_RI_IMAG_CENTER'], Width=AGG_Info['AGG_RI_IMAG_STANDARD_DEVIATION'], Number=AGG_Info['MONTECARLO_ARRAY_SIZE'], digit=2)
    arrAgg_WLength = FN.createRandomNormalArr(Center=AGG_Info['AGG_WLENGTH_CENTER'], Width=AGG_Info['AGG_WLENGTH_STANDARD_DEVIATION'], Number=AGG_Info['MONTECARLO_ARRAY_SIZE'], digit=1)
    arrAgg_Primary_Diameter = FN.createRandomNormalArr(Center=AGG_Info['AGG_PRIMARY_DIAMETER_CENTER'], Width=AGG_Info['AGG_PRIMARY_DIAMETER_STANDARD_DEVIATION'], Number=AGG_Info['MONTECARLO_ARRAY_SIZE'], digit=1)
    arrAgg_Polydispersity_Sigma_Within = FN.createRandomNormalArr(Center=AGG_Info['AGG_POLYDISPERSITY_SIGMA_WITHIN_CENTER'], Width=AGG_Info['AGG_POLYDISPERSITY_SIGMA_WITHIN_STANDARD_DEVIATION'], Number=AGG_Info['MONTECARLO_ARRAY_SIZE'], digit=2)
    arrAgg_Polydispersity_Sigma_Each_Mobility = FN.createRandomNormalArr(Center=AGG_Info['AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER'], Width=AGG_Info['AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_STANDARD_DEVIATION'], Number=AGG_Info['MONTECARLO_ARRAY_SIZE'], digit=2)
    arrAgg_Prefactor_Projected_Area = FN.createRandomNormalArr(Center=AGG_Info['AGG_PREFACTOR_PROJECTED_AREA_COEFFICIENT_CENTER'], Width=AGG_Info['AGG_PREFACTOR_PROJECTED_AREA_COEFFICIENT_STANDARD_DEVIATION'], Number=AGG_Info['MONTECARLO_ARRAY_SIZE'], digit=2)
    arrAgg_Exponent_Projected_Area = FN.createRandomNormalArr(Center=AGG_Info['AGG_EXPONENT_PROJECTED_AREA_COEFFICIENT_CENTER'], Width=AGG_Info['AGG_EXPONENT_PROJECTED_AREA_COEFFICIENT_STANDARD_DEVIATION'], Number=AGG_Info['MONTECARLO_ARRAY_SIZE'], digit=2)
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

    arrDf_Bound = [1.801, 2.799]
    arrKf_Bound = [1.19, 1.21]
    arrRI_Real_Bound = [1.06, 1.99]
    arrRI_Imag_Bound = [0.01, 0.99]
    arrNp_Bound = [2, 2999]
    arrMonomerParameter_Bound = [0.06, 0.49]

    arrDf_Index = FN.getGoodIndexes(Array=arrAgg_Fractal_Dimension_Random, Bound=arrDf_Bound)
    arrKf_Index = FN.getGoodIndexes(Array=arrAgg_Fractal_Prefactor_Random, Bound=arrKf_Bound)
    arrRI_Real_Index = FN.getGoodIndexes(Array=arrAgg_RI_Real_Random, Bound=arrRI_Real_Bound)
    arrRI_Imag_Index = FN.getGoodIndexes(Array=arrAgg_RI_Imag_Random, Bound=arrRI_Imag_Bound)
    arrNp_Index = FN.getGoodIndexes(Array=arrAgg_Primary_Number_Random, Bound=arrNp_Bound)
    arrMonomerParameter_Index = FN.getGoodIndexes(Array=arrAgg_Monomer_Parameter_Random, Bound=arrMonomerParameter_Bound)

    arrPossible_Indexes = FN.findCommonIndex(arrDf_Index, arrKf_Index, arrRI_Real_Index, arrRI_Imag_Index, arrNp_Index, arrMonomerParameter_Index)

    arrDf_Possible = FN.getPossibleArray(Array=arrAgg_Fractal_Dimension_Random, Indexes=arrPossible_Indexes)
    arrKf_Possible = FN.getPossibleArray(Array=arrAgg_Fractal_Prefactor_Random, Indexes=arrPossible_Indexes)
    arrRI_Real_Possible = FN.getPossibleArray(Array=arrAgg_RI_Real_Random, Indexes=arrPossible_Indexes)
    arrRI_Imag_Possible = FN.getPossibleArray(Array=arrAgg_RI_Imag_Random, Indexes=arrPossible_Indexes)
    arrNp_Possible = FN.getPossibleArray(Array=arrAgg_Primary_Number_Random, Indexes=arrPossible_Indexes)
    arrPrimary_Diameter_Possible = FN.getPossibleArray(Array=arrAgg_Primary_Diameter_Random, Indexes=arrPossible_Indexes)
    arrWLength_Possible = FN.getPossibleArray(Array=arrAgg_WLength_Random, Indexes=arrPossible_Indexes)
    #### to do add more

    ################################################################################################################
    # Check with Databases
    arrVersion_Array = FN.createConstantArray(Number=Version, Howmany=len(arrPossible_Indexes))
    arrIndex_Array = FN.createIndexArray(0, len=len(arrVersion_Array))
    # ['Df', 'kf', 'R_RI', 'I_RI', 'WaveL', 'dp', 'Np', 'Version']
    TMatrix_Main_Input_Array_Indexed = FN.joinArray(arrIndex_Array, arrDf_Possible, arrKf_Possible, arrRI_Real_Possible, arrRI_Imag_Possible, arrWLength_Possible, arrPrimary_Diameter_Possible, arrNp_Possible, arrVersion_Array)
    TMatrix_Main_Input_Array = FN.joinArray(arrDf_Possible, arrKf_Possible, arrRI_Real_Possible, arrRI_Imag_Possible, arrWLength_Possible, arrPrimary_Diameter_Possible, arrNp_Possible, arrVersion_Array)
    # ['Df', 'kf', 'R_RI', 'I_RI', 'WaveL', 'dp', 'Np', 'Version']
    RDG_Main_Input_Array_Indexed = FN.joinArray(arrIndex_Array, arrDf_Possible, arrKf_Possible, arrRI_Real_Possible, arrRI_Imag_Possible, arrWLength_Possible, arrPrimary_Diameter_Possible, arrNp_Possible, arrVersion_Array)
    RDG_Main_Input_Array = FN.joinArray(arrDf_Possible, arrKf_Possible, arrRI_Real_Possible, arrRI_Imag_Possible, arrWLength_Possible, arrPrimary_Diameter_Possible, arrNp_Possible, arrVersion_Array)
    '''
    TMatrix_Main_Input_Array.append([2.3, 1.2, 1.6, 0.6, 860, 33, 50, 0.1])
    TMatrix_Main_Input_Array.append([2.4, 1.2, 1.6, 0.6, 860, 33, 50, 0.1])
    TMatrix_Main_Input_Array.append([2.6, 1.2, 1.6, 0.6, 860, 33, 50, 0.1])
    TMatrix_Main_Input_Array.append([2.1, 1.2, 1.6, 0.6, 860, 33, 50, 0.1])
    '''
    TMatrix_DB_Input_Found, TMatrix_DB_Output_Found, TMatrix_Planned_Input, TMatrix_New = FN.checkMethodDBforTMatrixIndexes(INFO=DB_Info, TableName=TMatrix_Table_Name, Header=TMatrix_Table_Input_Headers, Array=TMatrix_Main_Input_Array)
    # storedInput_RDG, storedOutput_RDG, plannedInput_RDG = FN.checkMethodDBforTMatrixIndexes(INFO=DB_Info, TableName=RDG_Table_Name, Header=RDG_Table_Input_Headers, Array=RDG_Main_Input_Array)
    ################################################################################################################

    TMatrix_DB_Main_TableName = 'Raw_V1'

    TMatrix_DB_Main_Column_Name, TMatrix_DB_Main_Data_Full = DB.readAllRowsfromTable(INFO=DB_Info, TableName=TMatrix_DB_Main_TableName)
    TMatrix_DB_Main_Column_Name = GF.selectColumnsList(ColumnIndex=[1, 2, 3, 4, 5], List=TMatrix_DB_Main_Column_Name, Dimension=1)

    TMatrix_DB_Main_SCT_Coeff_Full = GF.selectColumnsList(ColumnIndex=[7], List=TMatrix_DB_Main_Data_Full)
    TMatrix_DB_Main_ABS_Coeff_Full = GF.selectColumnsList(ColumnIndex=[8], List=TMatrix_DB_Main_Data_Full)
    TMatrix_DB_Main_Data_Full = GF.selectColumnsList(ColumnIndex=[1, 2, 3, 4, 5], List=TMatrix_DB_Main_Data_Full)

    TMatrix_DB_Main_Unique_Values = FN.uniqueEntry(TMatrix_DB_Main_Data_Full)
    ################################################################################################################
    division = 4
    TMatrix_Planned_Input_Chopped = GF.divideArray(NumberofDivisions=division, List=TMatrix_Planned_Input)
    func = partial(TM.TmatrixInterpolator, TMatrix_DB_Main_Data_Full, TMatrix_DB_Main_Unique_Values, TMatrix_DB_Main_ABS_Coeff_Full, TMatrix_DB_Main_SCT_Coeff_Full)
    arrInTMatrix = []
    arrOutTMatrix = []
    with ThreadPoolExecutor() as executor:
        for inputs, outputs in zip(TMatrix_Planned_Input_Chopped, executor.map(func, TMatrix_Planned_Input_Chopped)):
            arrInTMatrix.append(inputs)
            arrOutTMatrix.append(outputs)
    ####################################
    TMatrix_Interpolation_Input = []
    TMatrix_Interpolation_Output = []
    for i in range(division):
        for j in range(len(arrInTMatrix[i])):
            TMatrix_Interpolation_Input.append(arrInTMatrix[i][j])
        for j in range(len(arrOutTMatrix[i][0])):
            arrT = []
            arrT.append(arrOutTMatrix[i][0][j])
            arrT.append(arrOutTMatrix[i][1][j])
            TMatrix_Interpolation_Output.append(arrT)
    ####################################
    if TMatrix_Interpolation_Input and TMatrix_Interpolation_Output:
        LastID = DB.insertArrayIntoTable(INFO=DB_Info, TableName=TMatrix_Table_Name, NameArray=TMatrix_Table_Input_Headers, giveID=True, Array=TMatrix_Interpolation_Input)
        TMatrix_Interpolation_Output_M = FN.addMlinkToArray(Array=TMatrix_Interpolation_Output, LastID=LastID)
        TMatrix_Table_Output_Headers.append('MLink')
        DB.insertArrayIntoTable(INFO=DB_Info, TableName=TMatrix_Table_Name + "_Out", NameArray=TMatrix_Table_Output_Headers, Array=TMatrix_Interpolation_Output_M)
    DB.showAllTablesInDBSummary(DB_Info)
    ################################################################################################################
    TMatrix_Final_Input = []
    TMatrix_Final_Output = []
    old_Counter = 0
    new_Counter = 0

    if TMatrix_Interpolation_Input == TMatrix_Planned_Input:
        for i in range(len(TMatrix_New)):
            if TMatrix_New[i] == -1:
                TMatrix_Final_Input.append(TMatrix_Interpolation_Input[new_Counter])
                TMatrix_Final_Output.append(TMatrix_Interpolation_Output[new_Counter])
                new_Counter += 1
            elif TMatrix_New[i] == 1:
                TMatrix_Final_Input.append(TMatrix_DB_Input_Found[old_Counter])
                TMatrix_Final_Output.append(TMatrix_DB_Output_Found[old_Counter])
                old_Counter += 1
    else:
        raise Exception('change in TMatrix input: Interpolation!')

    if TMatrix_Main_Input_Array != TMatrix_Final_Input:
        raise Exception('change in TMatrix input: Database!')

    ################################################################################################################
    A = 51
    ################################################################################################################
    # TMatrix_ABS_CS, TMatrix_SCA_CS = TM.TmatrixInterpolator(FullMainDB=TMatrixDBInputData_Full, MainDBUniques=TMatrixDBUniqueColumn, ABS_MainDB=TMatrixDBABS_Coeff_Full, SCA_MainDB=TMatrixDBSCT_Coeff_Full, TargetArray=plannedInput_TMatrix)
    ################################################################################################################
