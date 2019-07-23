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
import FSAC_RDG as FRDG
from decimal import Decimal
import pandas as pd

if __name__ == "__main__":
    CP.readLogConfig()
    logging.info("App Started.")
    DB_Info = CP.readConfigToDict(SectionName="DatabaseInfo")
    FF_Info = CP.readConfigToDict(SectionName="FilesFoldersInfo")
    AGG_Info = CP.readConfigToDict(SectionName="AggregateDetails", ConvertParseTo='float', hasComment=True)
    logging.info("config retrieved.")
    #######################
    appDirectory = GF.getRootDirectory()
    #######################
    saveHistogram = True
    ##############################################
    RDG_Table_Name = "RDG_V1"
    RDG_Table_Input_Headers = ['Df', 'kf', 'R_RI', 'I_RI', 'WaveL', 'dp', 'Np', 'Version']
    RDG_Table_Output_Headers = ['RDG_ABS_CRS', 'RDG_SCA_CRS']
    DB.createTable(INFO=DB_Info, TableName=RDG_Table_Name, arrHeaderNamesInput=RDG_Table_Input_Headers, arrHeaderNamesOutput=RDG_Table_Output_Headers)
    ##############################################
    TMatrix_Table_Name = "TMatrix_V1"
    TMatrix_Table_Input_Headers = ['Df', 'kf', 'R_RI', 'I_RI', 'WaveL', 'dp', 'Np', 'Version']
    TMatrix_Table_Output_Headers = ['TMatrix_ABS_CRS', 'TMatrix_SCA_CRS']
    DB.createTable(INFO=DB_Info, TableName=TMatrix_Table_Name, arrHeaderNamesInput=TMatrix_Table_Input_Headers, arrHeaderNamesOutput=TMatrix_Table_Output_Headers)
    ##############################################
    Error_Table_Name = "Error_V1"
    Error_Table_Input_Headers = ['Df', 'kf', 'R_RI', 'I_RI', 'WaveL', 'dp', 'Np', 'Version']
    Error_Table_Output_Headers = ['ERR_RDG_M_TMatrix_ABS_CRS', 'ERR_RDG_M_TMatrix_SCA_CRS', 'ERR_RDG_M_TMatrix_ABS_CRS_Percent', 'ERR_RDG_M_TMatrix_SCA_CRS_Percent']
    DB.createTable(INFO=DB_Info, TableName=Error_Table_Name, arrHeaderNamesInput=Error_Table_Input_Headers, arrHeaderNamesOutput=Error_Table_Output_Headers)
    ################################################################################################################ DB Control
    # DB.createDB(INFO=DB_Info)
    DB.showAllTablesInDBSummary(DB_Info)
    # DB.dropTableSet(DB_Info, RDG_Table_Name)
    # DB.dropTableSet(DB_Info, TMatrix_Table_Name)
    # DB.dropTableSet(DB_Info, Error_Table_Name)
    # DB.reinitializeDB(DB_Info)
    # DB.dumpDB(INFO=DB_Info, FileAddress=GF.getAddressTo(appDirectory, FF_Info['FOLDER_NAME_DATABASE'], FileName=GF.getDateandTimeUTC(), Extension="sql.gz"))
    # DB.loadDB(INFO=DB_Info, FileAddress=GF.findLatestFile(GF.getFilesNameAddressinFolder(GF.getAddressTo(appDirectory, FF_Info['FOLDER_NAME_DATABASE']), Extension="sql.gz")))
    # DB.dumpTableSetCSV(INFO=DB_Info, TableName=RDG_Table_Name, AddressMain=appDirectory)
    # DB.dumpTableSetCSV(INFO=DB_Info, TableName=TMatrix_Table_Name, AddressMain=appDirectory)
    # DB.dumpTableSetCSV(INFO=DB_Info, TableName=Error_Table_Name, AddressMain=appDirectory)
    ################################################################################################################
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
    logging.info("Random parameter generated.")
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
    #### to do add more if needed
    logging.info("Boundary issued.")
    ################################################################################################################
    # Check with Databases
    arrVersion_Array = FN.createConstantArray(Number=Version, Howmany=len(arrPossible_Indexes))
    TMatrix_Main_Input_Array = FN.joinArray(arrDf_Possible, arrKf_Possible, arrRI_Real_Possible, arrRI_Imag_Possible, arrWLength_Possible, arrPrimary_Diameter_Possible, arrNp_Possible, arrVersion_Array)
    RDG_Main_Input_Array = FN.joinArray(arrDf_Possible, arrKf_Possible, arrRI_Real_Possible, arrRI_Imag_Possible, arrWLength_Possible, arrPrimary_Diameter_Possible, arrNp_Possible, arrVersion_Array)

    TMatrix_DB_Input_Found, TMatrix_DB_Output_Found, TMatrix_Planned_Input, TMatrix_New = FN.checkMethodDBforTMatrixIndexes(INFO=DB_Info, TableName=TMatrix_Table_Name, Header=TMatrix_Table_Input_Headers, Array=TMatrix_Main_Input_Array)
    RDG_DB_Input_Found, RDG_DB_Output_Found, RDG_Planned_Input, RDG_New = FN.checkMethodDBforRDGIndexes(INFO=DB_Info, TableName=RDG_Table_Name, Header=RDG_Table_Input_Headers, Array=RDG_Main_Input_Array)
    logging.info("DB old and new ones retrieved.")
    ################################################################################################################
    ################################################################################################################
    ################################################################################################################
    ################################################################################################################
    ################################################################################################################
    logging.info("RDG started.")
    RDG_Planned_Output = []
    for i in range(len(RDG_Planned_Input)):
        arrRDG = []
        rDG_ABS, rDG_SCA = FRDG.FSAC_RDG(Df=RDG_Planned_Input[i][0], kf=RDG_Planned_Input[i][1], R_RI=RDG_Planned_Input[i][2], I_RL=RDG_Planned_Input[i][3], WaveL_nm=RDG_Planned_Input[i][4], dp=RDG_Planned_Input[i][5], Np=RDG_Planned_Input[i][6])
        arrRDG.append(rDG_ABS)
        arrRDG.append(rDG_SCA)
        RDG_Planned_Output.append(arrRDG)

    ####################################
    if RDG_Planned_Input and RDG_Planned_Output:
        LastID = DB.insertArrayIntoTable(INFO=DB_Info, TableName=RDG_Table_Name, NameArray=RDG_Table_Input_Headers, giveID=True, Array=RDG_Planned_Input)
        RDG_Planned_Output_M = FN.addMlinkToArray(Array=RDG_Planned_Output, LastID=LastID)
        RDG_Table_Output_Headers.append('MLink')
        DB.insertArrayIntoTable(INFO=DB_Info, TableName=RDG_Table_Name + "_Out", NameArray=RDG_Table_Output_Headers, Array=RDG_Planned_Output_M)
    DB.showAllTablesInDBSummary(DB_Info)
    ################################################################################################################
    RDG_Final_Input = []
    RDG_Final_Output = []
    old_Counter = 0
    new_Counter = 0
    for i in range(len(RDG_New)):
        if RDG_New[i] == -1:
            RDG_Final_Input.append(RDG_Planned_Input[new_Counter])
            RDG_Final_Output.append(RDG_Planned_Output[new_Counter])
            new_Counter += 1
        elif RDG_New[i] == 1:
            RDG_Final_Input.append(RDG_DB_Input_Found[old_Counter])
            RDG_Final_Output.append(RDG_DB_Output_Found[old_Counter])
            old_Counter += 1

    if RDG_Main_Input_Array != RDG_Final_Input:
        logging.error("change in RDG input: Database: " + "\n" + "RDG_Main_Input_Array---" + str(RDG_Main_Input_Array) + "\n" + "RDG_Final_Input---" + str(RDG_Final_Input))
        raise Exception('change in RDG input: Database')
    logging.info("RDG finished retrieved.")
    ################################################################################################################
    ################################################################################################################
    ################################################################################################################
    ################################################################################################################
    ################################################################################################################
    logging.info("T-Matrix started.")
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
        logging.error("change in T-Matrix input: Interpolation: " + "\n" + "TMatrix_Interpolation_Input---" + str(TMatrix_Interpolation_Input) + "\n" + "TMatrix_Planned_Input---" + str(TMatrix_Planned_Input))
        raise Exception('change in TMatrix input: Interpolation')

    if TMatrix_Main_Input_Array != TMatrix_Final_Input:
        logging.error("change in T-Matrix input: Database: " + "\n" + "TMatrix_Main_Input_Array---" + str(TMatrix_Main_Input_Array) + "\n" + "TMatrix_Final_Input---" + str(TMatrix_Final_Input))
        raise Exception('change in T-Matrix input: Database')
    logging.info("T-Matrix Finished.")
    ################################################################################################################
    ################################################################################################################
    ################################################################################################################
    ################################################################################################################
    ################################################################################################################
    if RDG_Final_Input != TMatrix_Final_Input:
        logging.error("discrepancy between RDG and T-Matrix inputs: Array " + "\n" + "RDG_Final_Input---" + str(RDG_Final_Input) + "\n" + "TMatrix_Final_Input---" + str(TMatrix_Final_Input))
        raise Exception('discrepancy between RDG and T-Matrix inputs: Array')
    ################################################################################################################
    ################################################################################################################
    ################################################################################################################
    ################################################################################################################
    ################################################################################################################
    logging.info("Error calculation started.")
    Error_Final_Input = []
    Error_Final_Output = []
    for i in range(len(RDG_Final_Output)):
        arrErr = []
        arrErr.append(RDG_Final_Output[i][0] - TMatrix_Final_Output[i][0])
        arrErr.append(RDG_Final_Output[i][1] - TMatrix_Final_Output[i][1])
        arrErr.append((RDG_Final_Output[i][0] - TMatrix_Final_Output[i][0]) * Decimal(100) / TMatrix_Final_Output[i][0])
        arrErr.append((RDG_Final_Output[i][1] - TMatrix_Final_Output[i][1]) * Decimal(100) / TMatrix_Final_Output[i][1])
        Error_Final_Output.append(arrErr)
        Error_Final_Input.append(RDG_Final_Input[i][:])
    ####################################
    if Error_Final_Input and Error_Final_Output:
        LastID = DB.insertArrayIntoTable(INFO=DB_Info, TableName=Error_Table_Name, NameArray=Error_Table_Input_Headers, giveID=True, Array=Error_Final_Input)
        Error_Final_Output_M = FN.addMlinkToArray(Array=Error_Final_Output, LastID=LastID)
        Error_Table_Output_Headers.append('MLink')
        DB.insertArrayIntoTable(INFO=DB_Info, TableName=Error_Table_Name + "_Out", NameArray=Error_Table_Output_Headers, Array=Error_Final_Output_M)
    DB.showAllTablesInDBSummary(DB_Info)
    logging.info("Error calculation finished.")
    ################################################################################################################
    ################################################################################################################
    ################################################################################################################
    ################################################################################################################
    ################################################################################################################
    logging.info("Plotting error calculation started.")

    # Selected Error array
    address_Graph_Real_Selected = GF.getAddressTo(Main=appDirectory, FolderName=FF_Info['FOLDER_NAME_GRAPH'], FileName="Selected_Points_Real_Error", Extension="jpg")
    address_Graph_Percent_Selected = GF.getAddressTo(Main=appDirectory, FolderName=FF_Info['FOLDER_NAME_GRAPH'], FileName="Selected_Points_Percent_Error", Extension="jpg")
    address_Graph_Real_Selected_Scatter = GF.getAddressTo(Main=appDirectory, FolderName=FF_Info['FOLDER_NAME_GRAPH'], FileName="Selected_Points_Real_Error_ScatterMatrix", Extension="jpg")
    address_Graph_Percent_Selected_Scatter = GF.getAddressTo(Main=appDirectory, FolderName=FF_Info['FOLDER_NAME_GRAPH'], FileName="Selected_Points_Percent_Error_ScatterMatrix", Extension="jpg")

    Data_Frame_Error_Real_Selected_Arr, n1 = FN.joinColumnsToArray(Array=Error_Final_Input, ArrtobeJoined=Error_Final_Output, ColumnIndexes=[0, 1])
    Data_Frame_Error_Percent_Selected_Arr, n2 = FN.joinColumnsToArray(Array=Error_Final_Input, ArrtobeJoined=Error_Final_Output, ColumnIndexes=[2, 3])
    n1.append("Absorption Difference(um" + "$^{}$".format(2) + ")")
    n1.append("Scattering Difference(um" + "$^{}$".format(2) + ")")
    n2.append("Absorption Difference(%)")
    n2.append("Scattering Difference(%)")
    Data_Frame_Error_Real_Selected_Object = pd.DataFrame(Data_Frame_Error_Real_Selected_Arr, columns=n1, dtype=float)
    Data_Frame_Error_Percent_Selected_Object = pd.DataFrame(Data_Frame_Error_Percent_Selected_Arr, columns=n2, dtype=float)
    FN.Fig_Plot_Save_Scatterplot_Matrix(Address=address_Graph_Real_Selected_Scatter, Dataframe=Data_Frame_Error_Real_Selected_Object)
    FN.Fig_Plot_Save_Scatterplot_Matrix(Address=address_Graph_Percent_Selected_Scatter, Dataframe=Data_Frame_Error_Percent_Selected_Object)

    ABS_Error_Real_Selected = GF.selectColumnsList(ColumnIndex=[0], List=Error_Final_Output)
    SCA_Error_Real_Selected = GF.selectColumnsList(ColumnIndex=[1], List=Error_Final_Output)
    ABS_Error_Percent_Selected = GF.selectColumnsList(ColumnIndex=[2], List=Error_Final_Output)
    SCA_Error_Percent_Selected = GF.selectColumnsList(ColumnIndex=[3], List=Error_Final_Output)
    FN.Fig_Plot_Save_Scatter_X_Linear_Y_Linear(Address=address_Graph_Real_Selected, X_Array=ABS_Error_Real_Selected, Y_array=SCA_Error_Real_Selected, X_Label="Absorption Difference(um" + "$^{}$".format(2) + ")", Y_label="Scattering Difference(um" + "$^{}$".format(2) + ")", Plot_Title="RDG and T-Matrix Real Difference")
    FN.Fig_Plot_Save_Scatter_X_Linear_Y_Linear(Address=address_Graph_Percent_Selected, X_Array=ABS_Error_Percent_Selected, Y_array=SCA_Error_Percent_Selected, tickLabelStyle='plain', X_Label="Absorption Difference (%)", Y_label="Scattering Difference (%)", Plot_Title="RDG and T-Matrix Percentage Difference")
    ################################################
    # Full Error Array
    address_Graph_Real_Full = GF.getAddressTo(Main=appDirectory, FolderName=FF_Info['FOLDER_NAME_GRAPH'], FileName="Full_Points_Real_Error", Extension="jpg")
    address_Graph_Percent_Full = GF.getAddressTo(Main=appDirectory, FolderName=FF_Info['FOLDER_NAME_GRAPH'], FileName="Full_Points_Percent_Error", Extension="jpg")
    address_Graph_Real_Full_Scatter = GF.getAddressTo(Main=appDirectory, FolderName=FF_Info['FOLDER_NAME_GRAPH'],
                                                      FileName="Full_Points_Real_Error_ScatterMatrix", Extension="jpg")
    address_Graph_Percent_Full_Scatter = GF.getAddressTo(Main=appDirectory, FolderName=FF_Info['FOLDER_NAME_GRAPH'],
                                                         FileName="Full_Points_Percent_Error_ScatterMatrix", Extension="jpg")

    Error_Table_Df = DB.readColwithColumnName(INFO=DB_Info, TableName=Error_Table_Name, ColumnName="Df")
    Error_Table_kf = DB.readColwithColumnName(INFO=DB_Info, TableName=Error_Table_Name, ColumnName="kf")
    Error_Table_R_RI = DB.readColwithColumnName(INFO=DB_Info, TableName=Error_Table_Name, ColumnName="R_RI")
    Error_Table_I_RI = DB.readColwithColumnName(INFO=DB_Info, TableName=Error_Table_Name, ColumnName="I_RI")
    Error_Table_WaveL = DB.readColwithColumnName(INFO=DB_Info, TableName=Error_Table_Name, ColumnName="WaveL")
    Error_Table_dp = DB.readColwithColumnName(INFO=DB_Info, TableName=Error_Table_Name, ColumnName="dp")
    Error_Table_Np = DB.readColwithColumnName(INFO=DB_Info, TableName=Error_Table_Name, ColumnName="Np")
    Error_In = FN.joinArray(Error_Table_Df, Error_Table_kf, Error_Table_R_RI, Error_Table_I_RI, Error_Table_WaveL, Error_Table_dp, Error_Table_Np)

    ABS_Error_Real_Full = DB.readColwithColumnName(INFO=DB_Info, TableName=Error_Table_Name + "_Out", ColumnName="ERR_RDG_M_TMatrix_ABS_CRS")
    SCA_Error_Real_Full = DB.readColwithColumnName(INFO=DB_Info, TableName=Error_Table_Name + "_Out", ColumnName="ERR_RDG_M_TMatrix_SCA_CRS")
    ABS_Error_Percent_Full = DB.readColwithColumnName(INFO=DB_Info, TableName=Error_Table_Name + "_Out", ColumnName="ERR_RDG_M_TMatrix_ABS_CRS_Percent")
    SCA_Error_Percent_Full = DB.readColwithColumnName(INFO=DB_Info, TableName=Error_Table_Name + "_Out", ColumnName="ERR_RDG_M_TMatrix_SCA_CRS_Percent")
    Error_Out = FN.joinArray(ABS_Error_Real_Full, SCA_Error_Real_Full, ABS_Error_Percent_Full, SCA_Error_Percent_Full)

    Data_Frame_Error_Real_Full_Arr, n1 = FN.joinColumnsToArray(Array=Error_In, ArrtobeJoined=Error_Out, ColumnIndexes=[0, 1])
    Data_Frame_Error_Percent_Full_Arr, n2 = FN.joinColumnsToArray(Array=Error_In, ArrtobeJoined=Error_Out, ColumnIndexes=[2, 3])
    n1.append("Absorption Difference(um" + "$^{}$".format(2) + ")")
    n1.append("Scattering Difference(um" + "$^{}$".format(2) + ")")
    n2.append("Absorption Difference(%)")
    n2.append("Scattering Difference(%)")
    Data_Frame_Error_Real_Full_Object = pd.DataFrame(Data_Frame_Error_Real_Full_Arr, columns=n1, dtype=float)
    Data_Frame_Error_Percent_Full_Object = pd.DataFrame(Data_Frame_Error_Percent_Full_Arr, columns=n2, dtype=float)
    FN.Fig_Plot_Save_Scatterplot_Matrix(Address=address_Graph_Real_Full_Scatter, Dataframe=Data_Frame_Error_Real_Full_Object)
    FN.Fig_Plot_Save_Scatterplot_Matrix(Address=address_Graph_Percent_Full_Scatter, Dataframe=Data_Frame_Error_Percent_Full_Object)

    FN.Fig_Plot_Save_Scatter_X_Linear_Y_Linear(Address=address_Graph_Real_Full, X_Array=ABS_Error_Real_Full, Y_array=SCA_Error_Real_Full, X_Label="Absorption Difference(um" + "$^{}$".format(2) + ")", Y_label="Scattering Difference(um" + "$^{}$".format(2) + ")", Plot_Title="RDG and T-Matrix Real Difference")
    FN.Fig_Plot_Save_Scatter_X_Linear_Y_Linear(Address=address_Graph_Percent_Full, X_Array=ABS_Error_Percent_Full, Y_array=SCA_Error_Percent_Full, tickLabelStyle='plain', X_Label="Absorption Difference (%)", Y_label="Scattering Difference (%)", Plot_Title="RDG and T-Matrix Percentage Difference")
    logging.info("Plotting error calculation finished.")
    ################################################################################################################
    ################################################################################################################
    ################################################################################################################
    ################################################################################################################
    ################################################################################################################
    logging.info("App finished.")
    A = 51
