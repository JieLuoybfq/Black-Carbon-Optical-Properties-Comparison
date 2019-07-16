'''
TMatrix_Table_Input_Headers1 = ['Df', 'kf', 'R_RI', 'I_RI', 'WaveL', 'dp', 'Np', 'Version']
TMatrix_Table_Output_Headers1 = ['TMatrix_ABS_CRS', 'TMatrix_SCA_CRS', 'MLink']
Array1 = [[2.3, 1.2, 1.6, 0.6, 860, 33, 50, 0.1], [2.4, 1.2, 1.6, 0.6, 860, 33, 50, 0.1], [2.6, 1.2, 1.6, 0.6, 860, 33, 50, 0.1]]
Array2 = [[0.056, 0.666, 20], [0.053, 0.666, 21], [0.051, 0.666, 22]]
DB.insertArrayIntoTable(INFO=DB_Info, TableName=TMatrix_Table_Name, NameArray=TMatrix_Table_Input_Headers1, Array=Array1)
DB.insertArrayIntoTable(INFO=DB_Info, TableName=TMatrix_Table_Name + "_out", NameArray=TMatrix_Table_Output_Headers1, Array=Array2)
DB.showAllTablesInDBSummary(DB_Info)
'''
'''
    TMatrix_Main_Input_Array.append([2.3, 1.2, 1.6, 0.6, 860, 33, 50, 0.1])
    TMatrix_Main_Input_Array.append([2.4, 1.2, 1.6, 0.6, 860, 33, 50, 0.1])
    TMatrix_Main_Input_Array.append([2.6, 1.2, 1.6, 0.6, 860, 33, 50, 0.1])
    TMatrix_Main_Input_Array.append([2.1, 1.2, 1.6, 0.6, 860, 33, 50, 0.1])
    '''
