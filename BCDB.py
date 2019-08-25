class BCDBManagement:
    def __init__(self):
        self.RDG_Table_Name = "RDG_V1"
        self.RDG_Table_Input_Headers = ['Df', 'kf', 'R_RI', 'I_RI', 'WaveL', 'dp', 'Np', 'Version']
        self.RDG_Table_Output_Headers = ['RDG_ABS_CRS', 'RDG_SCA_CRS']
        DB.CreateTable(INFO=DB_Info, tableName=RDG_Table_Name, arrHeaderNamesInput=RDG_Table_Input_Headers, arrHeaderNamesOutput=RDG_Table_Output_Headers)
        ##############################################
        self.TMatrix_Table_Name = "TMatrix_V1"
        self.TMatrix_Table_Input_Headers = ['Df', 'kf', 'R_RI', 'I_RI', 'WaveL', 'dp', 'Np', 'Version']
        self.TMatrix_Table_Output_Headers = ['TMatrix_ABS_CRS', 'TMatrix_SCA_CRS']
        DB.CreateTable(INFO=DB_Info,
                       tableName=TMatrix_Table_Name,
                       arrHeaderNamesInput=TMatrix_Table_Input_Headers,
                       arrHeaderNamesOutput=TMatrix_Table_Output_Headers)
        ##############################################
        self.Error_Table_Name = "Error_V1"
        self.Error_Table_Input_Headers = ['Df', 'kf', 'R_RI', 'I_RI', 'WaveL', 'dp', 'Np', 'Version']
        self.Error_Table_Output_Headers = ['ERR_RDG_M_TMatrix_ABS_CRS', 'ERR_RDG_M_TMatrix_SCA_CRS', 'ERR_RDG_M_TMatrix_ABS_CRS_Percent',
                                           'ERR_RDG_M_TMatrix_SCA_CRS_Percent']
        DB.CreateTable(INFO=DB_Info, tableName=Error_Table_Name, arrHeaderNamesInput=Error_Table_Input_Headers, arrHeaderNamesOutput=Error_Table_Output_Headers)
