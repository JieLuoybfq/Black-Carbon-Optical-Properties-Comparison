import ConfigParserM as CP
import DBMethods as DB
from ConfigParserM import logging
import GeneralFunctions as GF
from decimal import Decimal
import AppFunctions as FN

if __name__ == "__main__":
    CP.readLogConfig()
    DB_Info = CP.readConfigToDict(SectionName="DatabaseInfo")
    FF_Info = CP.readConfigToDict(SectionName="FilesFoldersInfo")
    ##############################################################################
    DB.createDB(INFO=DB_Info)
    # DB.showDBs(DB_Info)
    # Name='bc_db_liu_v1'
    DB.showAllTablesInDBSummary(DB_Info)
    # DB.dropDB(DB_Info,Name)
    # DB.reinitializeDB(DB_Info)
    # DB.showAllTablesInDBSummary(DB_Info,Name)
    ##############################################################################
    logging.info("Import Started!")
    appDirectory = GF.getRootDirectory()
    Address1p8 = GF.getAddressTo(Main=appDirectory, FolderName=FF_Info['FOLDER_NAME_BC_DATABASE'], FileName=FF_Info['FILE_NAME_BC_DATABASE_1p8'], Extension="DAT")
    Address2p3 = GF.getAddressTo(Main=appDirectory, FolderName=FF_Info['FOLDER_NAME_BC_DATABASE'], FileName=FF_Info['FILE_NAME_BC_DATABASE_2p3'], Extension="DAT")
    Address2p8 = GF.getAddressTo(Main=appDirectory, FolderName=FF_Info['FOLDER_NAME_BC_DATABASE'], FileName=FF_Info['FILE_NAME_BC_DATABASE_2p8'], Extension="DAT")

    Header = ['Df', 'R_RI', 'I_RI', 'MP', 'Np', 'Ext_CS_M', 'Sca_CS_M', 'Abs_CS_M']
    tableName = 'Raw_V1'
    DB.createTable(INFO=DB_Info, TableName=tableName, arrHeaderNamesInput=Header)

    ##############################################################################
    List2D1p8 = GF.DAT_To2DList(Address=Address1p8)
    List2D1p8 = GF.changeColumnOf2D(List=List2D1p8, ColumnNumber=0, ChangeTo="1.8")
    for i in range(len(List2D1p8)):
        inputDict = {
            'Df': Decimal(List2D1p8[i][0]),
            'R_RI': Decimal(List2D1p8[i][1]),
            'I_RI': Decimal(List2D1p8[i][2]),
            'MP': Decimal(List2D1p8[i][3]),
            'Np': Decimal(List2D1p8[i][4]),
            'Ext_CS_M': Decimal(List2D1p8[i][6]),
            'Sca_CS_M': Decimal(List2D1p8[i][7]),
            'Abs_CS_M': Decimal(List2D1p8[i][8])
        }
        DB.insertRowIntoTable(INFO=DB_Info, TableName=tableName, ValuesDictionary=inputDict)
    ##############################################################################
    List2D2p3 = GF.DAT_To2DList(Address=Address2p3)
    List2D2p3 = GF.changeColumnOf2D(List=List2D2p3, ColumnNumber=0, ChangeTo="2.3")
    for i in range(len(List2D2p3)):
        inputDict = {
            'Df': Decimal(List2D2p3[i][0]),
            'R_RI': Decimal(List2D2p3[i][1]),
            'I_RI': Decimal(List2D2p3[i][2]),
            'MP': Decimal(List2D2p3[i][3]),
            'Np': Decimal(List2D2p3[i][4]),
            'Ext_CS_M': Decimal(List2D2p3[i][6]),
            'Sca_CS_M': Decimal(List2D2p3[i][7]),
            'Abs_CS_M': Decimal(List2D2p3[i][8])
        }
        DB.insertRowIntoTable(INFO=DB_Info, TableName=tableName, ValuesDictionary=inputDict)
    ##############################################################################
    List2D2p8 = GF.DAT_To2DList(Address=Address2p8)
    List2D2p8 = GF.changeColumnOf2D(List=List2D2p8, ColumnNumber=0, ChangeTo="2.8")
    for i in range(len(List2D2p8)):
        inputDict = {
            'Df': Decimal(List2D2p8[i][0]),
            'R_RI': Decimal(List2D2p8[i][1]),
            'I_RI': Decimal(List2D2p8[i][2]),
            'MP': Decimal(List2D2p8[i][3]),
            'Np': Decimal(List2D2p8[i][4]),
            'Ext_CS_M': Decimal(List2D2p8[i][6]),
            'Sca_CS_M': Decimal(List2D2p8[i][7]),
            'Abs_CS_M': Decimal(List2D2p8[i][8])
        }
        DB.insertRowIntoTable(INFO=DB_Info, TableName=tableName, ValuesDictionary=inputDict)
    ##############################################################################

    DB.showAllTablesInDBSummary(DB_Info)

    B = 3
