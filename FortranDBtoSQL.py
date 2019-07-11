import ConfigParserM as CP
import DBMethods as DB
from ConfigParserM import logging
import GeneralFunctions as GF

if __name__ == "__main__":
    CP.readLogConfig()
    DB_Info = CP.readConfigToDict(SectionName="DatabaseInfo")
    FF_Info = CP.readConfigToDict(SectionName="FilesFoldersInfo")
    logging.info("Import Started!")
    appDirectory = GF.getRootDirectory()
    ##############################################################################
    DB.createDB(INFO=DB_Info)
    # DB.dumpDB(INFO=DB_Info, FileAddress=GF.getAddressTo(appDirectory, FF_Info['FOLDER_NAME_DATABASE'], FileName=GF.getDateandTimeUTC(), Extension="sql.gz"))
    # DB.loadDB(INFO=DB_Info, FileAddress=GF.findLatestFile(GF.getFilesNameAddressinFolder(GF.getAddressTo(appDirectory, FF_Info['FOLDER_NAME_DATABASE']), Extension="sql.gz")))
    # DB.showAllTablesInDBSummary(DB_Info)
    # DB.reinitializeDB(DB_Info)
    ##############################################################################
    Address1p8 = GF.getAddressTo(Main=appDirectory, FolderName=FF_Info['FOLDER_NAME_BC_DATABASE'], FileName=FF_Info['FILE_NAME_BC_DATABASE_1p8'], Extension="DAT")
    Address2p3 = GF.getAddressTo(Main=appDirectory, FolderName=FF_Info['FOLDER_NAME_BC_DATABASE'], FileName=FF_Info['FILE_NAME_BC_DATABASE_2p3'], Extension="DAT")
    Address2p8 = GF.getAddressTo(Main=appDirectory, FolderName=FF_Info['FOLDER_NAME_BC_DATABASE'], FileName=FF_Info['FILE_NAME_BC_DATABASE_2p8'], Extension="DAT")

    Header = ['Df', 'R_RI', 'I_RI', 'MP', 'Np', 'Ext_CS_M', 'Sca_CS_M', 'Abs_CS_M']
    tableName = 'Raw_V1'
    DB.createTable(INFO=DB_Info, TableName=tableName, arrHeaderNamesInput=Header)

    ##############################################################################
    List2D1p8 = GF.DAT_To2DList(Address=Address1p8)
    List2D1p8 = GF.changeColumnOf2D(List=List2D1p8, ColumnNumber=0, ChangeTo="1.8")
    List2D1p8 = GF.selectColumnsList(ColumnIndex=[0, 1, 2, 3, 4, 6, 7, 8], List=List2D1p8)
    numberDivisions = 4
    List2D1p8_M = GF.divideArray(numberDivisions, List2D1p8)
    for i in range(numberDivisions):
        DB.insertArrayIntoTable(INFO=DB_Info, TableName=tableName, NameArray=Header, Array=List2D1p8_M[i])
    ##############################################################################
    List2D2p3 = GF.DAT_To2DList(Address=Address2p3)
    List2D2p3 = GF.changeColumnOf2D(List=List2D2p3, ColumnNumber=0, ChangeTo="2.3")
    List2D2p3 = GF.selectColumnsList(ColumnIndex=[0, 1, 2, 3, 4, 6, 7, 8], List=List2D2p3)
    numberDivisions = 4
    List2D2p3_M = GF.divideArray(numberDivisions, List2D2p3)
    for i in range(numberDivisions):
        DB.insertArrayIntoTable(INFO=DB_Info, TableName=tableName, NameArray=Header, Array=List2D2p3_M[i])
    ##############################################################################
    List2D2p8 = GF.DAT_To2DList(Address=Address2p8)
    List2D2p8 = GF.changeColumnOf2D(List=List2D2p8, ColumnNumber=0, ChangeTo="2.8")
    List2D2p8 = GF.selectColumnsList(ColumnIndex=[0, 1, 2, 3, 4, 6, 7, 8], List=List2D2p8)
    numberDivisions = 4
    List2D2p8_M = GF.divideArray(numberDivisions, List2D2p8)
    for i in range(numberDivisions):
        DB.insertArrayIntoTable(INFO=DB_Info, TableName=tableName, NameArray=Header, Array=List2D2p8_M[i])

    # CSV_Address = GF.getAddressTo(Main=appDirectory, FolderName=FF_Info['FOLDER_NAME_CSV'], FileName=GF.getDateandTimeUTC(), Extension="csv")
    # DB.dumpTableCSV(INFO=DB_Info, TableName=tableName, Address=CSV_Address)
    # DB.showAllTablesInDBSummary(DB_Info)
