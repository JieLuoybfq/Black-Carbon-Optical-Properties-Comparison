# Keyhan Babaee, https://github.com/KeyhanB
# V1.1
# July 2019
import mysql.connector
from ConfigParserM import logging
from mysql.connector import errorcode
import os
import pandas
import GeneralFunctions as GF


############################################    Server Level
def connectServer(INFO):
    try:
        mySQLDB = mysql.connector.connect(
            host=INFO["Host"],
            port=INFO["Port"],
            user=INFO["Username"],
            passwd=INFO["Password"])
        return mySQLDB

    except Exception as e:
        logging.exception(e)
        raise


########################################################################################

############################################    Database Level
def connectDB(INFO, DBName=None):
    try:
        if DBName:
            DBnametoconnect = DBName
        else:
            DBnametoconnect = INFO["DBName"]

        if DBnametoconnect:
            mySQLDB = mysql.connector.connect(
                host=INFO["Host"],
                port=INFO["Port"],
                user=INFO["Username"],
                passwd=INFO["Password"],
                database=DBnametoconnect)
            return mySQLDB
        else:
            raise Exception(f'DB Name is not provided={DBnametoconnect}')

    except Exception as e:
        logging.exception(e)
        raise


def createDB(INFO, DBName=None):
    try:
        mySQLDB = connectServer(INFO)
        if DBName:
            DBnametoconnect = DBName
        else:
            DBnametoconnect = INFO["DBName"]

        cursor = mySQLDB.cursor()
        cursor.execute("CREATE DATABASE IF NOT EXISTS " + str(DBnametoconnect))
        mySQLDB.commit()
        cursor.close()
        mySQLDB.close()

    except Exception as e:
        logging.exception(e)
        raise


def reinitializeDB(INFO, DBName=None):
    try:
        if DBName:
            DBnametoconnect = DBName
        else:
            DBnametoconnect = INFO["DBName"]

        dropDB(INFO, DBnametoconnect)
        createDB(INFO, DBnametoconnect)

    except Exception as e:
        logging.exception(e)
        raise


def dropDB(INFO, DBName=None):
    try:
        if DBName:
            DBnametoconnect = DBName
        else:
            DBnametoconnect = INFO["DBName"]

        mySQLDB = connectServer(INFO)
        cursor = mySQLDB.cursor()
        cursor.execute("DROP DATABASE IF EXISTS " + str(DBnametoconnect))
        mySQLDB.commit()
        cursor.close()
        mySQLDB.close()

    except Exception as e:
        logging.exception(e)
        raise


def showDBs(INFO):
    try:
        mySQLDB = connectServer(INFO)
        cursor = mySQLDB.cursor()
        cursor.execute("SHOW DATABASES")
        print("========================")
        for x in cursor:
            print(x)
        print("========================")

    except Exception as e:
        logging.exception(e)
        raise


def loadDB(INFO, FileAddress, DBName=None):
    try:
        if DBName:
            DBnametoconnect = DBName

        else:
            DBnametoconnect = INFO["DBName"]
        createDB(INFO, DBnametoconnect)
        os.system("gzip -dc %s | mysql --host=%s --port=%s --user=%s --password=%s %s" % (FileAddress.resolve(), INFO["Host"], INFO["Port"], INFO["Username"], INFO["Password"], DBnametoconnect))


    except Exception as e:
        logging.exception(e)
        raise


########################################################################################

############################################    Table Level
def createTable(INFO, TableName, arrHeaderNamesInput, DataType="DECIMAL(40,30)", HashType="SHA1", arrHeaderNamesOutput=None, OutputTableName="Out", withShadow=True, ShadowTableName="Shadow", DBName=None):
    try:
        mySQLDB = connectDB(INFO, DBName)
        cursor = mySQLDB.cursor()
        S_in = ""
        S_out = ""
        S_hash = ""
        if arrHeaderNamesInput:
            for name in arrHeaderNamesInput:
                S_in += str(name) + " " + str(DataType) + " NOT NULL , "
                S_hash += str(name) + ", "
            S_hash = S_hash[:-2]
            query = ("CREATE TABLE IF NOT EXISTS " + str(TableName)
                     + " (ID_" + str(TableName) + " INT NOT NULL AUTO_INCREMENT, "
                     + S_in
                     + "Hash varchar(255) AS (" + str(HashType) + "(CONCAT(" + S_hash + "))), "
                     + "created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP, "
                     + "updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP, "
                     + "PRIMARY KEY (ID_" + str(TableName) + "))"
                     + " ENGINE=InnoDB"
                     )
            cursor.execute(query)
            mySQLDB.commit()

        if arrHeaderNamesOutput:
            for name in arrHeaderNamesOutput:
                S_out += str(name) + " " + str(DataType) + " NOT NULL , "
            outputTableName = str(TableName) + "_" + str(OutputTableName)
            query = ("CREATE TABLE IF NOT EXISTS " + outputTableName
                     + " (ID_" + outputTableName + " INT NOT NULL AUTO_INCREMENT, "
                     + S_out
                     + "created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP, "
                     + "updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP, "
                     + " MLink INT NOT NULL REFERENCES " + str(TableName) + "(" + "ID_" + str(TableName) + "), "
                     + "PRIMARY KEY (ID_" + outputTableName + "))"
                     + " ENGINE=InnoDB"
                     )
            cursor.execute(query)
            mySQLDB.commit()

        if withShadow:
            shadowTableName = str(TableName) + "_" + str(ShadowTableName)
            query = ("CREATE TABLE IF NOT EXISTS " + shadowTableName
                     + " (ID_" + shadowTableName + " INT NOT NULL AUTO_INCREMENT, "
                     + S_in
                     + "Hash varchar(255) AS (SHA1(CONCAT(" + S_hash + "))), "
                     + "PRIMARY KEY (ID_" + shadowTableName + "))"
                     + " ENGINE=InnoDB"
                     )
            cursor.execute(query)
            mySQLDB.commit()

        cursor.close()
        mySQLDB.close()

    except Exception as e:
        logging.exception(e)
        raise


def showTablesInDB(INFO, DBName=None):
    try:
        mySQLDB = connectDB(INFO, DBName)
        cursor = mySQLDB.cursor()
        cursor.execute("SHOW TABLES")
        print("========================")
        for x in cursor:
            print(x)
        print("========================")

    except Exception as e:
        logging.exception(e)
        raise


def showAllTablesInDBSummary(INFO, DBName=None):
    try:
        mySQLDB = connectDB(INFO, DBName)
        cursor = mySQLDB.cursor()
        cursor.execute("SHOW TABLES")
        print("========================")
        for x in cursor:
            showTableSummary(INFO, x[0], DBName=DBName)
        print("========================")

    except Exception as e:
        logging.exception(e)
        raise


def dropTable(INFO, TableName, DBName=None):
    try:
        mySQLDB = connectDB(INFO, DBName)
        cursor = mySQLDB.cursor()
        cursor.execute("DROP TABLE IF EXISTS " + str(TableName))
        mySQLDB.commit()
        cursor.close()
        mySQLDB.close()

    except Exception as e:
        logging.exception(e)
        raise


def dropTableSet(INFO, TableName, DBName=None):
    try:
        mySQLDB = connectDB(INFO, DBName)
        cursor = mySQLDB.cursor()
        cursor.execute("DROP TABLE IF EXISTS " + str(TableName))
        mySQLDB.commit()
        cursor.execute("DROP TABLE IF EXISTS " + str(TableName) + "_Out")
        mySQLDB.commit()
        cursor.execute("DROP TABLE IF EXISTS " + str(TableName) + "_Shadow")
        mySQLDB.commit()
        cursor.close()
        mySQLDB.close()

    except Exception as e:
        logging.exception(e)
        raise


########################################################################################

############################################    insert in table

def insertRowIntoTable(INFO, TableName, ValuesDictionary, giveID=False, DBName=None):
    try:
        mySQLDB = connectDB(INFO, DBName)
        cursor = mySQLDB.cursor()
        S1 = " ("
        S2 = " ("
        for key in ValuesDictionary:
            S1 += str(key) + ", "
            S2 += "%(" + str(key) + ")s, "
        S1 = S1[:-2]
        S2 = S2[:-2]
        S1 += ") "
        S2 += ")"
        query = ("INSERT INTO " + str(TableName) + S1 + "VALUES" + S2)

        cursor.execute(query, ValuesDictionary)
        mySQLDB.commit()
        ID = cursor.lastrowid
        cursor.close()
        mySQLDB.close()
        if giveID:
            return ID

    except Exception as e:
        logging.exception(e)
        raise


def insertArrayIntoTable(INFO, TableName, NameArray, Array, DBName=None):
    try:
        mySQLDB = connectDB(INFO, DBName)
        cursor = mySQLDB.cursor()
        S1 = " ("
        S2 = " ("
        for key in NameArray:
            S1 += str(key) + ", "
            S2 += "%s, "
        S1 = S1[:-2]
        S2 = S2[:-2]
        S1 += ") "
        S2 += ")"
        query = ("INSERT INTO " + str(TableName) + S1 + "VALUES" + S2)

        cursor.executemany(query, Array)
        mySQLDB.commit()
        cursor.close()
        mySQLDB.close()


    except Exception as e:
        logging.exception(e)
        raise


########################################################################################

############################################    read from table
def readRowfromTablewithID(INFO, TableName, ID=1, DBName=None):
    try:
        mySQLDB = connectDB(INFO, DBName)
        cursor = mySQLDB.cursor()
        cursor.execute("SELECT * FROM " + str(TableName) + " WHERE " + "ID_" + str(TableName) + "=%s", (ID,))
        Row = cursor.fetchall()
        tableDict = {}
        i = 0
        for cd in cursor.description:
            tableDict[cd[0]] = Row[0][i]
            i += 1
        return tableDict

    except Exception as e:
        logging.exception(e)
        raise


def readAllRowsfromTable(INFO, TableName, DBName=None):
    try:
        mySQLDB = connectDB(INFO, DBName)
        cursor = mySQLDB.cursor()
        cursor.execute("SELECT * FROM " + str(TableName))
        Row = cursor.fetchall()
        columnName = []
        i = 0
        for cd in cursor.description:
            columnName.append(cd[0])

        return columnName, Row

    except Exception as e:
        logging.exception(e)
        raise


def readRowwithColumnNameandValue(INFO, TableName, ColumnName, Value, isRowCount=True, DBName=None):
    try:
        mySQLDB = connectDB(INFO, DBName)
        cursor = mySQLDB.cursor()
        cursor.execute("SELECT * FROM " + str(TableName) + " WHERE " + str(ColumnName) + "=%s", (Value,))

        Rows = cursor.fetchall()
        RowCount = cursor.rowcount
        tableDict = {}
        if RowCount:
            i = 0
            for cd in cursor.description:
                tableDict[cd[0]] = Rows[0][i]
                i += 1
        if isRowCount:
            return tableDict, RowCount
        else:
            return tableDict
    except Exception as e:
        logging.exception(e)
        raise


def readColwithColumnName(INFO, TableName, ColumnName, DBName=None):
    try:
        mySQLDB = connectDB(INFO, DBName)
        cursor = mySQLDB.cursor(dictionary=True)
        cursor.execute("SELECT " + str(ColumnName) + " FROM " + str(TableName))
        # sql = ('select field1, field2 from table')
        List = []
        for rows in cursor:
            List.append(rows[str(ColumnName)])
        return List

    except Exception as e:
        logging.exception(e)
        raise


########################################################################################


############################################    delete from table
def deleteRowfromTablewithID(INFO, TableName, ID=1, DBName=None):
    try:
        mySQLDB = connectDB(INFO, DBName)
        cursor = mySQLDB.cursor()
        cursor.execute("DELETE FROM " + str(TableName) + " WHERE " + "ID_" + str(TableName) + "=%s", (ID,))
        mySQLDB.commit()
        cursor.close()
        mySQLDB.close()

    except Exception as e:
        logging.exception(e)
        raise


def deleteAllRowsfromTable(INFO, TableName, DBName=None):
    try:
        mySQLDB = connectDB(INFO, DBName)
        cursor = mySQLDB.cursor()
        cursor.execute("DELETE FROM " + str(TableName))
        mySQLDB.commit()
        cursor.close()
        mySQLDB.close()

    except Exception as e:
        logging.exception(e)
        raise


########################################################################################


############################################    dump files
def dumpTableCSV(INFO, TableName, Address, DBName=None):
    try:
        mySQLDB = connectDB(INFO, DBName)
        query = "SELECT * FROM " + str(TableName)

        results = pandas.read_sql_query(query, mySQLDB)
        results.to_csv(Address, index=False)

    except Exception as e:
        logging.exception(e)
        raise


def dumpTableSetCSV(INFO, TableName, AddressMain, DBName=None):
    try:
        mySQLDB = connectDB(INFO, DBName)
        query = "SELECT * FROM " + str(TableName)

        A1 = GF.getAddressTo(AddressMain, "CSV", FileName=str(TableName) + "_" + GF.getDateandTimeUTC(), Extension="csv")
        results = pandas.read_sql_query(query, mySQLDB)
        results.to_csv(A1, index=False)

        query = "SELECT * FROM " + str(TableName) + "_Out"
        A1 = GF.getAddressTo(AddressMain, "CSV", FileName=str(TableName) + "_Out" + "_" + GF.getDateandTimeUTC(), Extension="csv")
        results = pandas.read_sql_query(query, mySQLDB)
        results.to_csv(A1, index=False)

        query = "SELECT * FROM " + str(TableName) + "_Shadow"
        A1 = GF.getAddressTo(AddressMain, "CSV", FileName=str(TableName) + "_Shadow" + "_" + GF.getDateandTimeUTC(), Extension="csv")
        results = pandas.read_sql_query(query, mySQLDB)
        results.to_csv(A1, index=False)

    except Exception as e:
        logging.exception(e)
        raise


def dumpDB(INFO, FileAddress, DBName=None):
    try:
        if DBName:
            DBnametoconnect = DBName
        else:
            DBnametoconnect = INFO["DBName"]

        os.system("mysqldump --host=%s --port=%s --user=%s --password=%s  -e --opt -c %s | gzip -c > %s" % (INFO["Host"], INFO["Port"], INFO["Username"], INFO["Password"], DBnametoconnect, FileAddress.resolve()))


    except Exception as e:
        logging.exception(e)
        raise


########################################################################################


############################################    Hash Section
def getHashofRow(INFO, TableName, ValuesDictionary, DBName=None):
    try:

        insertRowIntoTable(INFO, TableName + "_Shadow", ValuesDictionary, DBName=DBName)
        shadowTableRow = readAllRowsfromTable(INFO, TableName + "_Shadow", DBName=DBName)

        Hash = shadowTableRow["Hash"]
        deleteAllRowsfromTable(INFO, TableName + "_Shadow", DBName=DBName)
        return Hash

    except Exception as e:
        logging.exception(e)
        raise


def getHashofArray(INFO, TableName, Header, Array, DBName=None):
    try:

        insertArrayIntoTable(INFO, TableName + "_Shadow", Header, Array, DBName=DBName)
        HashCol = readColwithColumnName(INFO, TableName + "_Shadow", "Hash", DBName=DBName)
        deleteAllRowsfromTable(INFO, TableName + "_Shadow", DBName=DBName)
        return HashCol

    except Exception as e:
        logging.exception(e)
        raise


def checkHashValue(INFO, TableName, Hash, DBName=None):
    try:
        Row, Rowcount = readRowwithColumnNameandValue(INFO=INFO, TableName=TableName, ColumnName="Hash", Value=Hash, DBName=DBName)
        if Rowcount == 0:
            return -1, Rowcount
        elif Rowcount > 0:
            return Row["ID_" + str(TableName)], Rowcount


    except Exception as e:
        logging.exception(e)
        raise


def checkHashArray(INFO, TableName, HashArray, DBName=None):
    try:
        IndexFind = []
        HashMain = readColwithColumnName(INFO=INFO, TableName=TableName, ColumnName="Hash", DBName=DBName)
        for i in HashArray:
            if i in HashMain:
                IndexFind.append(i)
            else:
                IndexFind.append(-1)
        return IndexFind


    except Exception as e:
        logging.exception(e)
        raise


########################################################################################

############################################    reporting
def showTableSummary(INFO, TableName, Number_Entry=3, Length=20, DBName=None):
    try:
        mySQLDB = connectDB(INFO, DBName)
        cursor = mySQLDB.cursor()
        cursor.execute("SELECT * FROM " + str(TableName) + " ORDER BY " + "ID_" + str(TableName) + " DESC LIMIT " + str(Number_Entry))
        results = cursor.fetchmany(Number_Entry)

        widths = []
        columns = []
        tavnit = '|'
        separator = '+'
        # max_col_length = max(list(map(lambda x: len(str(x[4])), results)))
        max_col_length = Length
        for cd in cursor.description:
            widths.append(max(max_col_length, len(cd[0])))
            columns.append(cd[0])
        for w in widths:
            tavnit += " %-" + "%s.%ss |" % (w, w)
            separator += '-' * w + '--+'
        print(separator)
        print(tavnit % tuple(columns))
        print(separator)
        for row in results:
            print(tavnit % row)
        print(separator)


    except Exception as e:
        logging.exception(e)
        raise
########################################################################################
