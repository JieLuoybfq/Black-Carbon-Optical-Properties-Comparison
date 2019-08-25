# Keyhan Babaee, https://github.com/KeyhanB
# V1.2
# Aug 2019
import mysql.connector
from mysql.connector import Error as SQLError
from ConfigParserM import logging
import os
import pandas as pd


class MySQLManagement:

    def __init__(self, infoDic):
        try:
            self.__host = infoDic["Host"]
            self.__port = infoDic["Port"]
            self.__user = infoDic["Username"]
            self.__password = infoDic["Password"]
            self.__DBNameToConnect = infoDic["DBName"]

        except Exception as e:
            logging.exception(e)
            raise

    def CheckDBName(self, DBName):
        try:
            if DBName:
                DBnameToConnect = DBName
            elif self.__DBNameToConnect:
                DBnameToConnect = self.__DBNameToConnect
            else:
                raise Exception("DB Name is not provided")
            return DBnameToConnect

        except Exception as e:
            logging.exception(e)
            raise

    def ConnectServer(self):
        try:
            serverConnection = mysql.connector.connect(
                host=self.__host,
                port=self.__port,
                user=self.__user,
                passwd=self.__password)
            return serverConnection

        except SQLError as error:
            serverConnection.rollback()
            logging.exception(f"Connection to the server is failed: {error}")
            raise

        except Exception as e:
            logging.exception(e)
            raise

    ########################################################################################
    ############################################    Database Level
    def ConnectDB(self, DBName=None):
        try:
            DBnameToConnect = self.CheckDBName(DBName=DBName)
            DBConnection = mysql.connector.connect(
                host=self.__host,
                port=self.__port,
                user=self.__user,
                passwd=self.__password,
                database=DBnameToConnect)
            return DBConnection

        except SQLError as error:
            DBConnection.rollback()
            logging.exception(f"Connection to the server is failed: {error}")
            raise

        except Exception as e:
            logging.exception(e)
            raise

    def CreateDB(self, DBName=None):
        try:
            # don't expose this to public: SQL injection danger
            serverConnection = self.ConnectServer()
            DBnameToConnect = self.CheckDBName(DBName=DBName)
            cursor = serverConnection.cursor()
            query = f"CREATE DATABASE IF NOT EXISTS {DBnameToConnect}"
            cursor.execute(query)
            serverConnection.commit()

        except SQLError as error:
            serverConnection.rollback()
            logging.exception(f"Connection to the server failed: {error}")
            raise

        except Exception as e:
            logging.exception(e)
            raise

        finally:
            if serverConnection.is_connected():
                cursor.close()
                serverConnection.close()
                logging.info("MySQL connection is closed")

    def ReInitializeDB(self, DBName=None):
        try:
            DBnameToConnect = self.CheckDBName(DBName=DBName)
            self.DropDB(DBName=DBnameToConnect)
            self.CreateDB(DBName=DBnameToConnect)

        except Exception as e:
            logging.exception(e)
            raise

    def DropDB(self, DBName=None):
        try:
            # don't expose this to public: SQL injection danger
            DBnameToConnect = self.CheckDBName(DBName=DBName)
            serverConnection = self.ConnectServer()
            cursor = serverConnection.cursor()
            cursor.execute(f"DROP DATABASE IF EXISTS {DBnameToConnect}")
            serverConnection.commit()

        except SQLError as error:
            serverConnection.rollback()
            logging.exception(f"Connection to the server failed: {error}")

        except Exception as e:
            logging.exception(e)
            raise

        finally:
            if serverConnection.is_connected():
                cursor.close()
                serverConnection.close()
                logging.info("MySQL connection is closed")

    def ShowDBs(self):
        try:
            serverConnection = self.ConnectServer()
            cursor = serverConnection.cursor()
            cursor.execute("SHOW DATABASES")
            print("======================== Databases ========================")
            for x in cursor:
                print(x[0])
            print("================================================")

        except SQLError as error:
            serverConnection.rollback()
            logging.exception(f"Connection to the server failed: {error}")

        except Exception as e:
            logging.exception(e)
            raise

        finally:
            if serverConnection.is_connected():
                cursor.close()
                serverConnection.close()
                logging.info("MySQL connection is closed")

    def LoadDB(self, fileAddress, DBName=None):
        try:
            DBnameToConnect = self.CheckDBName(DBName=DBName)
            self.CreateDB(DBnameToConnect)
            os.system("gzip -dc %s | mysql --host=%s --port=%s --user=%s --password=%s %s" % (
                fileAddress.resolve(), self.__host, self.__port, self.__user, self.__password, DBnameToConnect))

        except Exception as e:
            logging.exception(e)
            raise

    ########################################################################################
    ############################################    Table Level
    def CreateTable(self, tableName, dictInputHeader, withHash=True, HashType="SHA1", withUnique=True, dictOutputHeader=None,
                    outputTableNameExt="Out", withCreated_UpdatedAt=True, withShadowTable=True,
                    shadowTableNameExt="Shadow", DBEngine="InnoDB", DBName=None):
        try:
            DBConnection = self.ConnectDB(DBName=DBName)
            cursor = DBConnection.cursor()

            if dictInputHeader:

                S_main = ""
                S_in = ""
                S_hash = ""
                S_unique = ""
                S_createdAt = ""
                S_updatedAt = ""
                S_primary = ""
                S_foreign = ""
                withPrimaryKey = True
                withForeignKey = False

                for key in dictInputHeader:

                    item = dictInputHeader[key]

                    S_in += f"{key} {item['dataType']}"
                    if item['notNull'] == 1:
                        S_in += " NOT NULL , "
                    elif item['notNull'] == 0:
                        S_in += " , "

                    if withHash:
                        if item['inHash'] == 1:
                            S_hash += f"{key} , "
                    if withUnique:
                        if item['inUnique'] == 1:
                            S_unique += f"{key} , "
                S_main += S_in

                if withHash:
                    S_hash = f"Hash varchar(255) AS ({HashType} (CONCAT({S_hash[:-2]}))),"
                    S_main += S_hash

                if withCreated_UpdatedAt:
                    S_createdAt = "created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP, "
                    S_updatedAt = "updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP, "
                    S_main += S_createdAt
                    S_main += S_updatedAt

                if withUnique:
                    if withHash:
                        S_unique += "Hash, "
                    S_unique = f"UNIQUE({S_unique[:-2]}), "
                    S_main += S_unique

                if withPrimaryKey:
                    S_primary += f"PRIMARY KEY (ID_{tableName}), "
                    S_main += S_primary

                if withForeignKey:
                    S_foreign += f"MLink INT NOT NULL REFERENCES {tableName}(ID_{tableName}), "
                    S_main += S_foreign

                query = f"CREATE TABLE IF NOT EXISTS {tableName} ( ID_{tableName} INT NOT NULL AUTO_INCREMENT , {S_main[:-2]}) ENGINE={DBEngine}"
                cursor.execute(query)
                DBConnection.commit()

            if dictOutputHeader:
                S_main = ""
                S_in = ""
                S_hash = ""
                S_unique = ""
                S_createdAt = ""
                S_updatedAt = ""
                S_primary = ""
                S_foreign = ""
                withPrimaryKeyOut = True
                withForeignKeyOut = True
                withUniqueOut = False
                withHashOut = False

                for key in dictOutputHeader:

                    item = dictOutputHeader[key]

                    S_in += f"{key} {item['dataType']}"
                    if item['notNull'] == 1:
                        S_in += " NOT NULL , "
                    elif item['notNull'] == 0:
                        S_in += " , "

                    if withHashOut:
                        if item['inHash'] == 1:
                            S_hash += f"{key} , "
                    if withUniqueOut:
                        if item['inUnique'] == 1:
                            S_unique += f"{key} , "
                S_main += S_in

                if withHashOut:
                    S_hash = f"Hash varchar(255) AS ({HashType} (CONCAT({S_hash[:-2]}))),"
                    S_main += S_hash

                if withCreated_UpdatedAt:
                    S_createdAt = "created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP, "
                    S_updatedAt = "updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP, "
                    S_main += S_createdAt
                    S_main += S_updatedAt

                if withUniqueOut:
                    if withHashOut:
                        S_unique += "Hash, "
                    S_unique = f"UNIQUE({S_unique[:-2]}), "
                    S_main += S_unique

                if withPrimaryKeyOut:
                    S_primary += f"PRIMARY KEY (ID_{tableName}_{outputTableNameExt}), "
                    S_main += S_primary

                if withForeignKeyOut:
                    S_foreign += f"MLink INT NOT NULL REFERENCES {tableName}(ID_{tableName}), "
                    S_main += S_foreign

                query = f"CREATE TABLE IF NOT EXISTS {tableName}_{outputTableNameExt} ( ID_{tableName}_{outputTableNameExt} INT NOT NULL AUTO_INCREMENT , {S_main[:-2]}) ENGINE={DBEngine}"
                cursor.execute(query)
                DBConnection.commit()

            if withShadowTable:
                S_main = ""
                S_in = ""
                S_hash = ""
                S_unique = ""
                S_createdAt = ""
                S_updatedAt = ""
                S_primary = ""
                S_foreign = ""
                withPrimaryKeyShadow = True
                withForeignKeyShadow = False
                withUniqueShadow = True
                withHashShadow = True
                withCreated_UpdatedAtShadow = False
                for key in dictInputHeader:

                    item = dictInputHeader[key]

                    S_in += f"{key} {item['dataType']}"
                    if item['notNull'] == 1:
                        S_in += " NOT NULL , "
                    elif item['notNull'] == 0:
                        S_in += " , "

                    if withHashShadow:
                        if item['inHash'] == 1:
                            S_hash += f"{key} , "
                    if withUniqueShadow:
                        if item['inUnique'] == 1:
                            S_unique += f"{key} , "
                S_main += S_in

                if withHashShadow:
                    S_hash = f"Hash varchar(255) AS ({HashType} (CONCAT({S_hash[:-2]}))),"
                    S_main += S_hash

                if withCreated_UpdatedAtShadow:
                    S_createdAt = "created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP, "
                    S_updatedAt = "updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP, "
                    S_main += S_createdAt
                    S_main += S_updatedAt

                if withUniqueShadow:
                    if withHashShadow:
                        S_unique += "Hash, "
                    S_unique = f"UNIQUE({S_unique[:-2]}), "
                    S_main += S_unique

                if withPrimaryKeyShadow:
                    S_primary += f"PRIMARY KEY (ID_{tableName}_{shadowTableNameExt}), "
                    S_main += S_primary

                if withForeignKeyShadow:
                    S_foreign += f"MLink INT NOT NULL REFERENCES {tableName}(ID_{tableName}), "
                    S_main += S_foreign

                query = f"CREATE TABLE IF NOT EXISTS {tableName}_{shadowTableNameExt} ( ID_{tableName}_{shadowTableNameExt} INT NOT NULL AUTO_INCREMENT , {S_main[:-2]}) ENGINE={DBEngine}"
                cursor.execute(query)
                DBConnection.commit()

        except SQLError as error:
            DBConnection.rollback()
            logging.exception(f"Connection to the database failed: {error}")

        except Exception as e:
            logging.exception(e)
            raise

        finally:
            if DBConnection.is_connected():
                cursor.close()
                DBConnection.close()
                logging.info("MySQL connection is closed")

    def ShowTablesInDB(self, DBName=None):
        try:
            DBConnection = self.ConnectDB(DBName=DBName)
            cursor = DBConnection.cursor()
            cursor.execute("SHOW TABLES")
            print(f"======================== {cursor.column_names[0]} ========================")
            for x in cursor:
                print(x[0])
            print("================================================")

        except Exception as e:
            logging.exception(e)
            raise

    def ShowSummaryAllTablesInDB(self, DBName=None):
        try:
            DBConnection = self.ConnectDB(DBName=DBName)
            cursor = DBConnection.cursor()
            cursor.execute("SHOW TABLES")
            print(f"======================== {cursor.column_names[0]} ========================")
            for x in cursor:
                self.ShowTableSummary(x[0], DBName=DBName)
            print("================================================")

        except Exception as e:
            logging.exception(e)
            raise

    def DropTable(self, TableName, DBName=None):
        try:
            DBConnection = self.ConnectDB(DBName=DBName)
            cursor = DBConnection.cursor()
            cursor.execute(f"DROP TABLE IF EXISTS {TableName}")
            DBConnection.commit()

        except SQLError as error:
            DBConnection.rollback()
            logging.exception(f"Connection to the database failed: {error}")

        except Exception as e:
            logging.exception(e)
            raise

        finally:
            if DBConnection.is_connected():
                cursor.close()
                DBConnection.close()
                logging.info("MySQL connection is closed")

    def DropTableSet(self, TableName, DBName=None):
        try:

            self.DropTable(TableName=TableName, DBName=DBName)
            self.DropTable(TableName=TableName + "_Out", DBName=DBName)
            self.DropTable(TableName=TableName + "_Shadow", DBName=DBName)

        except Exception as e:
            logging.exception(e)
            raise

    ########################################################################################
    ############################################    insert in table

    def InsertRowIntoTable(self, TableName, ValuesDictionary, giveID=False, DBName=None):
        try:
            DBConnection = self.ConnectDB(DBName=DBName)
            cursor = DBConnection.cursor()
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
            DBConnection.commit()
            ID = cursor.lastrowid

            if giveID:
                return ID

        except SQLError as error:
            DBConnection.rollback()
            logging.exception(f"Connection to the database failed: {error}")


        except Exception as e:
            logging.exception(e)
            raise


        finally:
            if DBConnection.is_connected():
                cursor.close()
                DBConnection.close()
                logging.info("MySQL connection is closed")

    def InsertArrayIntoTable(self, TableName, NameArray, Array, giveID=False, DBName=None):
        try:
            DBConnection = self.ConnectDB(DBName=DBName)
            cursor = DBConnection.cursor()
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
            ID = cursor.lastrowid
            DBConnection.commit()

            if giveID:
                return ID

        except SQLError as error:
            DBConnection.rollback()
            logging.exception(f"Connection to the database failed: {error}")

        except Exception as e:
            logging.exception(e)
            raise

        finally:
            if DBConnection.is_connected():
                cursor.close()
                DBConnection.close()
                logging.info("MySQL connection is closed")

    ########################################################################################

    ############################################    read from table
    def ReadRowfromTablewithID(self, TableName, ID=1, DBName=None):
        try:
            DBConnection = self.ConnectDB(DBName=DBName)
            cursor = DBConnection.cursor()
            cursor.execute("SELECT * FROM " + str(TableName) + " WHERE " + "ID_" + str(TableName) + "=%s", (ID,))
            Row = cursor.fetchall()
            tableDict = {}
            i = 0
            for cd in cursor.description:
                tableDict[cd[0]] = Row[0][i]
                i += 1
            return tableDict


        except SQLError as error:
            DBConnection.rollback()
            logging.exception(f"Connection to the database failed: {error}")


        except Exception as e:
            logging.exception(e)
            raise


        finally:
            if DBConnection.is_connected():
                cursor.close()
                DBConnection.close()
                logging.info("MySQL connection is closed")

    def ReadAllRowsfromTable(self, TableName, DBName=None):
        try:
            DBConnection = self.ConnectDB(DBName=DBName)
            cursor = DBConnection.cursor()
            cursor.execute("SELECT * FROM " + str(TableName))
            Row = cursor.fetchall()
            columnName = []
            i = 0
            for cd in cursor.description:
                columnName.append(cd[0])

            return columnName, Row


        except SQLError as error:
            DBConnection.rollback()
            logging.exception(f"Connection to the database failed: {error}")


        except Exception as e:
            logging.exception(e)
            raise


        finally:

            if DBConnection.is_connected():
                cursor.close()
                DBConnection.close()
                logging.info("MySQL connection is closed")

    def ReadRowwithColumnNameandValue(self, TableName, ColumnName, Value, isRowCount=True, DBName=None):
        try:
            DBConnection = self.ConnectDB(DBName=DBName)
            cursor = DBConnection.cursor()
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

        except SQLError as error:
            DBConnection.rollback()
            logging.exception(f"Connection to the database failed: {error}")

        except Exception as e:
            logging.exception(e)
            raise

        finally:
            if DBConnection.is_connected():
                cursor.close()
                DBConnection.close()
                logging.info("MySQL connection is closed")

    def ReadColwithColumnName(self, TableName, ColumnName, DBName=None):
        try:
            def sortList(Target, Index):
                zipped_pairs = zip(Index, Target)
                z = [x for _, x in sorted(zipped_pairs)]
                return z

            DBConnection = self.ConnectDB(DBName=DBName)
            cursor = DBConnection.cursor(dictionary=True)
            cursor.execute("SELECT " + str(ColumnName) + ", " + "ID_" + str(TableName) + " FROM " + str(TableName))
            # sql = ('select field1, field2 from table')
            List1 = []
            ID = []
            for rows in cursor:
                List1.append(rows[str(ColumnName)])
                ID.append(rows["ID_" + str(TableName)])
            List = sortList(List1, ID)

            return List

        except SQLError as error:
            DBConnection.rollback()
            logging.exception(f"Connection to the database failed: {error}")

        except Exception as e:
            logging.exception(e)
            raise

        finally:
            if DBConnection.is_connected():
                cursor.close()
                DBConnection.close()
                logging.info("MySQL connection is closed")

    def GetOutputRowByHash(self, TableName, Hash):
        try:
            Main = self.ReadRowwithColumnNameandValue(TableName=TableName, ColumnName="Hash", Value=Hash, isRowCount=False)
            Output = self.ReadRowwithColumnNameandValue(TableName=TableName + "_out", ColumnName="MLink", Value=Main['ID_' + str(TableName)],
                                                        isRowCount=False)
            return Output

        except Exception as e:
            logging.exception(e)
            raise

    ########################################################################################

    ############################################    delete from table
    def DeleteRowfromTablewithID(self, TableName, ID=1, DBName=None):
        try:
            DBConnection = self.ConnectDB(DBName=DBName)
            cursor = DBConnection.cursor()
            cursor.execute("DELETE FROM " + str(TableName) + " WHERE " + "ID_" + str(TableName) + "=%s", (ID,))
            DBConnection.commit()



        except SQLError as error:
            DBConnection.rollback()
            logging.exception(f"Connection to the database failed: {error}")


        except Exception as e:
            logging.exception(e)
            raise


        finally:
            if DBConnection.is_connected():
                cursor.close()
                DBConnection.close()
                logging.info("MySQL connection is closed")

    def DeleteAllRowsfromTable(self, TableName, DBName=None):
        try:
            DBConnection = self.ConnectDB(DBName=DBName)
            cursor = DBConnection.cursor()
            cursor.execute("DELETE FROM " + str(TableName))
            DBConnection.commit()


        except SQLError as error:
            DBConnection.rollback()
            logging.exception(f"Connection to the database failed: {error}")


        except Exception as e:
            logging.exception(e)
            raise


        finally:
            if DBConnection.is_connected():
                cursor.close()
                DBConnection.close()
                logging.info("MySQL connection is closed")

    ########################################################################################
    ############################################    dump files
    def DumpTableCSV(self, TableName, Address, DBName=None):
        try:
            DBConnection = self.ConnectDB(DBName=DBName)
            query = "SELECT * FROM " + str(TableName)

            results = pd.read_sql_query(query, DBConnection)
            results.to_csv(Address, index=False)

        except Exception as e:
            logging.exception(e)
            raise

    def DumpTableSetCSV(self, TableName, AddressMain, DBName=None):
        try:
            self.DumpTableCSV(TableName=TableName, Address=AddressMain, DBName=DBName)
            self.DumpTableCSV(TableName=TableName + "_out", Address=AddressMain + "_out.csv", DBName=DBName)
            self.DumpTableCSV(TableName=TableName + "_shadow", Address=AddressMain + "_shadow.csv", DBName=DBName)

        except Exception as e:
            logging.exception(e)
            raise

    def DumpDB(self, FileAddress, DBName=None):
        try:
            DBnameToConnect = self.CheckDBName(DBName=DBName)

            os.system("mysqldump --host=%s --port=%s --user=%s --password=%s  -e --opt -c %s | gzip -c > %s" % (
                self.__host, self.__port, self.__user, self.__password, DBnameToConnect, FileAddress.resolve()))


        except Exception as e:
            logging.exception(e)
            raise

    ########################################################################################
    ############################################    Hash Section
    def GetHashofRow(self, TableName, ValuesDictionary, DBName=None):
        try:

            self.InsertRowIntoTable(TableName + "_Shadow", ValuesDictionary, DBName=DBName)
            shadowTableRow = self.ReadAllRowsfromTable(TableName + "_Shadow", DBName=DBName)
            Hash = shadowTableRow["Hash"]
            self.DeleteAllRowsfromTable(TableName + "_Shadow", DBName=DBName)
            return Hash

        except Exception as e:
            logging.exception(e)
            raise

    def GetHashofArray(self, TableName, Header, Array, DBName=None):
        try:

            self.InsertArrayIntoTable(TableName + "_Shadow", Header, Array, DBName=DBName)
            HashCol = self.ReadColwithColumnName(TableName + "_Shadow", "Hash", DBName=DBName)
            self.DeleteAllRowsfromTable(TableName + "_Shadow", DBName=DBName)
            return HashCol

        except Exception as e:
            logging.exception(e)
            raise

    def CheckHashValue(self, TableName, Hash, DBName=None):
        try:
            Row, Rowcount = self.ReadRowwithColumnNameandValue(TableName=TableName, ColumnName="Hash", Value=Hash, DBName=DBName)
            if Rowcount == 0:
                return -1, Rowcount
            elif Rowcount > 0:
                return Row["ID_" + str(TableName)], Rowcount


        except Exception as e:
            logging.exception(e)
            raise

    def CheckHashArray(self, TableName, HashArray, DBName=None):
        try:
            IndexFind = []
            HashMain = self.ReadColwithColumnName(TableName=TableName, ColumnName="Hash", DBName=DBName)
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
    def ShowTableSummary(self, TableName, Number_Entry=3, Length=20, DBName=None):
        try:
            DBConnection = self.ConnectDB(DBName=DBName)
            cursor = DBConnection.cursor()
            cursor.execute(f"SELECT * FROM {TableName} ORDER BY ID_{TableName} DESC LIMIT {Number_Entry}")
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
