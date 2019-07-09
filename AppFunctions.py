import DBMethods as DB
from ConfigParserM import logging


def isNewInputDB(INFO, TableName, InputDictionary):
    try:

        # DB.createTable(INFO=INFO, TableName=TableName, arrHeaderNamesInput=InputHeader, arrHeaderNamesOutput=OutputHeader)
        hash_Database = DB.getHashofRow(INFO=INFO, TableName=TableName, ValuesDictionary=InputDictionary)
        ID_Database, rowCount_Database = DB.checkHashValue(INFO=INFO, TableName=TableName, Hash=hash_Database)

        return ID_Database, rowCount_Database

    except Exception as e:
        logging.exception(e)
        raise
