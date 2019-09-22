# Keyhan Babaee, https://github.com/KeyhanB
# V1.0
# July 2019
import configparser
import logging
import GeneralFunctions as GF


def ReadLogConfig():
    try:

        appDirectory = GF.GetRootDirectory()
        file_Address_Private_Config = GF.GetAddressTo(main=appDirectory, folderName="Config", fileName="Private", extension="cnf")
        file_Address_Public_Config = GF.GetAddressTo(main=appDirectory, folderName="Config", fileName="Public", extension="cnf")

        if GF.IsFileExist(file_Address_Private_Config):
            ReadLoggerConfig(file_Address_Private_Config)

        elif GF.IsFileExist(file_Address_Public_Config):
            ReadLoggerConfig(file_Address_Public_Config)

        else:
            raise Exception("Couldn't find config file")

    except Exception as e:
        logging.exception(e)
        raise


def ReadConfigToDict(sectionName, convertParseTo='string', hasComment=False):
    try:

        appDirectory = GF.GetRootDirectory()
        file_Address_Private_Config = GF.GetAddressTo(main=appDirectory, folderName="Config", fileName="Private", extension="cnf")
        file_Address_Public_Config = GF.GetAddressTo(main=appDirectory, folderName="Config", fileName="Public", extension="cnf")

        if GF.IsFileExist(file_Address_Private_Config):
            cnf = ReadDatabaseRaw(address=file_Address_Private_Config, sectionName=sectionName, hasComment=hasComment)

        elif GF.IsFileExist(file_Address_Public_Config):
            cnf = ReadDatabaseRaw(address=file_Address_Public_Config, sectionName=sectionName, hasComment=hasComment)

        else:
            raise Exception("Couldn't find config file")

        dict = {}
        if convertParseTo == 'string':
            for key, val in cnf.items():
                dict[key] = str(val)

        elif convertParseTo == 'float':
            for key, val in cnf.items():
                dict[key] = float(val)

        return dict

    except Exception as e:
        logging.exception(e)
        raise


def ReadDatabaseRaw(address, sectionName, hasComment=False):
    try:

        if hasComment:
            config = configparser.ConfigParser(inline_comment_prefixes="#")
            config.optionxform = str
            config.read(address)
        else:
            config = configparser.RawConfigParser()
            config.optionxform = str
            config.read(address)

        sectionDictionary = {}
        for key, val in config.items(str(sectionName)):
            sectionDictionary[key] = val

        return sectionDictionary

    except Exception as e:
        logging.exception(e)
        raise


def ReadLoggerConfig(address):
    try:
        config = configparser.RawConfigParser()
        config.optionxform = str
        config.read(address)
        logFormat = config.get('LoggingInfo', 'format')
        logDate = config.get('LoggingInfo', 'datefmt')
        logFileName = config.get('LoggingInfo', 'filename')

        LEVELS = {
            'debug': logging.DEBUG,
            'info': logging.INFO,
            'warning': logging.WARNING,
            'error': logging.ERROR,
            'critical': logging.CRITICAL}

        logLevel = config.get('LoggingInfo', 'level')
        logFileMode = config.get('LoggingInfo', 'filemode')
        logging.basicConfig(format=logFormat, datefmt=logDate, filename=logFileName, level=LEVELS[logLevel], filemode=logFileMode)

    except Exception as e:
        logging.exception(e)
        raise
