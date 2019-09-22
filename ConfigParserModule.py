# Keyhan Babaee, https://github.com/KeyhanB
# V1.0
# July 2019
import configparser
import logging
import GeneralFunctions as GF


def readLogConfig():
    try:

        appDirectory = GF.GetRootDirectory()
        file_Address_Private_Config = GF.GetAddressTo(main=appDirectory, folderName="Config", fileName="Private", extension="cnf")
        file_Address_Public_Config = GF.GetAddressTo(main=appDirectory, folderName="Config", fileName="Public", extension="cnf")

        if GF.IsFileExist(file_Address_Private_Config):
            readLoggerConfig(file_Address_Private_Config)

        elif GF.IsFileExist(file_Address_Public_Config):
            readLoggerConfig(file_Address_Public_Config)

        else:
            raise Exception("Couldn't find config file")

    except Exception as e:
        logging.exception(e)
        raise


def readConfigToDict(SectionName, ConvertParseTo='string', hasComment=False):
    try:

        appDirectory = GF.GetRootDirectory()
        file_Address_Private_Config = GF.GetAddressTo(main=appDirectory, folderName="Config", fileName="Private", extension="cnf")
        file_Address_Public_Config = GF.GetAddressTo(main=appDirectory, folderName="Config", fileName="Public", extension="cnf")

        if GF.IsFileExist(file_Address_Private_Config):
            cnf = readDatabaseRaw(Address=file_Address_Private_Config, SectionName=SectionName, hasComment=hasComment)

        elif GF.IsFileExist(file_Address_Public_Config):
            cnf = readDatabaseRaw(Address=file_Address_Public_Config, SectionName=SectionName, hasComment=hasComment)

        else:
            raise Exception("Couldn't find config file")

        dict = {}
        if ConvertParseTo == 'string':
            for key, val in cnf.items():
                dict[key] = str(val)

        elif ConvertParseTo == 'float':
            for key, val in cnf.items():
                dict[key] = float(val)

        return dict

    except Exception as e:
        logging.exception(e)
        raise


def readDatabaseRaw(Address, SectionName, hasComment=False):
    try:

        if hasComment:
            config = configparser.ConfigParser(inline_comment_prefixes="#")
            config.optionxform = str
            config.read(Address)
        else:
            config = configparser.RawConfigParser()
            config.optionxform = str
            config.read(Address)

        sectionDictionary = {}
        for key, val in config.items(str(SectionName)):
            sectionDictionary[key] = val

        return sectionDictionary

    except Exception as e:
        logging.exception(e)
        raise


def readLoggerConfig(Address):
    try:
        config = configparser.RawConfigParser()
        config.optionxform = str
        config.read(Address)
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
