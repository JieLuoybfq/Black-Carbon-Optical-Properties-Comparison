# Keyhan Babaee, https://github.com/KeyhanB
# V1.2
# July 2019
import configparser
import csv
import datetime
import glob
import os
from pathlib import Path
from dateutil.parser import parse
from ConfigParserModule import logging


####### List of Functions (updated Sep 21,2019)
# - WriteFile(address, string, mode='w')
# - RunFile(address)
# - ReadFile(fileAddress)
# - FilterListFromNONE(list)
# - CSV_To2DList(address)
# - DAT_To2DList(address, getRidOfNone=True)
# - FindLatestFile(folderAddressDictionary, dateTimeFormat='%Y-%m-%dT%H%M%S%z')
# - SaveDictionaryAsConfig(address, dictionary, sectionName, mode='w')
# - ReportDictionary(dictionary)
# - ChangeColumnOf2D(list, columnNumber, changeTo)
# - ChangeArrayType(list, columnNumber, changeTo)
# - SelectColumnsList(columnIndex, list, dimension=2)
# - DivideArray(numberOfDivisions, List)
# - GetFilesNameAddressInFolder(folderAddress, extension=None)
# - ReportVariable(name, value)
# - IsDate(string, fuzzy=False)
# - GetDateAndTimeUTCNow()
# - GetAddressTo(main=None, folderName=None, fileName=None, extension=None)
# - GetVarName(obj, namespace)
# - GetRootDirectory()
# - IsFileExist(address)
#######

def WriteFile(address, string, mode='w'):
    try:
        f = open(address, mode)
        f.write(string)
        f.close()

    except Exception as e:
        logging.exception(e)
        raise


def RunFile(address):
    try:
        # I changed the Fortran code so it can be used here.
        # p1 = subprocess.Popen([str(Address.resolve())])
        os.system(str(address.resolve()))

    except Exception as e:
        logging.exception(e)
        raise


def ReadFile(fileAddress):
    try:

        A = open(fileAddress).read()
        return A

    except Exception as e:
        logging.exception(e)
        raise


def FilterListFromNONE(list):
    try:

        return list(filter(None, list))

    except Exception as e:
        logging.exception(e)
        raise


def CSV_To2DList(address):
    try:
        List2D = []
        with open(address) as csvFile:
            readCSV = csv.reader(csvFile, delimiter=',')
            for row in readCSV:
                List2D.append(row)
        return List2D

    except Exception as e:
        logging.exception(e)
        raise


def DAT_To2DList(address, getRidOfNone=True):
    try:
        List2D = []
        with open(address) as datFile:
            readDAT = csv.reader(datFile, delimiter=' ')
            for row in readDAT:
                if getRidOfNone:
                    List2D.append(FilterListFromNONE(row))
                else:
                    List2D.append(row)
        return List2D

    except Exception as e:
        logging.exception(e)
        raise


def FindLatestFile(folderAddressDictionary, dateTimeFormat='%Y-%m-%dT%H%M%S%z'):
    try:
        dT = {}

        for key in folderAddressDictionary:
            if IsDate(key):
                dT[key] = datetime.datetime.strptime(key, dateTimeFormat)

        C = max(dT, key=dT.get)

        return folderAddressDictionary[C]

    except Exception as e:
        logging.exception(e)
        raise


def SaveDictionaryAsConfig(address, dictionary, sectionName, mode='w'):
    try:
        parser = configparser.ConfigParser()
        parser.optionxform = str
        parser.add_section(str(sectionName))
        for key in dictionary.keys():
            parser.set(str(sectionName), key, str(dictionary[key]))

        with open(address, mode) as f:
            parser.write(f)

    except Exception as e:
        logging.exception(e)
        raise


def ReportDictionary(dictionary):
    try:

        for key in dictionary.keys():
            logging.info(f'{key}: {dictionary[key]}')
            print(f'{key}: {dictionary[key]}')

    except Exception as e:
        logging.exception(e)
        raise


def ChangeColumnOf2D(list, columnNumber, changeTo):
    try:

        rows = len(list)
        for i in range(rows):
            list[i][columnNumber] = changeTo

        return list

    except Exception as e:
        logging.exception(e)
        raise


def ChangeArrayType(list, columnNumber, changeTo):
    try:

        rows = len(list)
        for i in range(rows):
            list[i][columnNumber] = changeTo

        return list

    except Exception as e:
        logging.exception(e)
        raise


def SelectColumnsList(columnIndex, list, dimension=2):
    try:
        B = []
        if dimension == 1:
            for i in columnIndex:
                B.append(list[i])

        elif dimension == 2:
            for i in range(len(list)):
                A = []
                if len(columnIndex) == 1:
                    for j in columnIndex:
                        B.append(list[i][j])
                else:
                    for j in columnIndex:
                        A.append(list[i][j])
                    B.append(A)

        return B

    except Exception as e:
        logging.exception(e)
        raise


def DivideArray(numberOfDivisions, List):
    try:
        B = []
        Counter = 0
        rows = len(List)
        Section = int(rows / numberOfDivisions)
        for i in range(numberOfDivisions - 1):
            C = []
            for j in range(Section):
                C.append(List[Counter])
                Counter += 1
            B.append(C)

        D = []
        for k in range(Counter, rows):
            D.append(List[Counter])
            Counter += 1
        B.append(D)

        return B


    except Exception as e:
        logging.exception(e)
        raise


def GetFilesNameAddressInFolder(folderAddress, extension=None):
    try:
        if extension:
            S = str(folderAddress.resolve()) + "\\*." + str(extension)
        else:
            S = str(folderAddress.resolve()) + "\\*.*"
        dict = {}
        for file in glob.glob(S):
            if extension:
                Name = os.path.basename(file).replace("." + extension, "")
            else:
                Name = os.path.basename(file)[:-3]
            dict[Name] = Path(file)

        return dict

    except Exception as e:
        logging.exception(e)
        raise


def ReportVariable(name, value):
    try:

        logging.info(f'{name}: {value}')
        print(f'{name}: {value}')

    except Exception as e:
        logging.exception(e)
        raise


def IsDate(string, fuzzy=False):
    try:

        try:
            parse(string, fuzzy=fuzzy)
            return True

        except ValueError:
            return False


    except Exception as e:
        logging.exception(e)
        raise


def GetDateAndTimeUTCNow():
    try:

        return str(datetime.datetime.utcnow().replace(tzinfo=datetime.timezone.utc).isoformat("T", "seconds")).replace(":", "")


    except Exception as e:
        logging.exception(e)
        raise


def GetAddressTo(main=None, folderName=None, fileName=None, extension=None):
    try:
        if main is None:
            main = GetRootDirectory()
        if folderName:
            Path1 = Path(main) / Path(folderName)
        else:
            Path1 = Path(main)

        if not os.path.exists(Path1):
            os.makedirs(Path1)
        if fileName:
            if extension:
                File_Address = Path1 / Path(fileName + "." + extension)
            else:
                File_Address = Path1 / Path(fileName)
        else:
            File_Address = Path1
        return File_Address

    except Exception as e:
        logging.exception(e)
        raise


def GetVarName(obj, namespace):
    try:
        A = [name for name in namespace if namespace[name] is obj]
        return A[0]

    except Exception as e:
        logging.exception(e)
        raise


def GetRootDirectory():
    try:
        return Path(os.path.dirname(os.path.realpath('__file__')))

    except Exception as e:
        logging.exception(e)
        raise


def IsFileExist(address):
    try:
        return os.path.isfile(address)

    except Exception as e:
        logging.exception(e)
        raise
