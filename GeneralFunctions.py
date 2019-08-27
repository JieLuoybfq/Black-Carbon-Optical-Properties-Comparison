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

from ConfigParserM import logging


def writeFile(Address, String, Mode='w'):
    try:
        f = open(Address, Mode)
        f.write(String)
        f.close()

    except Exception as e:
        logging.exception(e)
        raise


def runFile(Address):
    try:
        # I changed the Fortran code so it can be used here.
        # p1 = subprocess.Popen([str(Address.resolve())])
        os.system(str(Address.resolve()))

    except Exception as e:
        logging.exception(e)
        raise


def readFile(FileAddress):
    try:

        A = open(FileAddress).read()
        return A

    except Exception as e:
        logging.exception(e)
        raise


def filterListFromNONE(List):
    try:

        return list(filter(None, List))

    except Exception as e:
        logging.exception(e)
        raise


def CSV_To2DList(Address):
    try:
        List2D = []
        with open(Address) as csvFile:
            readCSV = csv.reader(csvFile, delimiter=',')
            for row in readCSV:
                List2D.append(row)
        return List2D

    except Exception as e:
        logging.exception(e)
        raise


def DAT_To2DList(Address, getridofNone=True):
    try:
        List2D = []
        with open(Address) as datFile:
            readDAT = csv.reader(datFile, delimiter=' ')
            for row in readDAT:
                if getridofNone:
                    List2D.append(filterListFromNONE(row))
                else:
                    List2D.append(row)
        return List2D

    except Exception as e:
        logging.exception(e)
        raise


def findLatestFile(FolderAddressDictionary, DateTimeFormat='%Y-%m-%dT%H%M%S%z'):
    try:
        dT = {}

        for key in FolderAddressDictionary:
            if isDate(key):
                dT[key] = datetime.datetime.strptime(key, DateTimeFormat)

        C = max(dT, key=dT.get)

        return FolderAddressDictionary[C]

    except Exception as e:
        logging.exception(e)
        raise


def saveDictionaryAsConfig(Address, Dictionary, SectionName, mode='w'):
    try:
        parser = configparser.ConfigParser()
        parser.optionxform = str
        parser.add_section(str(SectionName))
        for key in Dictionary.keys():
            parser.set(str(SectionName), key, str(Dictionary[key]))

        with open(Address, mode) as f:
            parser.write(f)

    except Exception as e:
        logging.exception(e)
        raise


def reportDictionary(Dictionary):
    try:

        for key in Dictionary.keys():
            logging.info(f'{key}: {Dictionary[key]}')
            print(f'{key}: {Dictionary[key]}')

    except Exception as e:
        logging.exception(e)
        raise


def changeColumnOf2D(List, ColumnNumber, ChangeTo):
    try:

        rows = len(List)
        for i in range(rows):
            List[i][ColumnNumber] = ChangeTo

        return List

    except Exception as e:
        logging.exception(e)
        raise


def changeArrayType(List, ColumnNumber, ChangeTo):
    try:

        rows = len(List)
        for i in range(rows):
            List[i][ColumnNumber] = ChangeTo

        return List

    except Exception as e:
        logging.exception(e)
        raise


def selectColumnsList(ColumnIndex, List, Dimension=2):
    try:
        B = []
        if Dimension == 1:
            for i in ColumnIndex:
                B.append(List[i])

        elif Dimension == 2:
            for i in range(len(List)):
                A = []
                if len(ColumnIndex) == 1:
                    for j in ColumnIndex:
                        B.append(List[i][j])
                else:
                    for j in ColumnIndex:
                        A.append(List[i][j])
                    B.append(A)

        return B

    except Exception as e:
        logging.exception(e)
        raise


def divideArray(NumberofDivisions, List):
    try:
        B = []
        Counter = 0
        rows = len(List)
        Section = int(rows / NumberofDivisions)
        for i in range(NumberofDivisions - 1):
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


def getFilesNameAddressinFolder(FolderAddress, Extension=None):
    try:
        if Extension:
            S = str(FolderAddress.resolve()) + "\\*." + str(Extension)
        else:
            S = str(FolderAddress.resolve()) + "\\*.*"
        dict = {}
        for file in glob.glob(S):
            if Extension:
                Name = os.path.basename(file).replace("." + Extension, "")
            else:
                Name = os.path.basename(file)[:-3]
            dict[Name] = Path(file)

        return dict

    except Exception as e:
        logging.exception(e)
        raise


def reportVariable(Name, Value):
    try:

        logging.info(f'{Name}: {Value}')
        print(f'{Name}: {Value}')

    except Exception as e:
        logging.exception(e)
        raise


def isDate(string, fuzzy=False):
    try:

        try:
            parse(string, fuzzy=fuzzy)
            return True

        except ValueError:
            return False


    except Exception as e:
        logging.exception(e)
        raise


def getDateandTimeUTCNow():
    try:

        return str(datetime.datetime.utcnow().replace(tzinfo=datetime.timezone.utc).isoformat("T", "seconds")).replace(":", "")


    except Exception as e:
        logging.exception(e)
        raise


def getAddressTo(Main=None, FolderName=None, FileName=None, Extension=None):
    try:
        if Main is None:
            Main = getRootDirectory()
        if FolderName:
            Path1 = Path(Main) / Path(FolderName)
        else:
            Path1 = Path(Main)

        if not os.path.exists(Path1):
            os.makedirs(Path1)
        if FileName:
            if Extension:
                File_Address = Path1 / Path(FileName + "." + Extension)
            else:
                File_Address = Path1 / Path(FileName)
        else:
            File_Address = Path1
        return File_Address

    except Exception as e:
        logging.exception(e)
        raise


def getVarName(obj, namespace):
    try:
        A = [name for name in namespace if namespace[name] is obj]
        return A[0]

    except Exception as e:
        logging.exception(e)
        raise



def getRootDirectory():
    try:
        return Path(os.path.dirname(os.path.realpath('__file__')))

    except Exception as e:
        logging.exception(e)
        raise


def isFileExist(Address):
    try:
        return os.path.isfile(Address)

    except Exception as e:
        logging.exception(e)
        raise
