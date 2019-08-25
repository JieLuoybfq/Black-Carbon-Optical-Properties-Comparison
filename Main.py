# Programmed by Keyhan Babaee Under Prof. Steven Rogak supervision, https://github.com/KeyhanB
# Using Fortran Code from: C. Liu, X. Xu, Y. Yin, M. Schnaiter, Y. L. Yung, 2018: Black carbon aggregate: A database for optical properties
Version = 0.2
# Aug 2019
import ConfigParserM as CP
from ConfigParserM import logging

import RDGTMatrixInputGenerator as InputG

if __name__ == "__main__":
    CP.readLogConfig()
    logging.info("App Started.")
    DB_Info = CP.readConfigToDict(SectionName="DatabaseInfo")
    FF_Info = CP.readConfigToDict(SectionName="FilesFoldersInfo")
    AGG_Info = CP.readConfigToDict(SectionName="AggregateDetails", ConvertParseTo='float', hasComment=True)
    logging.info("config retrieved.")
    ##############################################
    inputV1 = InputG.KeyhanV1(AGG_Info)
    inputV1.KeyhanV1Calc(DB_Info)
    ################################################################################################################
    logging.info("App finished.")
    A = 51
