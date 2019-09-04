# Graphing using existing data
Version = 0.1
# Aug 2019
import ConfigParserM as CP
from ConfigParserM import logging
import GraphManagement as GM

if __name__ == "__main__":
    CP.readLogConfig()
    logging.info("App Started.")
    FF_Info = CP.readConfigToDict(SectionName="FilesFoldersInfo")
    logging.info("config retrieved.")
    ##############################################
    G1 = GM.GraphTools(FolderInfo=FF_Info)
    G1.RDG_TMatrixComparisonGraphs()
    ################################################################################################################
    logging.info("App finished.")
    A = 51
