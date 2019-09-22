# Graphing  Figures for T-matrix and RDG-FA comparison using simulated data
Version = 0.2
# Sep 2019
import ConfigParserModule as CPM
from ConfigParserModule import logging
import ComparisonGraphModule as CGM

if __name__ == "__main__":
    CPM.readLogConfig()
    logging.info("Graph App Started.")
    FF_Info = CPM.readConfigToDict(SectionName="FilesFoldersInfo")
    logging.info("config retrieved.")
    ##############################################
    C1 = CGM.GraphTools(FolderInfo=FF_Info)
    C1.RDG_TMatrixComparisonGraphs()
    ################################################################################################################
    logging.info("App finished.")
    A = 51
