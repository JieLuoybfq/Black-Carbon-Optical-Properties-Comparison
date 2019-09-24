# Graphing  figures for T-matrix and RDG-FA comparison using simulated data
Version = 0.2
# Sep 2019
import ConfigReaderModule as CRM
from ConfigReaderModule import logging
import ComparisonGraphModule as CGM

if __name__ == "__main__":
    CRM.ReadLogConfig()
    logging.info("Graph App Started.")
    FF_Info = CRM.ReadConfigToDict(sectionName="FilesFoldersInfo")
    logging.info("Config retrieved.")
    ##############################################
    C1 = CGM.GraphTools(folderInfo=FF_Info)
    C1.RDG_TMatrixComparisonGraphs()
    ################################################################################################################
    logging.info("Graph App finished.")
    A = 51
