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
    AGG_Info = CP.readConfigToDict(SectionName="AggregateDetails", ConvertParseTo='float', hasComment=True)
    logging.info("config retrieved.")
    ##############################################
    G1 = GM.GraphTools(FolderInfo=FF_Info)
    G1.Experiment_RDG_TMatrixComparisonGraphs(AGG_Info=AGG_Info)
    ################################################################################################################
    logging.info("App finished.")
    A = 51
