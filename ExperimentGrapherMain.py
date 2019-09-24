# Graphing figures to compare experiment data with RDG-FA and T-matrix using simulated data
Version = 0.2
# Sep 2019
import ConfigReaderModule as CPM
from ConfigReaderModule import logging
import ExperimentGraphModule as EGM

if __name__ == "__main__":
    CPM.ReadLogConfig()
    logging.info("Graph App Started.")
    FF_Info = CPM.ReadConfigToDict(sectionName="FilesFoldersInfo")
    # AGG_Info = CPM.ReadConfigToDict(sectionName="AggregateDetails", convertParseTo='float', hasComment=True)
    logging.info("config retrieved.")
    ##############################################
    G1 = EGM.GraphTools(folderInfo=FF_Info)
    G1.Experiment_RDG_TMatrixComparisonGraphs()
    ################################################################################################################
    logging.info("Graph App finished.")
    A = 51
