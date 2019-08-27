from ConfigParserM import logging
import GeneralFunctions as GF
import numpy as np
from matplotlib import rcParams
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd


class GraphTools:
    def __init__(self, FolderInfo):
        try:
            self.__folderNameGraph = FolderInfo['FOLDER_NAME_GRAPH']
            self.__folderNameData = FolderInfo['FOLDER_NAME_DATA']
        except Exception as e:
            logging.exception(e)
            raise

    def RDG_TMatrixComparisonGraphs(self):
        try:
            MainDF = pd.read_csv(GF.getAddressTo(FolderName=self.__folderNameData, FileName='beacon', Extension='csv'))
            plotSeries = MainDF.loc[MainDF['AA_Plot'] == 1, 'AA_FileName']
            fullTable = MainDF.loc[MainDF['AA_Plot'] == 1, 'AA_Plot':]
            for index, item in plotSeries.iteritems():
                print(index)
                print(item)
            # A=plotSeries.to_dict()
            a = 3


        except Exception as e:
            logging.exception(e)
            raise

    def toSaveHistogram(self, Name, Array, Figure_DPI=400):
        try:
            Arr = []
            for i in range(len(Array)):
                Arr.append(float(Array[i]))
            n, bins, patches = plt.hist(Arr, 50, density=True, facecolor='b', alpha=0.75)
            plt.ylabel('Probability')
            plt.title('Histogram of ' + str(Name))
            plt.grid(True)
            # Address = GF.getAddressTo(Folder, None, Name, "jpg")
            #plt.savefig(Address, format='jpg', dpi=Figure_DPI, bbox_inches='tight')
            plt.clf()
            plt.close()

        except Exception as e:
            logging.exception(e)
            raise
