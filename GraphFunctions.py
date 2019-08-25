from ConfigParserM import logging
import AppFunctions as FN
import numpy as np
from matplotlib import rcParams
from matplotlib import pyplot as plt
import seaborn as sns


class GraphMethods:
    def __init__(self):
        pass

    def toSaveHistogram(Folder, Name, Array, Figure_DPI=400):
        try:
            Arr = []
            for i in range(len(Array)):
                Arr.append(float(Array[i]))
            n, bins, patches = plt.hist(Arr, 50, density=True, facecolor='b', alpha=0.75)
            plt.ylabel('Probability')
            plt.title('Histogram of ' + str(Name))
            plt.grid(True)
            Address = GF.getAddressTo(Folder, None, Name, "jpg")
            plt.savefig(Address, format='jpg', dpi=Figure_DPI, bbox_inches='tight')
            plt.clf()
            plt.close()

        except Exception as e:
            logging.exception(e)
            raise
