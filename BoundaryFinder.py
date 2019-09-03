from ConfigParserM import logging
import AppFunctions as FN


class BCDBBoundaryCheck:

    def __init__(self):
        try:
            self._arrDf_Bound = [1.801, 2.799]
            self._arrKf_Bound = [1.199, 1.201]
            self._arrRI_Real_Bound = [1.051, 1.999]
            self._arrRI_Imag_Bound = [0.001, 0.999]
            self._arrNp_Bound = [1, 3000]
            self._arrMonomerParameter_Bound = [0.051, 0.499]
        except Exception as e:
            logging.exception(e)
            raise

    def CheckInputDict(self, inputDict):
        try:
            outputDict = {}
            arrDf_Index = FN.getGoodIndexes(Array=inputDict['Df'], Bound=self._arrDf_Bound)
            arrKf_Index = FN.getGoodIndexes(Array=inputDict['Kf'], Bound=self._arrKf_Bound)
            arrRI_Real_Index = FN.getGoodIndexes(Array=inputDict['RI_Real'], Bound=self._arrRI_Real_Bound)
            arrRI_Imag_Index = FN.getGoodIndexes(Array=inputDict['RI_Imag'], Bound=self._arrRI_Imag_Bound)
            arrNp_Index = FN.getGoodIndexes(Array=inputDict['Np'], Bound=self._arrNp_Bound)
            arrMonomerParameter_Index = FN.getGoodIndexes(Array=inputDict['MP'], Bound=self._arrMonomerParameter_Bound)

            arrPossible_Indexes = FN.findCommonIndex(arrDf_Index, arrKf_Index, arrRI_Real_Index, arrRI_Imag_Index, arrNp_Index, arrMonomerParameter_Index)

            outputDict['Df'] = FN.getPossibleArray(Array=inputDict['Df'], Indexes=arrPossible_Indexes)
            outputDict['Kf'] = FN.getPossibleArray(Array=inputDict['Kf'], Indexes=arrPossible_Indexes)
            outputDict['RI_Real'] = FN.getPossibleArray(Array=inputDict['RI_Real'], Indexes=arrPossible_Indexes)
            outputDict['RI_Imag'] = FN.getPossibleArray(Array=inputDict['RI_Imag'], Indexes=arrPossible_Indexes)
            outputDict['Np'] = FN.getPossibleArray(Array=inputDict['Np'], Indexes=arrPossible_Indexes)
            outputDict['dp'] = FN.getPossibleArray(Array=inputDict['dp'], Indexes=arrPossible_Indexes)
            outputDict['wL'] = FN.getPossibleArray(Array=inputDict['wL'], Indexes=arrPossible_Indexes)
            outputDict['chance'] = FN.getPossibleArray(Array=inputDict['chance'], Indexes=arrPossible_Indexes)
            sum = 0
            for i in outputDict['chance']:
                sum += i
            outputDict['chanceSum'] = sum
            #### to do add more if needed
            logging.info("Boundary issued.")

            return outputDict

        except Exception as e:
            logging.exception(e)
            raise
