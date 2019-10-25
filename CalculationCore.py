# Keyhan Babaee, https://github.com/KeyhanB
# V1.2
# Oct 2019
from ConfigReaderModule import logging

from TMatrix import TMatrixCalculation
from FSAC_RDG import RDGCalculation


class CalculationCore:
    def __init__(self, infoDict, calcDict):
        try:
            self.__info = infoDict
            self.__calcDict = calcDict
        except Exception as e:
            logging.exception(e)
            raise

    def Calc(self):
        try:
            T1 = TMatrixCalculation(DBInfo=self.__info)
            R1 = RDGCalculation()

            dictInputT1, dictOutputT1, dictInputR1, dictOutputR1 = {}, {}, {}, {}
            for dm in self.__calcDict:
                dictInputT1[dm], dictOutputT1[dm] = T1.TMatrixCalc(TMatrix_Planned_Input=self.__calcDict[dm])
                dictInputR1[dm], dictOutputR1[dm] = R1.RDGCalc(RDG_Planned_Input=self.__calcDict[dm])
                if dictInputT1[dm] != dictInputR1[dm]:
                    raise
                logging.info(f"Calculation for dm:{dm} was finished.")
            return dictOutputT1, dictOutputR1
        except Exception as e:
            logging.exception(e)
            raise
