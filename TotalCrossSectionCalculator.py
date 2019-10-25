# Keyhan Babaee, https://github.com/KeyhanB
# V1.2
# Oct 2019
from ConfigReaderModule import logging
from decimal import Decimal
import pandas as pd
from math import *


class TotalCrossSectionCalculator:
    def __init__(self, infoDict, outTMatrix, outRDG, checkedDict, dict_dpMedianNano):
        try:
            self.__outTMatrix = outTMatrix
            self.__outRDG = outRDG
            self.__checkedDict = checkedDict
            self.__AGG_MATERIAL_DENSITY_CENTER = infoDict['AGG_MATERIAL_DENSITY_CENTER']
            self.D_TEM = infoDict['D_TEM']
            self.dp100_nano = infoDict['dp100_nano']
            self.dict_dpMedianNano = dict_dpMedianNano
        except Exception as e:
            logging.exception(e)
            raise

    def Calc(self):
        try:

            df = pd.DataFrame(columns=['dm', 'ABS_TMatrix', 'SCA_TMatrix', 'ABS_RDG', 'SCA_RDG',
                                       'NumberOfCalcs', 'dp_median', 'Np_Ave', 'dp_Ave',
                                       'D_TEM', 'dp_100nm',
                                       'RealErrorABS', 'RealErrorSCA', 'AbsoluteErrorABS', 'AbsoluteErrorSCA',
                                       'RealPercentErrorABS', 'AbsolutePercentErrorABS', 'RatioABS',
                                       'RealPercentErrorSCA', 'AbsolutePercentErrorSCA', 'RatioSCA',
                                       'SSA_TMatrix', 'SSA_RDG',
                                       'MAC_TMatrix', 'MSC_TMatrix', 'MAC_RDG', 'MSC_RDG', 'Agg_Mass_gr'])

            for dm in self.__outTMatrix:

                Np_Ave = 0
                dp_Ave = 0
                ABS_TMatrix = 0
                SCA_TMatrix = 0
                ABS_RDG = 0
                SCA_RDG = 0
                AggregateMass = 0

                for i in range(len(self.__outTMatrix[dm])):
                    ABS_TMatrix += self.__outTMatrix[dm][i][0] * self.__checkedDict[dm]['chance'][i]
                    SCA_TMatrix += self.__outTMatrix[dm][i][1] * self.__checkedDict[dm]['chance'][i]
                    ABS_RDG += self.__outRDG[dm][i][0] * self.__checkedDict[dm]['chance'][i]
                    SCA_RDG += self.__outRDG[dm][i][1] * self.__checkedDict[dm]['chance'][i]
                    Np_Ave += self.__checkedDict[dm]['Np'][i] * self.__checkedDict[dm]['chance'][i]
                    dp_Ave += self.__checkedDict[dm]['dp'][i] * self.__checkedDict[dm]['chance'][i]
                    AggregateMass += self._convertDensityToMass(dp=self.__checkedDict[dm]['dp'][i], density=self.__AGG_MATERIAL_DENSITY_CENTER, Np=self.__checkedDict[dm]['Np'][i]) * \
                                     self.__checkedDict[dm]['chance'][i]
                ################################################
                ################################################
                ################################################ MAC and MSC
                if AggregateMass != 0:
                    MAC_TMatrix = ABS_TMatrix * (Decimal(10) ** (Decimal(-12))) / AggregateMass
                    MSC_TMatrix = SCA_TMatrix * (Decimal(10) ** (Decimal(-12))) / AggregateMass
                    MAC_RDG = ABS_RDG * (Decimal(10) ** (Decimal(-12))) / AggregateMass
                    MSC_RDG = SCA_RDG * (Decimal(10) ** (Decimal(-12))) / AggregateMass

                else:
                    MAC_TMatrix = 0
                    MSC_TMatrix = 0
                    MAC_RDG = 0
                    MSC_RDG = 0

                errorRealABS = ABS_TMatrix - ABS_RDG
                errorRealSCA = SCA_TMatrix - SCA_RDG
                errorAbsoluteABS = abs(ABS_TMatrix - ABS_RDG)
                errorAbsoluteSCA = abs(SCA_TMatrix - SCA_RDG)
                ################################################
                ################################################
                ################################################ SSA
                if (ABS_TMatrix + SCA_TMatrix) != 0:
                    SSA_TMatrix = SCA_TMatrix / (ABS_TMatrix + SCA_TMatrix)
                else:
                    SSA_TMatrix = 0

                if (ABS_RDG + SCA_RDG) != 0:
                    SSA_RDG = SCA_RDG / (ABS_RDG + SCA_RDG)
                else:
                    SSA_RDG = 0
                ################################################
                ################################################
                ################################################ Error percenet and ratio
                if ABS_RDG != 0:
                    errorRealPercentABS = 100 * (ABS_TMatrix - ABS_RDG) / ABS_RDG
                    errorAbsolutePercentABS = 100 * abs(ABS_TMatrix - ABS_RDG) / ABS_RDG
                    ratioABS = ABS_TMatrix / ABS_RDG
                else:
                    errorRealPercentABS = 0
                    errorAbsolutePercentABS = 0
                    ratioABS = 0

                if SCA_RDG != 0:
                    errorRealPercentSCA = 100 * (SCA_TMatrix - SCA_RDG) / SCA_RDG
                    errorAbsolutePercentSCA = 100 * abs(SCA_TMatrix - SCA_RDG) / SCA_RDG
                    ratioSCA = SCA_TMatrix / SCA_RDG
                else:
                    errorRealPercentSCA = 0
                    errorAbsolutePercentSCA = 0
                    ratioSCA = 0
                ################################################
                ################################################
                ################################################
                df.loc[len(df)] = {'dm': dm, 'ABS_TMatrix': ABS_TMatrix, 'SCA_TMatrix': SCA_TMatrix,
                                   'ABS_RDG': ABS_RDG, 'SCA_RDG': SCA_RDG, 'NumberOfCalcs': len(self.__outTMatrix[dm]),
                                   'dp_median': self.dict_dpMedianNano[dm], 'Np_Ave': Np_Ave, 'dp_Ave': dp_Ave,
                                   'RealErrorABS': errorRealABS, 'RealErrorSCA': errorRealSCA,
                                   'D_TEM': self.D_TEM, 'dp_100nm': self.dp100_nano,
                                   'AbsoluteErrorABS': errorAbsoluteABS, 'AbsoluteErrorSCA': errorAbsoluteSCA,
                                   'RealPercentErrorABS': errorRealPercentABS, 'AbsolutePercentErrorABS': errorAbsolutePercentABS, 'RatioABS': ratioABS,
                                   'RealPercentErrorSCA': errorRealPercentSCA, 'AbsolutePercentErrorSCA': errorAbsolutePercentSCA, 'RatioSCA': ratioSCA,
                                   'SSA_TMatrix': SSA_TMatrix, 'SSA_RDG': SSA_RDG,
                                   'MAC_TMatrix': MAC_TMatrix, 'MSC_TMatrix': MSC_TMatrix, 'MAC_RDG': MAC_RDG, 'MSC_RDG': MSC_RDG,
                                   'Agg_Mass_gr': AggregateMass}
            return df

        except Exception as e:
            logging.exception(e)
            raise

    def _convertDensityToMass(self, dp, density, Np):
        try:
            V = (Decimal(10) ** (Decimal(-27))) * Decimal(pi) * (dp ** Decimal(3)) / Decimal(6)
            mass = Np * V * Decimal(density) * Decimal(1000)  # to gram
            return mass

        except Exception as e:
            logging.exception(e)
            raise
