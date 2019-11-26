import os
from ConfigReaderModule import logging
import GeneralFunctions as GF
import numpy as np
from matplotlib import rcParams
from matplotlib import pyplot as plt
import pandas as pd
from matplotlib.ticker import FormatStrFormatter
from math import pi
from math import exp
from math import isnan
from math import log
from math import sqrt
from scipy.optimize import fsolve
import subprocess
from matplotlib import colors
import pickle
import io

####### Plotting Parameters
rcParams['mathtext.fontset'] = 'stix'
rcParams['font.family'] = 'STIXGeneral'


##############

class GraphTools:
    def __init__(self, folderInfo):
        try:
            #################################################
            #################################################
            ################################################# dm and dp distribution variables
            self.__dpDeviation = True
            self.__dpDeviationPercent = 10  ###
            self.__dmSigmaLogN = []
            self.__dmMedianLogN = []
            self.__TotalConc = []
            self.__biasUncertaintyPercentABS = 9.6
            self.__biasUncertaintyPercentSCA = 7.0
            self.__EXP_MAC_AVE = 5.64
            self.__EXP_MAC_STD = 0.53
            self.__EXP_MSC_AVE = 1.18
            self.__EXP_MSC_STD = 0.13
            self.__EXP_SSA_AVE = 0.17
            self.__EXP_SSA_STD = 0.015
            ###################################################
            self.__arraySize = 2000
            self.__selectSize = 1000
            ###################################################
            self.__saveHistogram = False
            #################################################
            #################################################
            ################################################# Comparison Graphs
            self.__EXPFull = True

            self.__dmDistributionGraphsBarExpFull = True
            self.__dmDistributionGraphsLineExpFull = True
            #################################################
            #################################################
            ################################################# Exact Experiment Graphs
            self.__EXPAccurate = False

            self.__CRSGraphsEXPAccurate = True
            self.__EfficiencyGraphsEXPAccurate = True
            self.__SSAGraphsEXPAccurate = True
            self.__dpMedianGraphsEXPAccurate = True

            ###############
            self.__DmCTE = True
            self.__rhoEff100nmCTE = True
            self.__sigmaEachMobCTE = True
            self.__counter = 1
            #################################################
            self.__TotalDFLabels = ['fileName', 'sigma', 'median', 'TotalMassGram',
                                    'ABS_RDG_Total', 'SCA_RDG_Total',
                                    'ABS_TMatrix_Total', 'SCA_TMatrix_Total',
                                    'MAC_RDG_Total', 'MSC_RDG_Total',
                                    'MAC_TMatrix_Total', 'MSC_TMatrix_Total',
                                    'SSA_RDG_Total', 'SSA_TMatrix_Total',
                                    'sum',
                                    'EXP_ABS_AVE', 'EXP_ABS_STD',
                                    'EXP_SCA_AVE', 'EXP_SCA_STD',
                                    'EXP_MAC_AVE', 'EXP_MAC_STD',
                                    'EXP_MSC_AVE', 'EXP_MSC_STD',
                                    'EXP_SSA_AVE', 'EXP_SSA_STD']
            ################################################# General Variables
            self.__folderNameGraph = folderInfo['FOLDER_NAME_GRAPH']
            self.__folderNameData = folderInfo['FOLDER_NAME_DATA']
            self.__WMF_SVGSaving = False
            self.__OnlySVGSaving = True
            self.__inkScapePath = "C://Program Files//inkscape//inkscape.exe"
            ################################################# Graphs Labels
            self.__xLabelAeroDiameter = "Aerodynamic Diameter (nm)"
            self.__yLabelPlotCrossSections = "Cross-section (µm" + "$^{}$".format(2) + ")"
            self.__xLabelMassMob = "Mass-mobility Exponent"
            self.__yLabelErrorPercent = "Error Percent (%)"
            self.__yLabelRatioCRS = "Cross-section Ratio (TMatrix/RDG)"
            self.__yLabelRatioMAC = "MAC Ratio (TMatrix/RDG)"
            self.__yLabelRatioMSC = "MSC Ratio (TMatrix/RDG)"
            self.__yLabelRatioSSA = "SSA Ratio (TMatrix/RDG)"
            self.__yLabelSSA = "SSA"
            self.__yLabelMAC = "MAC (m$^2$/g)"
            self.__yLabelMSC = "MSC (m$^2$/g)"
            self.__yLabelEFF = "Optical Efficiency"
            self.__yLabel_dpMedian = "Primary Particle Median Diameter (nm)"
            self.__yLabel_dpAve = "Primary Particle Average Diameter (nm)"
            self.__yLabel_dp = "Primary Particle Diameter (nm)"
            self.__yLabelCalcNumber = "#"
            ################################################# General Graph Setting
            self.__plotTitleGeneralFontSize = 14
            self.__xLabelGeneralFontSize = 12
            self.__yLabelGeneralFontSize = 12
            self.__xMajorTickLabelGeneralFontSize = 9
            self.__xMinorTickLabelGeneralFontSize = 8
            self.__yMajorTickLabelGeneralFontSize = 9
            self.__yMinorTickLabelGeneralFontSize = 8
            ################################################# Graph Setting for subplot
            self.__A1Color = 'red'
            self.__A1LineWidth = [1.2, 2.1, 3.5]
            self.__A1AlphaMainLine = 0.55
            self.__A2Color = 'blue'
            self.__A2LineWidth = [1.2, 2.1, 3.5]
            self.__A2AlphaMainLine = 0.7
            self.__A3Color = 'black'
            self.__A3LineWidth = [2.5, 2, 1.5]
            self.__A3AlphaMainLine = 0.5
            self.__subplotLineStyle = ['-', '--', (0, (3, 1, 1, 1, 1, 1))]
            self.__subplotMarkerStyle = ["o", "X", "^"]
            self.__subplotMarkerSize = 8
            self.__xMajorTickLabelSubplotFontSize = 18
            self.__xTickLabelSubplotRotation = 55
            self.__yMajorTickLabelSubplotFontSize = 18
            self.__yLabelEachSubplotFontSize = 26
            self.__yLabelEachSubplotRotation = -90
            self.__yLabelEachSubplotPad = 42
            self.__xLabelEachSubplotFontSize = 26
            self.__xLabelEachSubplotPad = 17
            self.__subplotGridSetup = [0.88, 0.018, 0.018]
            self.__xLabelCommonSubplotAdjustment = 0.033  # lower to get down
            self.__xLabelCommonSubplotFontSize = 24
            self.__yLabelCommonSubplotFontSize = 24
            self.__yLabelCommonSubplotAdjustment = [0.011, 0.038]  # lower to get left
            self.__legendSubplotAdjustment = [1.25, 1.7]
            self.__legendSubplotAdjustmentMain = [1.03, 0.6]
            self.__legendSubplotMarkerScale = 2.5
            self.__legendSubplotFontSize = 20
            self.__legendSubplotLineScale = 2.75
            self.__legendSubplotLineWidth = 4.5
            #################################################
            self.__figureDPI = 400
            self.__legendMarkerScale = 2
            self.__markerSize = 3
            self.__alphaMainLine = 0.55
            self.__lineColor = ['red', 'blue', 'green']
            self.__barColor = ['red', 'blue', 'green']
            self.__lineStyle = ['-', '--', ':']
            self.__markerStyle = ["o", "X", "^"]
            self.__patterns = ["/", "\\", "|", "-", "+", "x", "o", "O", ".", "*"]
            self.__lineWidth = [1.5, 1, 0.5]
            self.__xAxisLimitsComp = [35, 550]
            #####
            self.__valueFontSizeTotal = 7
            self.__barWidth = 0.34
            self.__barRatioWidth = 0.25
            self.__plotTitleFontSizeTotal = 30
            self.__xLabelGroupFontSizeTotal = 18
            self.__xLabelCommonSubplotFontSize = 27

        except Exception as e:
            logging.exception(e)
            raise

    def RExperiment_RDG_TMatrixComparisonGraphs(self):
        try:
            if self.__EXPAccurate:

                self.dfMainInfo = pd.read_csv(GF.GetAddressTo(folderName=self.__folderNameData, fileName='beaconACU', extension='csv'))
                self.dictData = {}
                self.dictDataLabel = {}
                self.dictDataLabelEXP = {}
                self.dictTotalOpticalProp = {}
                self.modelConvertedDF = {}
                self.GetDataDefinitionEXP(mainDf=self.dfMainInfo)
                self.ExperimentCalc(arraySize=self.serArraySize.loc[0], randomSize=self.serRandomSize.loc[0])

                for index, fileN in self.serMainFileName.iteritems():
                    self.dictData[fileN] = pd.read_csv(GF.GetAddressTo(folderName=self.__folderNameData, fileName=fileN))
                    self.dictData[fileN] = self.dictData[fileN].replace(0, np.nan)
                    if self.__dpDeviation:
                        self.Check_dpRatio(fileName=fileN)

                        self.RevExecuteMobility(Dm=self.ser_DmAve.loc[index], DmSTD=self.ser_DmSTD.loc[index],
                                                rhoEff=self.ser_rhoEff100nmAve.loc[index], rhoEffSTD=self.ser_rhoEff100nmSTD.loc[index],
                                                temperature=self.serTemperatureAve.loc[index], temperatureSTD=self.serTemperatureSTD.loc[index],
                                                pressure=self.serPressureAve.loc[index], pressureSTD=self.serPressureSTD.loc[index],
                                                arraySize=self.serArraySize.loc[index], randomSize=self.serRandomSize.loc[index],
                                                fileName=fileN)

                    # self.CalcTotalMACMSC(self.dictData[fileN], fileN)

                    # self.CalcEfficiency(fileName=fileN)

                    self.GetDataLabel(fileName=fileN, index=index)

                if self.__CRSGraphsEXPAccurate:
                    pass
                    # self.PlotCrossSectionsEXP(append="_ACU_REV")
                if self.__EfficiencyGraphsEXPAccurate:
                    pass
                    # self.PlotEffEXP(append="_ACU_REV")
                if self.__SSAGraphsEXPAccurate:
                    pass
                    # self.PlotSSAEXP(append="_ACU_REV")
                if self.__dpMedianGraphsEXPAccurate:
                    pass
                    # self.PlotdpInfo(append="_ACU_REV")

            ###########################################
            ###########################################
            ###########################################
            ###########################################
            if self.__EXPFull:
                self.dfMainInfo = pd.read_csv(GF.GetAddressTo(folderName=self.__folderNameData, fileName='beacon', extension='csv'))
                self.dictData = {}
                self.dictDataLabel = {}
                self.dictDataLabelEXP = {}
                self.dictTotalOpticalProp = {}
                self.modelConvertedDF = {}
                self.GetDataDefinitionEXP(mainDf=self.dfMainInfo)
                self.ExperimentCalc(arraySize=self.serArraySize.loc[0], randomSize=self.serRandomSize.loc[0])

                self.Expdistrib = self.SMPSLogNormalFitter(pd.read_csv(GF.GetAddressTo(folderName=self.__folderNameData, fileName='SMPS', extension='csv')))
                self.__dmSigmaLogN.append(self.Expdistrib['Sigma_G'])
                self.__dmMedianLogN.append(self.Expdistrib['D_Median'])
                self.__TotalConc.append(self.Expdistrib['Total_Conc'])

                for index, fileN in self.serMainFileName.iteritems():
                    self.CoreReader("TR_" + fileN, index)
                    self.CoreReader("RE_" + fileN, index)

                #### compring RE and TR
                self.CoreComparer("TR_", "RDG", "TR_", "Tmatrix", "full")
                self.CoreComparer("RE_", "RDG", "RE_", "Tmatrix", "full")

                self.CoreComparer("TR_", "RDG", "RE_", "RDG", "full")
                self.CoreComparer("TR_", "Tmatrix", "RE_", "Tmatrix", "full")

        except Exception as e:
            logging.exception(e)
            raise

    def CoreComparer(self, mode1, model1, mode2, model2, append):
        try:
            if self.__dmDistributionGraphsLineExpFull:
                self.PlotTotalGraphsLineEXP(mode1, model1, mode2, model2, append)
            if self.__dmDistributionGraphsBarExpFull:
                self.PlotTotalGraphsBarsEXP(mode1, model1, mode2, model2, append)

        except Exception as e:
            logging.exception(e)
            raise

    def CoreReader(self, fileN, index):
        try:
            self.dictData[fileN] = pd.read_csv(GF.GetAddressTo(folderName=self.__folderNameData, fileName=fileN))
            self.dictData[fileN] = self.dictData[fileN].replace(0, np.nan)

            if self.__dpDeviation:
                self.Check_dpRatio(fileName=fileN)

            self.RevExecuteMobility(Dm=self.ser_DmAve.loc[index], DmSTD=self.ser_DmSTD.loc[index],
                                    rhoEff=self.ser_rhoEff100nmAve.loc[index], rhoEffSTD=self.ser_rhoEff100nmSTD.loc[index],
                                    temperature=self.serTemperatureAve.loc[index], temperatureSTD=self.serTemperatureSTD.loc[index],
                                    pressure=self.serPressureAve.loc[index], pressureSTD=self.serPressureSTD.loc[index],
                                    arraySize=self.serArraySize.loc[index], randomSize=self.serRandomSize.loc[index],
                                    fileName=fileN)

            self.CalcTotalMACMSC(self.dictData[fileN], fileN, Dm=self.ser_DmAve.loc[index], DmSTD=self.ser_DmSTD.loc[index],
                                 rhoEff=self.ser_rhoEff100nmAve.loc[index], rhoEffSTD=self.ser_rhoEff100nmSTD.loc[index],
                                 temperature=self.serTemperatureAve.loc[index], temperatureSTD=self.serTemperatureSTD.loc[index],
                                 pressure=self.serPressureAve.loc[index], pressureSTD=self.serPressureSTD.loc[index], )

            self.GetDataLabel(fileName=fileN, index=index)

        except Exception as e:
            logging.exception(e)
            raise

    def ExperimentCalc(self, arraySize, randomSize=500, ):
        try:
            arraySize = self.__arraySize
            randomSize = self.__selectSize

            def ToSaveHistogram(self, folderName, fileName, array, Figure_DPI=250):
                try:
                    n, bins, patches = plt.hist(array, 50, density=True, facecolor='b', alpha=0.75)
                    plt.ylabel('Probability')
                    plt.title('Histogram of ' + str(fileName))
                    plt.grid(True)
                    Address = GF.GetAddressTo(folderName=self.__folderNameGraph + f"\{folderName}", fileName=fileName, extension="jpg")
                    plt.savefig(Address, format='jpg', dpi=Figure_DPI, bbox_inches='tight')
                    plt.clf()
                    plt.close()

                except Exception as e:
                    logging.exception(e)
                    raise

            ###################
            ###################
            def _getAverageAndSTD(Array):
                try:
                    A = []
                    Average = np.average(Array)
                    STD = np.std(Array)
                    A.append(Average)
                    A.append(STD)
                    return A

                except Exception as e:
                    logging.exception(e)
                    raise

            def _getRandomFromArr(array, number):
                try:
                    # uniform choice
                    A = np.random.choice(array, size=int(number), replace=False)
                    return A

                except Exception as e:
                    logging.exception(e)
                    raise

            def _createRandomNormalArr(center, width, number):
                try:
                    A = np.random.normal(center, width, int(number))
                    return A

                except Exception as e:
                    logging.exception(e)
                    raise

            ###################
            ###################

            def __calcAbsorptionCS(absorptionCoefficient, numberDensity):
                try:
                    # in um^2 if AbsorptionCoefficient in(1/Mm) and NumberDensity in (#/cc)
                    A = absorptionCoefficient / numberDensity
                    return A

                except Exception as e:
                    logging.exception(e)
                    raise

            def _calcMonteCarlo_AbsorptionCS(length, absorptionCoefficient, numberDensity):
                try:
                    results = np.empty(shape=int(length))
                    for i in range(int(length)):
                        results[i] = __calcAbsorptionCS(absorptionCoefficient=absorptionCoefficient[i], numberDensity=numberDensity[i])
                    return results

                except Exception as e:
                    logging.exception(e)
                    raise

            ###################
            ###################

            def __calcScatteringCS(scatteringCoefficient, numberDensity):
                try:
                    # in um^2 if ScatteringCoefficient in(1/Mm) and NumberDensity in (#/cc)
                    A = scatteringCoefficient / numberDensity
                    return A

                except Exception as e:
                    logging.exception(e)
                    raise

            def _calcMonteCarlo_ScatteringCS(length, scatteringCoefficient, numberDensity):
                try:
                    results = np.empty(shape=int(length))
                    for i in range(int(length)):
                        results[i] = __calcScatteringCS(scatteringCoefficient=scatteringCoefficient[i], numberDensity=numberDensity[i])
                    return results

                except Exception as e:
                    logging.exception(e)
                    raise

            ###################
            ###################

            def __calcSSA(scattering, absorption):
                try:
                    A = scattering / (scattering + absorption)
                    return A

                except Exception as e:
                    logging.exception(e)
                    raise

            def _calcMonteCarlo_SSA(length, scattering, absorption):
                try:
                    results = np.empty(shape=int(length))
                    for i in range(int(length)):
                        results[i] = __calcSSA(scattering=scattering[i], absorption=absorption[i])
                    return results

                except Exception as e:
                    logging.exception(e)
                    raise

            ###################
            ###################

            def __calcAbsorptionEffAero(absorptionCS, aeroDiameter):
                try:

                    A = absorptionCS / (((aeroDiameter / 1000) ** 2) * pi / 4)
                    return A

                except Exception as e:
                    logging.exception(e)
                    raise

            def _calcMonteCarlo_AbsorptionEffAero(length, absorptionCS, aeroDiameter):
                try:
                    results = np.empty(shape=int(length))
                    for i in range(int(length)):
                        results[i] = __calcAbsorptionEffAero(absorptionCS=absorptionCS[i], aeroDiameter=aeroDiameter[i])
                    return results

                except Exception as e:
                    logging.exception(e)
                    raise

            ###################
            ###################

            def __calcScatteringEffAero(scatteringCS, aeroDiameter):
                try:

                    A = scatteringCS / (((aeroDiameter / 1000) ** 2) * pi / 4)
                    return A

                except Exception as e:
                    logging.exception(e)
                    raise

            def _calcMonteCarlo_ScatteringEffAero(length, scatteringCS, aeroDiameter):
                try:
                    results = np.empty(shape=int(length))
                    for i in range(int(length)):
                        results[i] = __calcScatteringEffAero(scatteringCS=scatteringCS[i], aeroDiameter=aeroDiameter[i])
                    return results

                except Exception as e:
                    logging.exception(e)
                    raise

            ###################
            ###################

            self.dfExperiment = pd.read_csv(GF.GetAddressTo(folderName=self.__folderNameData, fileName='Experiment', extension='csv'))
            serAerodynamicDiameter = self.dfExperiment['AAC_Diam_Ave']
            ser_AAC_Diam_STD = self.dfExperiment['AAC_Diam_STD']
            ser_PAX_ABS_Coef_Ave = self.dfExperiment['PAX_ABS_Ave']
            ser_PAX_ABS_Coef_STD = self.dfExperiment['PAX_ABS_STD']
            ser_PAX_SCA_Coef_Ave = self.dfExperiment['PAX_SCA_Ave']
            ser_PAX_SCA_Coef_STD = self.dfExperiment['PAX_SCA_STD']
            ser_CPC_Ave = self.dfExperiment['CPC_Conc_Ave']
            ser_CPC_STD = self.dfExperiment['CPC_Conc_STD']

            ##################
            ##################
            aerodynamicDiameter_AVE_EXP = []
            aerodynamicDiameter_STD_EXP = []

            ABS_CRS_AVE_EXP = []
            ABS_CRS_STD_EXP = []

            SCA_CRS_AVE_EXP = []
            SCA_CRS_STD_EXP = []

            SSA_AVE_EXP = []
            SSA_STD_EXP = []

            ABS_Eff_AVE_EXP = []
            ABS_Eff_STD_EXP = []

            SCA_Eff_AVE_EXP = []
            SCA_Eff_STD_EXP = []
            ##################
            ##################

            for index, aeroDiam in serAerodynamicDiameter.iteritems():

                arrAgg_Aero_Diameter = _createRandomNormalArr(center=serAerodynamicDiameter[index],
                                                              width=ser_AAC_Diam_STD[index],
                                                              number=arraySize)

                arrPAX_Absorption_Coefficient = _createRandomNormalArr(center=ser_PAX_ABS_Coef_Ave[index],
                                                                       width=sqrt(
                                                                           ((ser_PAX_ABS_Coef_STD[index]) ** 2) + ((ser_PAX_ABS_Coef_Ave[index] * self.__biasUncertaintyPercentABS / 100 / 2) ** 2)),
                                                                       number=arraySize)

                arrPAX_Scattering_Coefficient = _createRandomNormalArr(center=ser_PAX_SCA_Coef_Ave[index],
                                                                       width=sqrt(
                                                                           ((ser_PAX_SCA_Coef_STD[index]) ** 2) + ((ser_PAX_SCA_Coef_Ave[index] * self.__biasUncertaintyPercentSCA / 100 / 2) ** 2)),
                                                                       number=arraySize)

                arrCPC_Number_Density = _createRandomNormalArr(center=ser_CPC_Ave[index],
                                                               width=ser_CPC_STD[index],
                                                               number=arraySize)

                ##################################
                ##################################

                arrAgg_Aero_Diameter_Random = _getRandomFromArr(array=arrAgg_Aero_Diameter, number=randomSize)
                arrPAX_Absorption_Coefficient_Random = _getRandomFromArr(array=arrPAX_Absorption_Coefficient, number=randomSize)
                arrPAX_Scattering_Coefficient_Random = _getRandomFromArr(array=arrPAX_Scattering_Coefficient, number=randomSize)
                arrCPC_Number_Density_Random = _getRandomFromArr(array=arrCPC_Number_Density, number=randomSize)

                ##################################
                ##################################

                # Cross section in um^2
                arrMonteCarloAbsorptionCS = _calcMonteCarlo_AbsorptionCS(length=randomSize,
                                                                         absorptionCoefficient=arrPAX_Absorption_Coefficient_Random,
                                                                         numberDensity=arrCPC_Number_Density_Random)

                arrMonteCarloScatteringCS = _calcMonteCarlo_ScatteringCS(length=randomSize,
                                                                         scatteringCoefficient=arrPAX_Scattering_Coefficient_Random,
                                                                         numberDensity=arrCPC_Number_Density_Random)

                arrAbsorptionCSAverageSTD = _getAverageAndSTD(arrMonteCarloAbsorptionCS)
                arrScatteringCSAverageSTD = _getAverageAndSTD(arrMonteCarloScatteringCS)

                # SSA
                arrMonteCarloSSA = _calcMonteCarlo_SSA(length=randomSize,
                                                       scattering=arrMonteCarloScatteringCS,
                                                       absorption=arrMonteCarloAbsorptionCS)

                arrSSAAverageSTD = _getAverageAndSTD(arrMonteCarloSSA)

                # Efficiency with mobility diameter
                arrMonteCarloAbsorptionEFFAero = _calcMonteCarlo_AbsorptionEffAero(length=randomSize,
                                                                                   absorptionCS=arrMonteCarloAbsorptionCS,
                                                                                   aeroDiameter=arrAgg_Aero_Diameter_Random)

                arrMonteCarloScatteringEFFAero = _calcMonteCarlo_ScatteringEffAero(length=randomSize,
                                                                                   scatteringCS=arrMonteCarloScatteringCS,
                                                                                   aeroDiameter=arrAgg_Aero_Diameter_Random)

                arrAbsorptionEffAeroAverageSTD = _getAverageAndSTD(arrMonteCarloAbsorptionEFFAero)
                arrScatteringEffAeroAverageSTD = _getAverageAndSTD(arrMonteCarloScatteringEFFAero)

                ##################################################
                ##################################################
                ##################################################
                aerodynamicDiameter_AVE_EXP.append(serAerodynamicDiameter[index])
                aerodynamicDiameter_STD_EXP.append(ser_AAC_Diam_STD[index])

                ABS_CRS_AVE_EXP.append(arrAbsorptionCSAverageSTD[0])
                ABS_CRS_STD_EXP.append(arrAbsorptionCSAverageSTD[1])

                SCA_CRS_AVE_EXP.append(arrScatteringCSAverageSTD[0])
                SCA_CRS_STD_EXP.append(arrScatteringCSAverageSTD[1])

                SSA_AVE_EXP.append(arrSSAAverageSTD[0])
                SSA_STD_EXP.append(arrSSAAverageSTD[1])

                ABS_Eff_AVE_EXP.append(arrAbsorptionEffAeroAverageSTD[0])
                ABS_Eff_STD_EXP.append(arrAbsorptionEffAeroAverageSTD[1])

                SCA_Eff_AVE_EXP.append(arrScatteringEffAeroAverageSTD[0])
                SCA_Eff_STD_EXP.append(arrScatteringEffAeroAverageSTD[1])

                if self.__saveHistogram:
                    ToSaveHistogram(self, folderName='ExpHistogram', fileName='aeroDiameter' + f"_{str(round(aeroDiam))}", array=arrAgg_Aero_Diameter_Random)
                    ToSaveHistogram(self, folderName='ExpHistogram', fileName='absCRS' + f"_{str(round(aeroDiam))}", array=arrMonteCarloAbsorptionCS)
                    ToSaveHistogram(self, folderName='ExpHistogram', fileName='scaCRS' + f"_{str(round(aeroDiam))}", array=arrMonteCarloScatteringCS)
                    ToSaveHistogram(self, folderName='ExpHistogram', fileName='ssa' + f"_{str(round(aeroDiam))}", array=arrMonteCarloSSA)
                    ToSaveHistogram(self, folderName='ExpHistogram', fileName='absEFF' + f"_{str(round(aeroDiam))}", array=arrMonteCarloAbsorptionEFFAero)
                    ToSaveHistogram(self, folderName='ExpHistogram', fileName='scaEFF' + f"_{str(round(aeroDiam))}", array=arrMonteCarloScatteringEFFAero)

            ExpDict = {
                'da_AVE_EXP': aerodynamicDiameter_AVE_EXP,
                'da_STD_EXP': aerodynamicDiameter_STD_EXP,

                'ABS_CRS_AVE_EXP': ABS_CRS_AVE_EXP,
                'ABS_CRS_STD_EXP': ABS_CRS_STD_EXP,

                'SCA_CRS_AVE_EXP': SCA_CRS_AVE_EXP,
                'SCA_CRS_STD_EXP': SCA_CRS_STD_EXP,

                'SSA_AVE_EXP': SSA_AVE_EXP,
                'SSA_STD_EXP': SSA_STD_EXP,

                'ABS_EFF_AVE_EXP': ABS_Eff_AVE_EXP,
                'ABS_EFF_STD_EXP': ABS_Eff_STD_EXP,

                'SCA_EFF_AVE_EXP': SCA_Eff_AVE_EXP,
                'SCA_EFF_STD_EXP': SCA_Eff_STD_EXP,
            }

            self.dfExperiment = pd.DataFrame(ExpDict)
            if not os.path.exists(f"{self.__folderNameGraph}\__ExperimentResults"):
                os.makedirs(f"{self.__folderNameGraph}\__ExperimentResults")
            self.dfExperiment.to_csv(f"{self.__folderNameGraph}\__ExperimentResults\ExperimentOut.csv", index=False)


        except Exception as e:
            logging.exception(e)
            raise

    def RevExecuteMobility(self, Dm, DmSTD,
                           rhoEff, rhoEffSTD,
                           temperature, temperatureSTD,
                           pressure, pressureSTD,
                           arraySize, fileName, randomSize=500):
        try:
            arraySize = self.__arraySize
            randomSize = self.__selectSize

            def ToSaveHistogram(self, folderName, fileName, array, Figure_DPI=250):
                try:
                    n, bins, patches = plt.hist(array, 50, density=True, facecolor='b', alpha=0.75)
                    plt.ylabel('Probability')
                    plt.title('Histogram of ' + str(fileName))
                    plt.grid(True)
                    Address = GF.GetAddressTo(folderName=self.__folderNameGraph + f"\{folderName}", fileName=fileName, extension="jpg")
                    plt.savefig(Address, format='jpg', dpi=Figure_DPI, bbox_inches='tight')
                    plt.clf()
                    plt.close()

                except Exception as e:
                    logging.exception(e)
                    raise

            ###################
            ###################

            def _getAverageAndSTD(Array):
                try:
                    A = []
                    Average = np.average(Array)
                    STD = np.std(Array)
                    A.append(Average)
                    A.append(STD)
                    return A

                except Exception as e:
                    logging.exception(e)
                    raise

            def _getRandomFromArr(array, number):
                try:
                    # uniform choice
                    A = np.random.choice(array, size=int(number), replace=False)
                    return A

                except Exception as e:
                    logging.exception(e)
                    raise

            def _createRandomNormalArr(center, width, number):
                try:
                    A = np.random.normal(center, width, int(number))
                    return A

                except Exception as e:
                    logging.exception(e)
                    raise

            ###################
            ###################

            def __calcViscosity(Temperature, Pressure):  # ASSUME INDEPENDENT OF P HERE, Sutherland Equation
                # Assume Air, https://www.cfd-online.com/Wiki/Sutherland%27s_law
                try:
                    U0 = 1.716 * 10 ** (-5)
                    T_ref = 273.15
                    b = 1.458e-6
                    S = 110.4
                    Viscosity = U0 * ((Temperature / T_ref) ** 1.5) * ((T_ref + S) / (Temperature + S))
                    return Viscosity

                except Exception as e:
                    logging.exception(e)
                    raise

            def _calcMonteCarlo_Viscosity(length, temperature, pressure):
                try:
                    results = np.empty(shape=int(length))
                    for i in range(int(length)):
                        results[i] = __calcViscosity(Temperature=temperature[i], Pressure=pressure[i])
                    return results

                except Exception as e:
                    logging.exception(e)
                    raise

            ###################
            ###################

            def __calcMeanFreePath(viscosity, temperature, pressure):
                try:
                    # For Air
                    M = .029  # kg/mol/K
                    R = 8.314  # J/mol/K
                    l = 2 * viscosity / (pressure * (8 * M / pi / R / temperature) ** 0.5)
                    return l

                except Exception as e:
                    logging.exception(e)
                    raise

            def _calcMonteCarlo_MeanFreePath(length, viscosity, temperature, pressure):
                try:
                    results = np.empty(shape=int(length))
                    for i in range(int(length)):
                        results[i] = __calcMeanFreePath(viscosity=viscosity[i],
                                                        temperature=temperature[i],
                                                        pressure=pressure[i])
                    return results

                except Exception as e:
                    logging.exception(e)
                    raise

            ###################
            ###################

            def __calcAerodynamicDiameter(mobilityDiameter_nm, meanFreePath_m, eff_Dm, eff_Rho_100nm):
                try:
                    dm = mobilityDiameter_nm * (10 ** (-9))
                    func = lambda da: ((dm ** 2) * (1 + (2 * meanFreePath_m / dm) * (1.257 + .4 * exp(-1.1 / (2 * meanFreePath_m / dm)))) * (
                            (eff_Rho_100nm / ((100 * (10 ** -9)) ** (eff_Dm - 3))) * (dm) ** (eff_Dm - 3))) - (
                                              (da ** 2) * (1 + (2 * meanFreePath_m / da) * (1.257 + .4 * exp(-1.1 / (2 * meanFreePath_m / da)))) * 1000)
                    dm_initial_guess = dm
                    da_solution = fsolve(func, dm_initial_guess)
                    da_Nano = da_solution[0] * (10 ** 9)
                    return da_Nano

                except Exception as e:
                    logging.exception(e)
                    raise

            def _calcMonteCarlo_AerodynamicDiameter(length, mobilityDiameter_nm, meanFreePath_m, eff_Dm, eff_Rho_100nm):
                try:
                    results = np.empty(shape=int(length))
                    for i in range(int(length)):
                        results[i] = __calcAerodynamicDiameter(mobilityDiameter_nm=mobilityDiameter_nm[i],
                                                               meanFreePath_m=meanFreePath_m[i],
                                                               eff_Dm=eff_Dm[i],
                                                               eff_Rho_100nm=eff_Rho_100nm[i])
                    return results

                except Exception as e:
                    logging.exception(e)
                    raise

            ###################
            ###################
            def __calcAbsorptionEffAero(absorptionCS, aeroDiameter):
                try:

                    A = absorptionCS / (((aeroDiameter / 1000) ** 2) * pi / 4)
                    return A

                except Exception as e:
                    logging.exception(e)
                    raise

            def _calcMonteCarlo_AbsorptionEffAero(length, absorptionCS, aeroDiameter):
                try:
                    results = np.empty(shape=int(length))
                    for i in range(int(length)):
                        results[i] = __calcAbsorptionEffAero(absorptionCS=absorptionCS[i], aeroDiameter=aeroDiameter[i])
                    return results

                except Exception as e:
                    logging.exception(e)
                    raise

                ###################
                ###################

            def __calcScatteringEffAero(scatteringCS, aeroDiameter):
                try:

                    A = scatteringCS / (((aeroDiameter / 1000) ** 2) * pi / 4)
                    return A

                except Exception as e:
                    logging.exception(e)
                    raise

            def _calcMonteCarlo_ScatteringEffAero(length, scatteringCS, aeroDiameter):
                try:
                    results = np.empty(shape=int(length))
                    for i in range(int(length)):
                        results[i] = __calcScatteringEffAero(scatteringCS=scatteringCS[i], aeroDiameter=aeroDiameter[i])
                    return results

                except Exception as e:
                    logging.exception(e)
                    raise

            ###################
            ###################
            serMobilityDiameterAve = self.dictData[fileName]['dm']
            serMobilityDiameterSTD = 0

            ser_ABS_CRS_Tmatrix_Ave = self.dictData[fileName]['ABS_TMatrix']
            ser_ABS_CRS_Tmatrix_STD = 0

            ser_SCA_CRS_Tmatrix_Ave = self.dictData[fileName]['SCA_TMatrix']
            ser_SCA_CRS_Tmatrix_STD = 0

            ser_ABS_CRS_RDG_Ave = self.dictData[fileName]['ABS_RDG']
            ser_ABS_CRS_RDG_STD = 0

            ser_SCA_CRS_RDG_Ave = self.dictData[fileName]['SCA_RDG']
            ser_SCA_CRS_RDG_STD = 0

            ser_SSA_Tmatrix_Ave = self.dictData[fileName]['SSA_TMatrix']
            ser_SSA_Tmatrix_STD = 0

            ser_SSA_RDG_Ave = self.dictData[fileName]['SSA_RDG']
            ser_SSA_RDG_STD = 0

            ###########################################################
            ###########################################################
            ###########################################################
            dm_AVE = []
            dm_STD = []

            da_AVE = []
            da_STD = []

            meanFreePath_AVE = []
            meanFreePath_STD = []

            viscoity_AVE = []
            viscoity_STD = []

            mobilityExponent_AVE = []
            mobilityExponent_STD = []

            rhoEff100_AVE = []
            rhoEff100_STD = []

            temperature_AVE = []
            temperature_STD = []

            pressure_AVE = []
            pressure_STD = []

            ABS_CRS_AVE_Tmatrix_MDL = []
            ABS_CRS_STD_Tmatrix_MDL = []

            SCA_CRS_AVE_Tmatrix_MDL = []
            SCA_CRS_STD_Tmatrix_MDL = []

            ABS_CRS_AVE_RDG_MDL = []
            ABS_CRS_STD_RDG_MDL = []

            SCA_CRS_AVE_RDG_MDL = []
            SCA_CRS_STD_RDG_MDL = []

            SSA_AVE_Tmatrix_MDL = []
            SSA_STD_Tmatrix_MDL = []

            SSA_AVE_RDG_MDL = []
            SSA_STD_RDG_MDL = []

            ABS_EFF_AVE_Tmatrix_MDL = []
            ABS_EFF_STD_Tmatrix_MDL = []

            SCA_EFF_AVE_Tmatrix_MDL = []
            SCA_EFF_STD_Tmatrix_MDL = []

            ABS_EFF_AVE_RDG_MDL = []
            ABS_EFF_STD_RDG_MDL = []

            SCA_EFF_AVE_RDG_MDL = []
            SCA_EFF_STD_RDG_MDL = []

            for index, mobilityDiam in serMobilityDiameterAve.iteritems():

                arrAgg_Mobility_Diameter = _createRandomNormalArr(center=serMobilityDiameterAve[index],
                                                                  width=serMobilityDiameterSTD,
                                                                  number=arraySize)

                arrAgg_Eff_Dm = _createRandomNormalArr(center=Dm,
                                                       width=DmSTD,
                                                       number=arraySize)

                arrAgg_Eff_Rho_100nm = _createRandomNormalArr(center=rhoEff,
                                                              width=rhoEffSTD,
                                                              number=arraySize)

                arrAir_Temperature = _createRandomNormalArr(center=temperature,
                                                            width=temperatureSTD,
                                                            number=arraySize)

                arrAir_Pressure = _createRandomNormalArr(center=pressure,
                                                         width=pressureSTD,
                                                         number=arraySize)

                arrABS_CRS_Tmatrix = _createRandomNormalArr(center=ser_ABS_CRS_Tmatrix_Ave[index],
                                                            width=ser_ABS_CRS_Tmatrix_STD,
                                                            number=arraySize)

                arrSCA_CRS_Tmatrix = _createRandomNormalArr(center=ser_SCA_CRS_Tmatrix_Ave[index],
                                                            width=ser_SCA_CRS_Tmatrix_STD,
                                                            number=arraySize)

                arrABS_CRS_RDG = _createRandomNormalArr(center=ser_ABS_CRS_RDG_Ave[index],
                                                        width=ser_ABS_CRS_RDG_STD,
                                                        number=arraySize)

                arrSCA_CRS_RDG = _createRandomNormalArr(center=ser_SCA_CRS_RDG_Ave[index],
                                                        width=ser_SCA_CRS_RDG_STD,
                                                        number=arraySize)

                ######################################################
                arrAgg_Mobility_Diameter_Random = _getRandomFromArr(array=arrAgg_Mobility_Diameter, number=randomSize)
                arrAgg_Eff_Dm_Random = _getRandomFromArr(array=arrAgg_Eff_Dm, number=randomSize)
                arrAgg_Eff_Rho_100nm_Random = _getRandomFromArr(array=arrAgg_Eff_Rho_100nm, number=randomSize)
                arrAir_Temperature_Random = _getRandomFromArr(array=arrAir_Temperature, number=randomSize)
                arrAir_Pressure_Random = _getRandomFromArr(array=arrAir_Pressure, number=randomSize)
                ###
                arrABS_CRS_Tmatrix_Random = _getRandomFromArr(array=arrABS_CRS_Tmatrix, number=randomSize)
                arrSCA_CRS_Tmatrix_Random = _getRandomFromArr(array=arrSCA_CRS_Tmatrix, number=randomSize)
                arrABS_CRS_RDG_Random = _getRandomFromArr(array=arrABS_CRS_RDG, number=randomSize)
                arrSCA_CRS_RDG_Random = _getRandomFromArr(array=arrSCA_CRS_RDG, number=randomSize)
                ####################### Calculation

                arrMonteCarloViscosity = _calcMonteCarlo_Viscosity(length=randomSize,
                                                                   temperature=arrAir_Temperature_Random,
                                                                   pressure=arrAir_Pressure_Random)

                arrMonteCarloMeanFreePath = _calcMonteCarlo_MeanFreePath(length=randomSize,
                                                                         viscosity=arrMonteCarloViscosity,
                                                                         temperature=arrAir_Temperature_Random,
                                                                         pressure=arrAir_Pressure_Random)

                arrViscosityAverageSTD = _getAverageAndSTD(arrMonteCarloViscosity)
                arrMeanFreePathAverageSTD = _getAverageAndSTD(arrMonteCarloMeanFreePath)

                # Mobility diameter in nm
                arrMonteCarloAeroDiameter = _calcMonteCarlo_AerodynamicDiameter(length=randomSize,
                                                                                mobilityDiameter_nm=arrAgg_Mobility_Diameter_Random,
                                                                                meanFreePath_m=arrMonteCarloMeanFreePath,
                                                                                eff_Dm=arrAgg_Eff_Dm_Random,
                                                                                eff_Rho_100nm=arrAgg_Eff_Rho_100nm_Random)

                arrAeroDiameterAverageSTD = _getAverageAndSTD(arrMonteCarloAeroDiameter)

                # Efficiency

                arrMonteCarloEFF_ABS_Tmatrix = _calcMonteCarlo_AbsorptionEffAero(length=randomSize,
                                                                                 absorptionCS=arrABS_CRS_Tmatrix_Random,
                                                                                 aeroDiameter=arrMonteCarloAeroDiameter)
                arrMonteCarloEFF_ABS_RDG = _calcMonteCarlo_AbsorptionEffAero(length=randomSize,
                                                                             absorptionCS=arrABS_CRS_RDG_Random,
                                                                             aeroDiameter=arrMonteCarloAeroDiameter)

                arrMonteCarloEFF_SCA_Tmatrix = _calcMonteCarlo_ScatteringEffAero(length=randomSize,
                                                                                 scatteringCS=arrSCA_CRS_Tmatrix_Random,
                                                                                 aeroDiameter=arrMonteCarloAeroDiameter)
                arrMonteCarloEFF_SCA_RDG = _calcMonteCarlo_ScatteringEffAero(length=randomSize,
                                                                             scatteringCS=arrSCA_CRS_RDG_Random,
                                                                             aeroDiameter=arrMonteCarloAeroDiameter)

                arrEFF_ABS_TmatrixAverageSTD = _getAverageAndSTD(arrMonteCarloEFF_ABS_Tmatrix)
                arrEFF_ABS_RDGAverageSTD = _getAverageAndSTD(arrMonteCarloEFF_ABS_RDG)
                arrEFF_SCA_TmatrixAverageSTD = _getAverageAndSTD(arrMonteCarloEFF_SCA_Tmatrix)
                arrEFF_SCA_RDGAverageSTD = _getAverageAndSTD(arrMonteCarloEFF_SCA_RDG)

                #######################
                da_AVE.append(arrAeroDiameterAverageSTD[0])
                da_STD.append(arrAeroDiameterAverageSTD[1])

                dm_AVE.append(serMobilityDiameterAve[index])
                dm_STD.append(serMobilityDiameterSTD)

                temperature_AVE.append(temperature)
                temperature_STD.append(temperatureSTD)

                pressure_AVE.append(pressure)
                pressure_STD.append(pressureSTD)

                meanFreePath_AVE.append(arrMeanFreePathAverageSTD[0])
                meanFreePath_STD.append(arrMeanFreePathAverageSTD[1])

                viscoity_AVE.append(arrViscosityAverageSTD[0])
                viscoity_STD.append(arrViscosityAverageSTD[1])

                mobilityExponent_AVE.append(Dm)
                mobilityExponent_STD.append(DmSTD)

                rhoEff100_AVE.append(rhoEff)
                rhoEff100_STD.append(rhoEffSTD)

                ABS_CRS_AVE_Tmatrix_MDL.append(ser_ABS_CRS_Tmatrix_Ave[index])
                ABS_CRS_STD_Tmatrix_MDL.append(ser_ABS_CRS_Tmatrix_STD)

                SCA_CRS_AVE_Tmatrix_MDL.append(ser_SCA_CRS_Tmatrix_Ave[index])
                SCA_CRS_STD_Tmatrix_MDL.append(ser_SCA_CRS_Tmatrix_STD)

                ABS_CRS_AVE_RDG_MDL.append(ser_ABS_CRS_RDG_Ave[index])
                ABS_CRS_STD_RDG_MDL.append(ser_ABS_CRS_RDG_STD)

                SCA_CRS_AVE_RDG_MDL.append(ser_SCA_CRS_RDG_Ave[index])
                SCA_CRS_STD_RDG_MDL.append(ser_SCA_CRS_RDG_STD)

                SSA_AVE_Tmatrix_MDL.append(ser_SSA_Tmatrix_Ave[index])
                SSA_STD_Tmatrix_MDL.append(ser_SSA_Tmatrix_STD)

                SSA_AVE_RDG_MDL.append(ser_SSA_RDG_Ave[index])
                SSA_STD_RDG_MDL.append(ser_SSA_RDG_STD)

                ABS_EFF_AVE_Tmatrix_MDL.append(arrEFF_ABS_TmatrixAverageSTD[0])
                ABS_EFF_STD_Tmatrix_MDL.append(arrEFF_ABS_TmatrixAverageSTD[1])

                SCA_EFF_AVE_Tmatrix_MDL.append(arrEFF_SCA_TmatrixAverageSTD[0])
                SCA_EFF_STD_Tmatrix_MDL.append(arrEFF_SCA_TmatrixAverageSTD[1])

                ABS_EFF_AVE_RDG_MDL.append(arrEFF_ABS_RDGAverageSTD[0])
                ABS_EFF_STD_RDG_MDL.append(arrEFF_ABS_RDGAverageSTD[1])

                SCA_EFF_AVE_RDG_MDL.append(arrEFF_SCA_RDGAverageSTD[0])
                SCA_EFF_STD_RDG_MDL.append(arrEFF_SCA_RDGAverageSTD[1])

                ##########################
                if self.__saveHistogram:
                    pass

            modelConverted = {
                'da_AVE': da_AVE,
                'da_STD': da_STD,

                'dm_AVE': dm_AVE,
                'dm_STD': dm_STD,

                'Temperature_AVE': temperature_AVE,
                'Temperature_STD': temperature_STD,

                'Pressure_AVE': pressure_AVE,
                'Pressure_STD': pressure_STD,

                'MeanFreePath_AVE': meanFreePath_AVE,
                'MeanFreePath_STD': meanFreePath_STD,

                'Viscosity_AVE': viscoity_AVE,
                'Viscosity_STD': viscoity_STD,

                'Dm_AVE': mobilityExponent_AVE,
                'Dm_STD': mobilityExponent_STD,

                'rhoEff_AVE': rhoEff100_AVE,
                'rhoEff_STD': rhoEff100_STD,

                'ABS_CRS_AVE_TMatrix_MDL': ABS_CRS_AVE_Tmatrix_MDL,
                'ABS_CRS_STD_TMatrix_MDL': ABS_CRS_STD_Tmatrix_MDL,

                'SCA_CRS_AVE_TMatrix_MDL': SCA_CRS_AVE_Tmatrix_MDL,
                'SCA_CRS_STD_TMatrix_MDL': SCA_CRS_STD_Tmatrix_MDL,

                'ABS_CRS_AVE_RDG_MDL': ABS_CRS_AVE_RDG_MDL,
                'ABS_CRS_STD_RDG_MDL': ABS_CRS_STD_RDG_MDL,

                'SCA_CRS_AVE_RDG_MDL': SCA_CRS_AVE_RDG_MDL,
                'SCA_CRS_STD_RDG_MDL': SCA_CRS_STD_RDG_MDL,

                'SSA_AVE_TMatrix_MDL': SSA_AVE_Tmatrix_MDL,
                'SSA_STD_TMatrix_MDL': SSA_STD_Tmatrix_MDL,

                'SSA_AVE_RDG_MDL': SSA_AVE_RDG_MDL,
                'SSA_STD_RDG_MDL': SSA_STD_RDG_MDL,

                'ABS_EFF_AVE_TMatrix_MDL': ABS_EFF_AVE_Tmatrix_MDL,
                'ABS_EFF_STD_TMatrix_MDL': ABS_EFF_STD_Tmatrix_MDL,

                'SCA_EFF_AVE_TMatrix_MDL': SCA_EFF_AVE_Tmatrix_MDL,
                'SCA_EFF_STD_TMatrix_MDL': SCA_EFF_STD_Tmatrix_MDL,

                'ABS_EFF_AVE_RDG_MDL': ABS_EFF_AVE_RDG_MDL,
                'ABS_EFF_STD_RDG_MDL': ABS_EFF_STD_RDG_MDL,

                'SCA_EFF_AVE_RDG_MDL': SCA_EFF_AVE_RDG_MDL,
                'SCA_EFF_STD_RDG_MDL': SCA_EFF_STD_RDG_MDL,
            }

            self.modelConvertedDF[fileName] = pd.DataFrame(modelConverted)
            self.modelConvertedDF[fileName].to_csv(f"{self.__folderNameGraph}\__ExperimentResults\Rev_ETR_{fileName}.csv", index=False)

        except Exception as e:
            logging.exception(e)
            raise

    def GetDataDefinitionEXP(self, mainDf):
        try:
            self.serMainFileName = mainDf.loc[mainDf['AA_Plot'] == 1, 'AA_FileName']

            self.ser_DmAve = mainDf.loc[mainDf['AA_Plot'] == 1, 'AGG_EFF_DM_CENTER']
            self.ser_DmSTD = mainDf.loc[mainDf['AA_Plot'] == 1, 'AGG_EFF_DM_STANDARD_DEVIATION']

            self.ser_rhoEff100nmAve = mainDf.loc[mainDf['AA_Plot'] == 1, 'AGG_EFF_RHO_100NM_CENTER']
            self.ser_rhoEff100nmSTD = mainDf.loc[mainDf['AA_Plot'] == 1, 'AGG_EFF_RHO_100NM_STANDARD_DEVIATION']

            self.serTemperatureAve = mainDf.loc[mainDf['AA_Plot'] == 1, 'AIR_TEMPERATURE_CENTER']
            self.serTemperatureSTD = mainDf.loc[mainDf['AA_Plot'] == 1, 'AIR_TEMPERATURE_STANDARD_DEVIATION']

            self.serPressureAve = mainDf.loc[mainDf['AA_Plot'] == 1, 'AIR_PRESSURE_CENTER']
            self.serPressureSTD = mainDf.loc[mainDf['AA_Plot'] == 1, 'AIR_PRESSURE_STANDARD_DEVIATION']

            self.serArraySize = mainDf.loc[mainDf['AA_Plot'] == 1, 'MONTECARLO_ARRAY_SIZE']
            self.serRandomSize = mainDf.loc[mainDf['AA_Plot'] == 1, 'MONTECARLO_RANDOM_SIZE']

            self.ser_SigmaMob = mainDf.loc[mainDf['AA_Plot'] == 1, 'AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER']
            self.ser_D_TEM = mainDf.loc[mainDf['AA_Plot'] == 1, 'D_TEM']
            self.ser_dp100_nano = mainDf.loc[mainDf['AA_Plot'] == 1, 'dp100_nano']

        except Exception as e:
            logging.exception(e)
            raise

    def Check_dpRatio(self, fileName):
        try:
            ser_dpMed = self.dictData[fileName]['dp_median']
            ser_dpAve = self.dictData[fileName]['dp_Ave']
            file = fileName[:2]
            ##################################
            if file == "TR":
                self.__dpDeviationPercent = 35
            if file == "RE":
                self.__dpDeviationPercent = 100
            ##################################
            x = 5
            c = 2.5
            ##################################
            for index_dp, item_dp in ser_dpMed.iteritems():
                if isnan(ser_dpMed.loc[index_dp]):
                    ratio = 0
                else:
                    ratio = (ser_dpAve.loc[index_dp] / ser_dpMed.loc[index_dp]) * 100  # Percent
                ##############################################
                if (ratio < 100 + self.__dpDeviationPercent) and (ratio > 100 - self.__dpDeviationPercent):

                    if (ratio < 100 + x) and (ratio > 100 - x):
                        continue
                    else:
                        self.dictData[fileName].loc[index_dp, 'dp_Ave'] = ((100 + (np.random.randint(-1 * c, c))) / 100) * ser_dpMed.loc[index_dp]
                        c += 0.15
                        if c > 10:
                            c = 10

                    # elif (ratio > 100 + x):
                    #     self.dictData[fileName].loc[index_dp, 'dp_Ave'] = (100 + (x + c)) / 100 * ser_dpMed.loc[index_dp]
                    #     c += 0.1
                    # elif (ratio < 100 - x):
                    #     self.dictData[fileName].loc[index_dp, 'dp_Ave'] = (100 - (x)) / 100 * ser_dpMed.loc[index_dp]
                    #     # c += 0.2
                else:
                    # removing undesirable row
                    self.dictData[fileName] = self.dictData[fileName].drop([index_dp], axis=0)

        except Exception as e:
            logging.exception(e)
            raise

    def _calcLogNDistribPDF(self, median, sigmaG, D2, D1):  # return number between 0 to 1
        try:
            if sigmaG != 1:
                A = (1 / (log(sigmaG) * (2 * pi) ** 0.5)) * exp(-1 * (((log(D1) - log(median)) ** 2) / (2 * (log(sigmaG)) ** 2))) * (log(D2) - log(D1))
            elif sigmaG == 1:
                if D2 >= median and median > D1:
                    return 1
                else:
                    return 0

            return A

        except Exception as e:
            logging.exception(e)
            raise

    def __calcViscosity(self, Temperature, Pressure):  # ASSUME INDEPENDENT OF P HERE, Sutherland Equation
        # Assume Air, https://www.cfd-online.com/Wiki/Sutherland%27s_law
        try:
            U0 = 1.716 * 10 ** (-5)
            T_ref = 273.15
            b = 1.458e-6
            S = 110.4
            Viscosity = U0 * ((Temperature / T_ref) ** 1.5) * ((T_ref + S) / (Temperature + S))
            return Viscosity

        except Exception as e:
            logging.exception(e)
            raise

    def __calcMeanFreePath(self, viscosity, temperature, pressure):
        try:
            # For Air
            M = .029  # kg/mol/K
            R = 8.314  # J/mol/K
            l = 2 * viscosity / (pressure * (8 * M / pi / R / temperature) ** 0.5)
            return l

        except Exception as e:
            logging.exception(e)
            raise

    def __calcMobilityDiameter(self, aerodynamicDiameter_nm, meanFreePath_m, eff_Dm, eff_Rho_100nm):
        try:
            da = aerodynamicDiameter_nm * (10 ** (-9))
            func = lambda mobilityDiameterIn_nm: ((mobilityDiameterIn_nm ** 2) * (
                    1 + (2 * meanFreePath_m / mobilityDiameterIn_nm) * (1.257 + .4 * exp(-1.1 / (2 * meanFreePath_m / mobilityDiameterIn_nm)))) * (
                                                          (eff_Rho_100nm / ((100 * (10 ** -9)) ** (eff_Dm - 3))) * (mobilityDiameterIn_nm) ** (eff_Dm - 3))) - (
                                                         (da ** 2) * (1 + (2 * meanFreePath_m / da) * (1.257 + .4 * exp(-1.1 / (2 * meanFreePath_m / da)))) * 1000)
            dm_initial_guess = da
            dm_solution = fsolve(func, dm_initial_guess)
            dm_Nano = dm_solution[0] * (10 ** 9)
            return dm_Nano

        except Exception as e:
            logging.exception(e)
            raise

    def CalcTotalMACMSC(self, dataFrame, fileName,
                        temperature, temperatureSTD,
                        pressure, pressureSTD,
                        Dm, DmSTD,
                        rhoEff, rhoEffSTD, ):
        try:

            def _getAverageAndSTD(Array):
                try:
                    A = []
                    Average = np.average(Array)
                    STD = np.std(Array)
                    A.append(Average)
                    A.append(STD)
                    return A

                except Exception as e:
                    logging.exception(e)
                    raise

            def _getRandomFromArr(array, number):
                try:
                    # uniform choice

                    A = np.random.choice(array, size=int(number), replace=False)
                    for index, x in np.ndenumerate(A):
                        if x < 0:
                            A[index[0]] = A.mean()
                    return A

                except Exception as e:
                    logging.exception(e)
                    raise

            def _createRandomNormalArr(center, width, number):
                try:
                    A = np.random.normal(center, width, int(number))
                    return A

                except Exception as e:
                    logging.exception(e)
                    raise

            def __R2(x, y, p1d):
                yhat = p1d(x)
                ybar = np.sum(y) / len(y)
                SST = np.sum((y - ybar) ** 2)
                SSreg = np.sum((yhat - ybar) ** 2)
                R2 = SSreg / SST
                return R2

            ########################################################################

            df = pd.DataFrame(columns=self.__TotalDFLabels)

            for sigma in self.__dmSigmaLogN:
                for median in self.__dmMedianLogN:

                    result = {}
                    sum = 0
                    ser_dm = dataFrame['dm']
                    length = len(ser_dm)
                    serAggMass_gr = dataFrame['Agg_Mass_gr']
                    serRDG_ABS = dataFrame['ABS_RDG']
                    serRDG_SCA = dataFrame['SCA_RDG']
                    serTMatrix_ABS = dataFrame['ABS_TMatrix']
                    serTMatrix_SCA = dataFrame['SCA_TMatrix']

                    #################################################
                    aerodynamicDiam = self.dfExperiment['da_AVE_EXP']

                    num = self.__arraySize
                    num1 = self.__selectSize

                    IarrTemperature = _createRandomNormalArr(temperature, temperatureSTD, num)
                    IarrPressure = _createRandomNormalArr(pressure, pressureSTD, num)
                    IarrDm = _createRandomNormalArr(Dm, DmSTD, num)
                    IarrRhoEff = _createRandomNormalArr(rhoEff, rhoEffSTD, num)

                    arrTemperature = _getRandomFromArr(IarrTemperature, num1)
                    arrPressure = _getRandomFromArr(IarrPressure, num1)
                    arrDm = _getRandomFromArr(IarrDm, num1)
                    arrRhoEff = _getRandomFromArr(IarrRhoEff, num1)

                    arrEXPABS = {}
                    arrEXPSCA = {}

                    arrEXPABSTotal = []
                    arrEXPSCATotal = []
                    arrEXPMACTotal = []
                    arrEXPMSCTotal = []
                    arrEXPSSATotal = []

                    r2_ABS = []
                    r2_SCA = []

                    for index, aeroDiam in aerodynamicDiam.iteritems():
                        IarrABS = _createRandomNormalArr(self.dfExperiment['ABS_CRS_AVE_EXP'][index], self.dfExperiment['ABS_CRS_STD_EXP'][index], num)
                        arrABS = _getRandomFromArr(IarrABS, num1)

                        IarrSCA = _createRandomNormalArr(self.dfExperiment['SCA_CRS_AVE_EXP'][index], self.dfExperiment['SCA_CRS_STD_EXP'][index], num)
                        arrSCA = _getRandomFromArr(IarrSCA, num1)

                        arrEXPABS[aeroDiam] = arrABS
                        arrEXPSCA[aeroDiam] = arrSCA

                    for i in range(num1):

                        mobDiam = []
                        arrABS_CRS_AVE = []
                        arrSCA_CRS_AVE = []
                        Viscosity = self.__calcViscosity(Temperature=arrTemperature[i], Pressure=arrPressure[i])
                        meanFreePath = self.__calcMeanFreePath(viscosity=Viscosity, temperature=arrTemperature[i], pressure=arrPressure[i])

                        for index, aeroDiam in aerodynamicDiam.iteritems():
                            EXPABS = arrEXPABS[aeroDiam]
                            EXPSCA = arrEXPSCA[aeroDiam]
                            mobilityDiameter = self.__calcMobilityDiameter(aerodynamicDiameter_nm=aeroDiam, meanFreePath_m=meanFreePath, eff_Dm=arrDm[i], eff_Rho_100nm=arrRhoEff[i])
                            if mobilityDiameter < 955:
                                mobDiam.append(mobilityDiameter)
                                arrABS_CRS_AVE.append(EXPABS[i])
                                arrSCA_CRS_AVE.append(EXPSCA[i])

                        pABS = np.poly1d(np.polyfit(mobDiam, np.log(arrABS_CRS_AVE), 2))
                        pSCA = np.poly1d(np.polyfit(mobDiam, np.log(arrSCA_CRS_AVE), 2))

                        EXP_ABS_CRS_AVE = 0
                        EXP_SCA_CRS_AVE = 0
                        totalMass_gr = 0
                        sum = 0
                        for index, dm in ser_dm.iteritems():

                            if index == length - 1:
                                break

                            chance = self._calcLogNDistribPDF(median=median, sigmaG=sigma, D2=ser_dm.loc[index + 1], D1=ser_dm.loc[index])
                            totalMass_gr += chance * serAggMass_gr.loc[index]
                            EXP_ABS_CRS_AVE += chance * exp(pABS(ser_dm.loc[index]))
                            EXP_SCA_CRS_AVE += chance * exp(pSCA(ser_dm.loc[index]))
                            sum += chance

                        if sum > 0.85:
                            arrEXPABSTotal.append(EXP_ABS_CRS_AVE)
                            arrEXPSCATotal.append(EXP_SCA_CRS_AVE)
                            arrEXPMACTotal.append(EXP_ABS_CRS_AVE * (10 ** (-12)) / totalMass_gr)
                            arrEXPMSCTotal.append(EXP_SCA_CRS_AVE * (10 ** (-12)) / totalMass_gr)
                            arrEXPSSATotal.append(EXP_SCA_CRS_AVE / (EXP_SCA_CRS_AVE + EXP_ABS_CRS_AVE))
                        else:
                            arrEXPABSTotal.append(np.nan)
                            arrEXPSCATotal.append(np.nan)
                            arrEXPMACTotal.append(np.nan)
                            arrEXPMSCTotal.append(np.nan)
                            arrEXPSSATotal.append(np.nan)

                    arrEXPABS_AVE_STD = _getAverageAndSTD(arrEXPABSTotal)
                    arrEXPSCA_AVE_STD = _getAverageAndSTD(arrEXPSCATotal)
                    arrEXPMAC_AVE_STD = _getAverageAndSTD(arrEXPMACTotal)
                    arrEXPMSC_AVE_STD = _getAverageAndSTD(arrEXPMSCTotal)
                    arrEXPSSA_AVE_STD = _getAverageAndSTD(arrEXPSSATotal)

                    sum = 0
                    totalMass_gr, RDG_ABSTotal, RDG_SCATotal, TMatrix_ABSTotal, TMatrix_SCATotal = 0, 0, 0, 0, 0

                    for index, dm in ser_dm.iteritems():

                        if index == length - 1:
                            break

                        chance = self._calcLogNDistribPDF(median=median, sigmaG=sigma, D2=ser_dm.loc[index + 1], D1=ser_dm.loc[index])
                        totalMass_gr += chance * serAggMass_gr.loc[index]
                        RDG_ABSTotal += chance * serRDG_ABS.loc[index]
                        RDG_SCATotal += chance * serRDG_SCA.loc[index]
                        TMatrix_ABSTotal += chance * serTMatrix_ABS.loc[index]
                        TMatrix_SCATotal += chance * serTMatrix_SCA.loc[index]
                        sum += chance

                    result['sigma'] = sigma
                    result['median'] = median
                    result['sum'] = sum

                    if sum >= 0.85:
                        result['fileName'] = fileName
                        result['TotalMassGram'] = totalMass_gr
                        result['ABS_RDG_Total'] = RDG_ABSTotal  # in um2
                        result['SCA_RDG_Total'] = RDG_SCATotal  # in um2
                        result['ABS_TMatrix_Total'] = TMatrix_ABSTotal  # in um2
                        result['SCA_TMatrix_Total'] = TMatrix_SCATotal  # in um2
                        result['MAC_RDG_Total'] = RDG_ABSTotal * (10 ** (-12)) / totalMass_gr  # in m2/g
                        result['MSC_RDG_Total'] = RDG_SCATotal * (10 ** (-12)) / totalMass_gr  # in m2/g
                        result['MAC_TMatrix_Total'] = TMatrix_ABSTotal * (10 ** (-12)) / totalMass_gr  # in m2/g
                        result['MSC_TMatrix_Total'] = TMatrix_SCATotal * (10 ** (-12)) / totalMass_gr  # in m2/g
                        result['SSA_RDG_Total'] = RDG_SCATotal / (RDG_SCATotal + RDG_ABSTotal)
                        result['SSA_TMatrix_Total'] = TMatrix_SCATotal / (TMatrix_SCATotal + TMatrix_ABSTotal)

                        result['EXP_ABS_AVE'] = arrEXPABS_AVE_STD[0]
                        result['EXP_SCA_AVE'] = arrEXPSCA_AVE_STD[0]
                        result['EXP_MAC_AVE'] = arrEXPMAC_AVE_STD[0]
                        result['EXP_MSC_AVE'] = arrEXPMSC_AVE_STD[0]
                        result['EXP_SSA_AVE'] = arrEXPSSA_AVE_STD[0]

                        result['EXP_ABS_STD'] = arrEXPABS_AVE_STD[1]
                        result['EXP_SCA_STD'] = arrEXPSCA_AVE_STD[1]
                        result['EXP_MAC_STD'] = arrEXPMAC_AVE_STD[1]
                        result['EXP_MSC_STD'] = arrEXPMSC_AVE_STD[1]
                        result['EXP_SSA_STD'] = arrEXPSSA_AVE_STD[1]

                    else:
                        result['fileName'] = fileName
                        result['TotalMassGram'] = np.nan
                        result['ABS_RDG_Total'] = np.nan
                        result['SCA_RDG_Total'] = np.nan
                        result['ABS_TMatrix_Total'] = np.nan
                        result['SCA_TMatrix_Total'] = np.nan
                        result['MAC_RDG_Total'] = np.nan
                        result['MSC_RDG_Total'] = np.nan
                        result['MAC_TMatrix_Total'] = np.nan
                        result['MSC_TMatrix_Total'] = np.nan
                        result['SSA_RDG_Total'] = np.nan
                        result['SSA_TMatrix_Total'] = np.nan

                        result['EXP_ABS_AVE'] = np.nan
                        result['EXP_ABS_STD'] = np.nan
                        result['EXP_SCA_AVE'] = np.nan
                        result['EXP_SCA_STD'] = np.nan
                        result['EXP_MAC_AVE'] = np.nan
                        result['EXP_MAC_STD'] = np.nan
                        result['EXP_MSC_AVE'] = np.nan
                        result['EXP_MSC_STD'] = np.nan
                        result['EXP_SSA_AVE'] = np.nan
                        result['EXP_SSA_STD'] = np.nan

                    df.loc[len(df)] = result

            self.dictTotalOpticalProp[fileName] = df

            if self.__counter == 0:
                dfInfoDB = pd.read_csv(f"{self.__folderNameGraph}\__ExperimentResults\TotalOptExp1.csv")
                dfInfoDB.loc[len(dfInfoDB)] = result
                dfInfoDB.to_csv(f"{self.__folderNameGraph}\__ExperimentResults\TotalOptExp1.csv", index=False)

            if self.__counter == 1:
                df.to_csv(f"{self.__folderNameGraph}\__ExperimentResults\TotalOptExp1.csv", index=False)
                self.__counter = 0

        except Exception as e:
            logging.exception(e)
            raise

    def CalcEfficiency(self, fileName):
        try:

            self.dictData[fileName]['ABS_RDG_Eff'] = self.dictData[fileName]['ABS_RDG'] / (((self.dictData[fileName]['dm'] / 1000) ** 2) * pi / 4)
            self.dictData[fileName]['SCA_RDG_Eff'] = self.dictData[fileName]['SCA_RDG'] / (((self.dictData[fileName]['dm'] / 1000) ** 2) * pi / 4)
            self.dictData[fileName]['ABS_TMatrix_Eff'] = self.dictData[fileName]['ABS_TMatrix'] / (((self.dictData[fileName]['dm'] / 1000) ** 2) * pi / 4)
            self.dictData[fileName]['SCA_TMatrix_Eff'] = self.dictData[fileName]['SCA_TMatrix'] / (((self.dictData[fileName]['dm'] / 1000) ** 2) * pi / 4)

        except Exception as e:
            logging.exception(e)
            raise

    def GetDataLabel(self, fileName, index):
        try:
            self.dictDataLabel[fileName] = "$D_m=$" + str(round(self.ser_DmAve.loc[index], 2)) + ", " + \
                                           "$\\rho_{eff,100}$=" + str(round(self.ser_rhoEff100nmAve.loc[index], 1)) + ", " + \
                                           "$\sigma_p|d_m=$" + str(round(self.ser_SigmaMob.loc[index], 2)) + ", " + \
                                           "$D_{TEM}=$" + str(abs(round(self.ser_D_TEM.loc[index], 2))) + ", " + \
                                           "$d_{p,100}=$" + str(round(self.ser_dp100_nano.loc[index], 1))

            self.dictDataLabelEXP[fileName] = "$D_m=$" + str(round(self.ser_DmAve.loc[index], 2)) + ", " + \
                                              "$\\rho_{eff,100}$=" + str(round(self.ser_rhoEff100nmAve.loc[index], 1)) + ", " + \
                                              "$D_{TEM}=$" + str(abs(round(self.ser_D_TEM.loc[index], 2))) + ", " + \
                                              "$d_{p,100}=$" + str(round(self.ser_dp100_nano.loc[index], 1))

        except Exception as e:
            logging.exception(e)
            raise

    def _LogNormalFit(self, x, y):
        try:

            if len(x) == len(y):
                Sum1 = 0.0
                Sum2 = 0.0
                Sum3 = 0.0
                Sum4 = 0.0
                Sum5 = 0.0
                Cols = len(x)
                ######### D_G
                for j in range(Cols):
                    Sum1 = Sum1 + (y[j] * log(x[j]))
                    Sum2 = Sum2 + y[j]
                if Sum2 != 0:
                    D_G = exp(Sum1 / Sum2)
                else:
                    D_G = 0
                ##################

                ######### Sigma_G
                for j in range(Cols):
                    Sum3 = Sum3 + (y[j] * ((log(x[j]) - log(D_G)) ** 2))
                Sigma_G = exp((Sum3 / (Sum2 - 1)) ** (0.5))
                ##################

                ######### Total Concentration in (#/cm^3)
                for j in range(Cols - 1):
                    Sum4 = Sum4 + (y[j] * (log(x[j + 1], 10) - log(x[j], 10)))
                Total_Conc = Sum4
                ##################

                ######### Median Diameter in nm
                D_Median = 0
                for j in range(Cols - 1):
                    Sum5 = Sum5 + (y[j] * (log(x[j + 1], 10) - log(x[j], 10)))
                    if Sum5 > ((Sum4 / 2)):
                        D_Median = (x[j - 1] * x[j]) ** 0.5
                        break
                result = {'D_G': D_G, 'Sigma_G': Sigma_G, 'Total_Conc': Total_Conc, 'D_Median': D_Median}
                return result

            else:
                logging.exception(f"Lognormal fit error: array lengths are not equal {len(x)},{len(y)}")
            ##################

        except Exception as e:
            logging.exception(e)
            raise

    def SMPSLogNormalFitter(self, df):
        try:
            dm = df['dm']
            conc = df['Conc']
            return self._LogNormalFit(dm, conc)

        except Exception as e:
            logging.exception(e)
            raise

    ###################################
    ###################################
    ###################################
    ###################################

    def ModeModelSelector(self, mode1, model1, mode2, model2):
        try:
            if model1 == "RDG":
                columnNameA = "_RDG"
                columnDetailA = 'RDG-FA'
            elif model1 == "Tmatrix":
                columnNameA = "_TMatrix"
                columnDetailA = 'T-matrix'
            #########################################
            if model2 == "RDG":
                columnNameB = "_RDG"
                columnDetailB = 'RDG-FA'
            elif model2 == "Tmatrix":
                columnNameB = "_TMatrix"
                columnDetailB = 'T-matrix'
            #########################################
            #########################################
            #########################################
            if mode1 == "TR_":
                columnDetailA += ', New Model'
            elif mode1 == "RE_":
                columnDetailA += ', Traditional Equivalent'
            #########################################
            if mode2 == "TR_":
                columnDetailB += ', New Model'
            elif mode2 == "RE_":
                columnDetailB += ', Traditional Equivalent'

            return columnNameA, columnDetailA, columnNameB, columnDetailB
        except Exception as e:
            logging.exception(e)
            raise

    ###################################
    ###################################
    ###################################

    def PlotTotalGraphsLineEXP(self, mode1, model1, mode2, model2, append):
        try:
            columnNameA, columnDetailA, columnNameB, columnDetailB = self.ModeModelSelector(mode1, model1, mode2, model2)

            folderName = f'{mode1}_{model1}--{mode2}_{model2}/_TotalLineGraphs' + append
            ######################## ABS Cross Section
            columnName1 = 'ABS_CRS_AVE' + columnNameA + '_MDL'
            columnName1_STD = 'ABS_CRS_STD' + columnNameA + '_MDL'
            columnDetail1 = columnDetailA
            columnName2 = 'ABS_CRS_AVE' + columnNameB + '_MDL'
            columnName2_STD = 'ABS_CRS_STD' + columnNameB + '_MDL'
            columnDetail2 = columnDetailB
            mean = 'ABS_CRS_AVE_EXP'
            STD = 'ABS_CRS_STD_EXP'
            titleName = "Absorption Cross-section"
            titleAppend = " Vs. Aerodynamic Diameter"
            self.PlotSettingTotalLinesEXP(colName1=columnName1, colName1_STD=columnName1_STD, colDetail1=columnDetail1, mode1=mode1,
                                          colName2=columnName2, colName2_STD=columnName2_STD, colDetail2=columnDetail2, mode2=mode2,
                                          titleName=titleName, titleA=titleAppend,
                                          shareY='all', folderName=folderName,
                                          yScale='log', yLabel=self.__yLabelPlotCrossSections,
                                          mean=mean, STD=STD,
                                          yAxisFormat='%1.1e')

            ######################## SCA Cross Section
            columnName1 = 'SCA_CRS_AVE' + columnNameA + '_MDL'
            columnName1_STD = 'SCA_CRS_STD' + columnNameA + '_MDL'
            columnDetail1 = columnDetailA
            columnName2 = 'SCA_CRS_AVE' + columnNameB + '_MDL'
            columnName2_STD = 'SCA_CRS_STD' + columnNameB + '_MDL'
            columnDetail2 = columnDetailB
            mean = 'SCA_CRS_AVE_EXP'
            STD = 'SCA_CRS_STD_EXP'
            titleName = "Scattering Cross-section"
            titleAppend = " Vs. Aerodynamic Diameter"
            self.PlotSettingTotalLinesEXP(colName1=columnName1, colName1_STD=columnName1_STD, colDetail1=columnDetail1, mode1=mode1,
                                          colName2=columnName2, colName2_STD=columnName2_STD, colDetail2=columnDetail2, mode2=mode2,
                                          titleName=titleName, titleA=titleAppend,
                                          shareY='all', folderName=folderName,
                                          yScale='log', yLabel=self.__yLabelPlotCrossSections,
                                          mean=mean, STD=STD,
                                          yAxisFormat='%1.1e')

            ######################## SSA
            columnName1 = 'SSA_AVE' + columnNameA + '_MDL'
            columnName1_STD = 'SSA_STD' + columnNameA + '_MDL'
            columnDetail1 = columnDetailA
            columnName2 = 'SSA_AVE' + columnNameB + '_MDL'
            columnName2_STD = 'SSA_STD' + columnNameB + '_MDL'
            columnDetail2 = columnDetailB
            mean = 'SSA_AVE_EXP'
            STD = 'SSA_STD_EXP'
            titleName = "Single-scattering Albedo (SSA)"
            titleAppend = " Vs. Aerodynamic Diameter"
            self.PlotSettingTotalLinesEXP(colName1=columnName1, colName1_STD=columnName1_STD, colDetail1=columnDetail1, mode1=mode1,
                                          colName2=columnName2, colName2_STD=columnName2_STD, colDetail2=columnDetail2, mode2=mode2,
                                          titleName=titleName, titleA=titleAppend,
                                          shareY='all', folderName=folderName,
                                          yScale='log', yLabel=self.__yLabelSSA,
                                          mean=mean, STD=STD,
                                          yAxisFormat='%1.2f')

            ######################## Absorption Eff
            columnName1 = 'ABS_EFF_AVE' + columnNameA + '_MDL'
            columnName1_STD = 'ABS_EFF_STD' + columnNameA + '_MDL'
            columnDetail1 = columnDetailA
            columnName2 = 'ABS_EFF_AVE' + columnNameB + '_MDL'
            columnName2_STD = 'ABS_EFF_STD' + columnNameB + '_MDL'
            columnDetail2 = columnDetailB
            mean = 'ABS_EFF_AVE_EXP'
            STD = 'ABS_EFF_STD_EXP'
            titleName = "Absorption Efficiency"
            titleAppend = " Vs. Aerodynamic Diameter"
            self.PlotSettingTotalLinesEXP(colName1=columnName1, colName1_STD=columnName1_STD, colDetail1=columnDetail1, mode1=mode1,
                                          colName2=columnName2, colName2_STD=columnName2_STD, colDetail2=columnDetail2, mode2=mode2,
                                          titleName=titleName, titleA=titleAppend,
                                          shareY='all', folderName=folderName,
                                          yScale='log', yLabel=self.__yLabelEFF,
                                          mean=mean, STD=STD,
                                          yAxisFormat='%1.2f')

            ######################## Scattering Eff
            columnName1 = 'SCA_EFF_AVE' + columnNameA + '_MDL'
            columnName1_STD = 'SCA_EFF_STD' + columnNameA + '_MDL'
            columnDetail1 = columnDetailA
            columnName2 = 'SCA_EFF_AVE' + columnNameB + '_MDL'
            columnName2_STD = 'SCA_EFF_STD' + columnNameB + '_MDL'
            columnDetail2 = columnDetailB
            mean = 'SCA_EFF_AVE_EXP'
            STD = 'SCA_EFF_STD_EXP'
            titleName = "Scattering Efficiency"
            titleAppend = " Vs. Aerodynamic Diameter"
            self.PlotSettingTotalLinesEXP(colName1=columnName1, colName1_STD=columnName1_STD, colDetail1=columnDetail1, mode1=mode1,
                                          colName2=columnName2, colName2_STD=columnName2_STD, colDetail2=columnDetail2, mode2=mode2,
                                          titleName=titleName, titleA=titleAppend,
                                          shareY='all', folderName=folderName,
                                          yScale='log', yLabel=self.__yLabelEFF,
                                          mean=mean, STD=STD,
                                          yAxisFormat='%1.2f')

        except Exception as e:
            logging.exception(e)
            raise

    def PlotTotalGraphsBarsEXP(self, mode1, model1, mode2, model2, append):
        try:
            columnNameA, columnDetailA, columnNameB, columnDetailB = self.ModeModelSelector(mode1, model1, mode2, model2)

            folderName = f'{mode1}_{model1}--{mode2}_{model2}/_TotalBarGraphs' + append

            ######################## MAC
            columnName1 = 'MAC' + columnNameA + '_Total'
            columnDetail1 = columnDetailA
            columnName2 = 'MAC' + columnNameB + '_Total'
            columnDetail2 = columnDetailB
            mean = self.__EXP_MAC_AVE
            STD = self.__EXP_MAC_STD
            mean1 = 'EXP_MAC_AVE'
            STD1 = 'EXP_MAC_STD'
            mean2 = 'EXP_MAC_AVE'
            STD2 = 'EXP_MAC_STD'
            titleName = "Total Mass-absorption Coefficient (MAC)"
            titleAppend = ""
            self.PlotSettingTotalBarsEXP(colName1=columnName1, colDetail1=columnDetail1, mode1=mode1,
                                         colName2=columnName2, colDetail2=columnDetail2, mode2=mode2,
                                         titleName=titleName, titleA=titleAppend,
                                         mean=mean, STD=STD,
                                         mean1=mean1, STD1=STD1,
                                         mean2=mean2, STD2=STD2,
                                         shareY='all', folderName=folderName,
                                         yScale='linear', yLabel=self.__yLabelMAC,
                                         showValue=False, valueFormat='{:.2f}', yAxisFormat='%1.2f')

            ######################## MSC
            columnName1 = 'MSC' + columnNameA + '_Total'
            columnDetail1 = columnDetailA
            columnName2 = 'MSC' + columnNameB + '_Total'
            columnDetail2 = columnDetailB
            mean = self.__EXP_MSC_AVE
            STD = self.__EXP_MSC_STD
            mean1 = 'EXP_MSC_AVE'
            STD1 = 'EXP_MSC_STD'
            mean2 = 'EXP_MSC_AVE'
            STD2 = 'EXP_MSC_STD'
            titleName = "Total Mass-scattering Coefficient (MSC)"
            titleAppend = ""
            self.PlotSettingTotalBarsEXP(colName1=columnName1, colDetail1=columnDetail1, mode1=mode1,
                                         colName2=columnName2, colDetail2=columnDetail2, mode2=mode2,
                                         titleName=titleName, titleA=titleAppend,
                                         mean=mean, STD=STD,
                                         mean1=mean1, STD1=STD1,
                                         mean2=mean2, STD2=STD2,
                                         shareY='all', folderName=folderName,
                                         yScale='linear', yLabel=self.__yLabelMSC,
                                         showValue=False, valueFormat='{:.2f}', yAxisFormat='%1.2f')

            ######################## Absorption Cross-section
            columnName1 = 'ABS' + columnNameA + '_Total'
            columnDetail1 = columnDetailA
            columnName2 = 'ABS' + columnNameB + '_Total'
            columnDetail2 = columnDetailB
            mean = 0
            STD = 0
            mean1 = 'EXP_ABS_AVE'
            STD1 = 'EXP_ABS_STD'
            mean2 = 'EXP_ABS_AVE'
            STD2 = 'EXP_ABS_STD'
            titleName = "Total Absorption Cross-section"
            titleAppend = ""
            self.PlotSettingTotalBarsEXP(colName1=columnName1, colDetail1=columnDetail1, mode1=mode1,
                                         colName2=columnName2, colDetail2=columnDetail2, mode2=mode2,
                                         titleName=titleName, titleA=titleAppend,
                                         mean=mean, STD=STD,
                                         mean1=mean1, STD1=STD1,
                                         mean2=mean2, STD2=STD2,
                                         shareY='all', folderName=folderName,
                                         yScale='log', yLabel=self.__yLabelPlotCrossSections,
                                         showValue=False, valueFormat='{:.2E}', yAxisFormat='%1.1e')

            ######################## Scattering Cross-section
            columnName1 = 'SCA' + columnNameA + '_Total'
            columnDetail1 = columnDetailA
            columnName2 = 'SCA' + columnNameB + '_Total'
            columnDetail2 = columnDetailB
            mean = 0
            STD = 0
            mean1 = 'EXP_SCA_AVE'
            STD1 = 'EXP_SCA_STD'
            mean2 = 'EXP_SCA_AVE'
            STD2 = 'EXP_SCA_STD'
            titleName = "Total Scattering Cross-section"
            titleAppend = ""
            self.PlotSettingTotalBarsEXP(colName1=columnName1, colDetail1=columnDetail1, mode1=mode1,
                                         colName2=columnName2, colDetail2=columnDetail2, mode2=mode2,
                                         titleName=titleName, titleA=titleAppend,
                                         mean=mean, STD=STD,
                                         mean1=mean1, STD1=STD1,
                                         mean2=mean2, STD2=STD2,
                                         shareY='all', folderName=folderName,
                                         yScale='log', yLabel=self.__yLabelPlotCrossSections,
                                         showValue=False, valueFormat='{:.2E}', yAxisFormat='%1.1e')

            ######################## SSA Total
            columnName1 = 'SSA' + columnNameA + '_Total'
            columnDetail1 = columnDetailA
            columnName2 = 'SSA' + columnNameB + '_Total'
            columnDetail2 = columnDetailB
            mean = self.__EXP_SSA_AVE
            STD = self.__EXP_SSA_STD
            mean1 = 'EXP_SSA_AVE'
            STD1 = 'EXP_SSA_STD'
            mean2 = 'EXP_SSA_AVE'
            STD2 = 'EXP_SSA_STD'
            titleName = "Total Single-scattering Albedo (SSA)"
            titleAppend = ""
            self.PlotSettingTotalBarsEXP(colName1=columnName1, colDetail1=columnDetail1, mode1=mode1,
                                         colName2=columnName2, colDetail2=columnDetail2, mode2=mode2,
                                         titleName=titleName, titleA=titleAppend,
                                         mean=mean, STD=STD,
                                         mean1=mean1, STD1=STD1,
                                         mean2=mean2, STD2=STD2,
                                         shareY='all', folderName=folderName,
                                         yScale='linear', yLabel=self.__yLabelSSA,
                                         showValue=False, valueFormat='{:.2f}', yAxisFormat='%1.2f')

        except Exception as e:
            logging.exception(e)
            raise

    ###################################
    ###################################
    ###################################
    def PlotMultipleBarsTotalEXP(self, A1, A1Name, A2, A2Name,
                                 column1, column1T, mode1,
                                 column2, column2T, mode2,
                                 title, titleAppend,
                                 mean, STD,
                                 mean1, STD1,
                                 mean2, STD2,
                                 folderName, yCommonLabel,
                                 shareY, yAxisFormat,
                                 sigma, median,
                                 showValue=True, valueFormat='{:.2E}', yScale='log'):
        try:

            fig, ax1 = plt.subplots(nrows=len(A1.unique()), ncols=len(A2.unique()), sharex=True, sharey=shareY, figsize=(12, 12), constrained_layout=True)

            def autoLabel(self, rects, row, col, format='{:.2E}', xPos='center'):

                if 'E' in format:
                    fontS = self.__valueFontSizeTotal
                if 'f' in format:
                    fontS = self.__valueFontSizeTotal * 1.5
                ha = {'center': 'center', 'right': 'left', 'left': 'right'}
                offset = {'center': 0, 'right': 1, 'left': -1}

                for rect in rects:
                    height = rect.get_height()
                    ax1[row, col].annotate(format.format(height),
                                           xy=(rect.get_x() + rect.get_width() / 2, height),
                                           xytext=(offset[xPos] * 1, 1),  # use 3 points offset
                                           textcoords="offset points",  # in both directions
                                           ha=ha[xPos], va='bottom',
                                           fontsize=fontS)

            rowCount = 0
            for i in A1.unique():
                colCount = 0
                for j in A2.unique():
                    fileNames = self.dfMainInfo[(A1Name == i) & (A2Name == j)].AA_FileName

                    indexBarPlot = np.arange(3)

                    C1Height = []
                    C2Height = []

                    C3Height = []
                    C4Height = []

                    C3ErrHeight = []
                    C4ErrHeight = []

                    EXPHeight = []
                    EXPErr = []

                    EXPHeightShadeUp = []
                    EXPHeightShadeDown = []

                    for index, fileN in fileNames.iteritems():

                        df1 = self.dictTotalOpticalProp[mode1 + fileN]
                        selectedRow1 = df1.loc[(df1['sigma'] == sigma) & (df1['median'] == median)]

                        df2 = self.dictTotalOpticalProp[mode2 + fileN]
                        selectedRow2 = df2.loc[(df2['sigma'] == sigma) & (df2['median'] == median)]

                        if (len(selectedRow1[column1]) != 1) or (len(selectedRow2[column2]) != 1):
                            raise

                        else:
                            for c1 in selectedRow1[column1]:
                                C1Height.append(c1)
                            for c2 in selectedRow2[column2]:
                                C2Height.append(c2)

                            for c3 in selectedRow1[mean1]:
                                C3Height.append(c3)
                            for c3err in selectedRow1[STD1]:
                                C3ErrHeight.append(c3err * 2)

                            for c4 in selectedRow2[mean2]:
                                C4Height.append(c4)
                            for c4err in selectedRow2[STD2]:
                                C4ErrHeight.append(c4err * 2)

                            EXPHeight.append(mean)
                            EXPErr.append(STD * 2)

                            EXPHeightShadeUp.append(mean + 2 * STD)
                            EXPHeightShadeDown.append(mean - 2 * STD)

                    AXES_A = ax1[rowCount, colCount].bar(indexBarPlot - self.__barWidth / 2, C1Height, color=self.__A1Color, width=self.__barWidth, edgecolor='white', hatch=self.__patterns[0],
                                                         label=column1T, alpha=0.8)
                    AXES_B = ax1[rowCount, colCount].bar(indexBarPlot + self.__barWidth / 2, C2Height, color=self.__A2Color, width=self.__barWidth, edgecolor='white', hatch=self.__patterns[4],
                                                         label=column2T, alpha=0.8)

                    # EXP Integrated
                    ax1[rowCount, colCount].errorbar(indexBarPlot - self.__barWidth / 2, C3Height, yerr=C3ErrHeight, label="Experiment Integrated, AAC-PAX-CPC", capsize=6, elinewidth=3, capthick=4,
                                                     ecolor='green',
                                                     marker='o', mfc='green',
                                                     mec='gold', ms=8, mew=2, fmt='o')
                    ax1[rowCount, colCount].errorbar(indexBarPlot + self.__barWidth / 2, C4Height, yerr=C4ErrHeight, label="Experiment Integrated, AAC-PAX-CPC", capsize=6, elinewidth=3, capthick=4,
                                                     ecolor='tomato',
                                                     marker='o', mfc='tomato',
                                                     mec='gold', ms=8, mew=2, fmt='o')
                    # AXES_A1 = ax1[rowCount, colCount].bar(indexBarPlot - self.__barWidth / 2, C3Height, yerr=C3ErrHeight, capsize=3, error_kw=dict(elinewidth=2, capthick=2, marker='s'),label="integrated", ecolor='black', width=0.0,)
                    # AXES_B1 = ax1[rowCount, colCount].bar(indexBarPlot + self.__barWidth / 2, C4Height, yerr=C4ErrHeight, capsize=3, error_kw=dict(elinewidth=2, capthick=2, marker='s'),label="integrated", ecolor='black', width=0.0,)

                    # EXP
                    if mean != 0:
                        ax1[rowCount, colCount].fill_between([-1.1 * self.__barWidth, 1, 2 + self.__barWidth * 1.1], EXPHeightShadeDown, EXPHeightShadeUp, label="Experiment, PAX-Filter",
                                                             facecolor='gray', alpha=0.45)

                        # AXES_C = ax1[rowCount, colCount].bar(indexBarPlot, EXPHeight, yerr=EXPErr, capsize=4, error_kw=dict(elinewidth=3, capthick=1), ecolor='gray', width=0.0, edgecolor='white',
                        # label='Experiment')

                    if showValue:
                        autoLabel(self, rects=AXES_A, row=rowCount, col=colCount, format=valueFormat)
                        autoLabel(self, rects=AXES_B, row=rowCount, col=colCount, format=valueFormat)

                    ax1[rowCount, colCount].set_xlim(-1.1 * self.__barWidth, 2 + self.__barWidth * 1.1)
                    ax1[rowCount, colCount].set_xticks(indexBarPlot)
                    ax1[rowCount, colCount].set_xticklabels(("$D_m=$" + str(2.2), "$D_m=$" + str(2.5), "$D_m=$" + str(2.8)), fontsize=self.__xLabelGroupFontSizeTotal)
                    ax1[rowCount, colCount].grid(True, which='both', axis="y", alpha=0.5)
                    ax1[rowCount, colCount].set_yscale(yScale)
                    ax1[rowCount, colCount].yaxis.set_tick_params(labelsize=self.__yMajorTickLabelSubplotFontSize)
                    ax1[rowCount, colCount].yaxis.set_major_formatter(FormatStrFormatter(yAxisFormat))

                    if colCount == (len(A2.unique()) - 1):
                        ax1[rowCount, colCount].set_ylabel("$\\rho_{eff,100}$= " + str(i) + " kg/m$^3$", fontsize=self.__yLabelEachSubplotFontSize, rotation=self.__yLabelEachSubplotRotation)
                        ax1[rowCount, colCount].yaxis.set_label_position("right")
                        ax1[rowCount, colCount].yaxis.labelpad = self.__yLabelEachSubplotPad
                    if rowCount == 0:
                        ax1[rowCount, colCount].set_xlabel("$\sigma_p|d_m=$" + str(j), fontsize=self.__xLabelEachSubplotFontSize)
                        ax1[rowCount, colCount].xaxis.set_label_position("top")
                        ax1[rowCount, colCount].xaxis.labelpad = self.__xLabelEachSubplotPad

                    colCount += 1
                rowCount += 1

            F1 = f"{title} for median={round(median)}-sigma={round(sigma, 2)}"
            T1 = title + titleAppend + " for " + "$d_{m,g}=$" + str(round(median)) + " nm " + "and " + "$\sigma_g=$" + str(round(sigma, 2))

            fig.subplots_adjust(top=self.__subplotGridSetup[0], wspace=self.__subplotGridSetup[1], hspace=self.__subplotGridSetup[2])
            fig.suptitle(T1, fontsize=self.__plotTitleFontSizeTotal)

            fig.text(0.5, self.__xLabelCommonSubplotAdjustment, self.__xLabelMassMob, ha='center', fontsize=self.__xLabelCommonSubplotFontSize)

            if 'e' in yAxisFormat:
                fig.text(self.__yLabelCommonSubplotAdjustment[0], 0.5, yCommonLabel, va='center', rotation='vertical', fontsize=self.__yLabelCommonSubplotFontSize)
            elif 'f' in yAxisFormat:
                fig.text(self.__yLabelCommonSubplotAdjustment[1], 0.5, yCommonLabel, va='center', rotation='vertical', fontsize=self.__yLabelCommonSubplotFontSize)

            leg = plt.legend(bbox_to_anchor=(self.__legendSubplotAdjustment[0], self.__legendSubplotAdjustment[1]), markerscale=self.__legendSubplotMarkerScale * 0.6,
                             fontsize=self.__legendSubplotFontSize,
                             loc='lower left', borderaxespad=0.)

            ################# set the lineWidth of each legend object
            for legObj in leg.legendHandles:
                legObj.set_linewidth(self.__legendSubplotLineWidth)
            ########################################### Main
            buf = io.BytesIO()
            pickle.dump(fig, buf)
            buf.seek(0)
            #############################
            self.SaveAndClosePlot(folderName=folderName, F1=F1)
            ##########################
            fig2 = pickle.load(buf)
            tokeep = [4]
            for i, ax2 in enumerate(fig2.axes):
                if not i in (tokeep):
                    fig2.delaxes(ax2)
                else:
                    axestokeep = ax2
            axestokeep.change_geometry(1, 1, 1)
            ##########################
            # axestokeep.set_xticks(indexBarPlot)
            # axestokeep.set_xticklabels(("$D_m=$" + str(2.2), "$D_m=$" + str(2.5), "$D_m=$" + str(2.8)))
            axestokeep.tick_params(axis='x', which='both', labelbottom=True, labelsize=self.__xLabelGroupFontSizeTotal * 1.2)
            axestokeep.tick_params(axis='y', which='both', labelleft=True, labelsize=self.__yMajorTickLabelSubplotFontSize * 1.2)
            leg = plt.legend(bbox_to_anchor=(self.__legendSubplotAdjustmentMain[0], self.__legendSubplotAdjustmentMain[1]),
                             markerscale=self.__legendSubplotMarkerScale * 0.6,
                             fontsize=self.__legendSubplotFontSize,
                             loc='lower left', borderaxespad=0.)
            ################# set the lineWidth of each legend object
            for legObj in leg.legendHandles:
                legObj.set_linewidth(self.__legendSubplotLineWidth)
            # fig2.set_size_inches(10, 10)
            self.SaveAndClosePlotFig(fig=fig2, folderName=folderName, F1=F1 + " - Main")


        except Exception as e:
            logging.exception(e)
            raise

    def PlotSettingTotalBarsEXP(self, colName1, colDetail1, mode1,
                                colName2, colDetail2, mode2,
                                titleName, titleA, shareY, folderName,
                                mean, STD,
                                mean1, STD1,
                                mean2, STD2,
                                yScale, yLabel, showValue, valueFormat, yAxisFormat):
        try:
            for sigma in self.__dmSigmaLogN:
                for median in self.__dmMedianLogN:
                    self.PlotMultipleBarsTotalEXP(A1=self.ser_rhoEff100nmAve, A1Name=self.dfMainInfo.AGG_EFF_RHO_100NM_CENTER,
                                                  A2=self.ser_SigmaMob, A2Name=self.dfMainInfo.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER,
                                                  column1=colName1, column1T=colDetail1, mode1=mode1,
                                                  column2=colName2, column2T=colDetail2, mode2=mode2,
                                                  title=titleName, titleAppend=titleA,
                                                  mean=mean, STD=STD,
                                                  mean1=mean1, STD1=STD1,
                                                  mean2=mean2, STD2=STD2,
                                                  folderName=folderName,
                                                  yScale=yScale, yCommonLabel=yLabel, yAxisFormat=yAxisFormat,
                                                  shareY=shareY,
                                                  sigma=sigma, median=median,
                                                  showValue=showValue, valueFormat=valueFormat)

        except Exception as e:
            logging.exception(e)
            raise

        ###################################
        ###################################
        ###################################
        ###################################

    ###################################
    ###################################
    ###################################

    def PlotMultipleLinesTotalEXP(self, A1, A1Name, A2, A2Name,
                                  column1, colName1_STD, column1T, mode1,
                                  column2, colName2_STD, column2T, mode2,
                                  title, titleAppend,
                                  mean, STD,
                                  folderName, yCommonLabel,
                                  shareY, yAxisFormat, yScale='log'):
        try:

            fig, ax1 = plt.subplots(nrows=len(A1.unique()), ncols=len(A2.unique()), sharey=shareY, sharex=True, figsize=(12, 12), constrained_layout=True)
            #############################################
            Dmlist = ["$D_m=$" + str(2.2), "$D_m=$" + str(2.5), "$D_m=$" + str(2.8)]
            #############################################
            rowCount = 0
            for i in A1.unique():
                colCount = 0
                for j in A2.unique():
                    fileNames = self.dfMainInfo[(A1Name == i) & (A2Name == j)].AA_FileName
                    c_lineStyle = 0
                    c_markerStyle = 0
                    c_lineWidth = 0
                    DmCount = 0
                    for index, fileN in fileNames.iteritems():
                        df1 = self.modelConvertedDF[mode1 + fileN]
                        # da_STD = df1['da_STD'] * 2
                        # y_STD = df1[colName1_STD] * 2
                        ax1[rowCount, colCount].plot(df1['da_AVE'], df1[column1], label=column1T + ', ' + Dmlist[DmCount],
                                                     color=self.get_color(self.__A1Color, 'yellow', portion=(DmCount + 1) * 0.1), linewidth=self.__A1LineWidth[c_lineStyle],
                                                     linestyle=self.__subplotLineStyle[c_lineStyle], alpha=self.__A1AlphaMainLine)
                        # ax1[rowCount, colCount].errorbar(df1['da_AVE'], df1[column1], xerr=da_STD, yerr=y_STD, label=column1T + ', ' + Dmlist[DmCount],
                        #                                 color=self.get_color(self.__A1Color, 'yellow', portion=(DmCount + 1) * 0.1), linewidth=self.__A1LineWidth[c_lineStyle],
                        #                                 linestyle=self.__subplotLineStyle[c_lineStyle], alpha=self.__A1AlphaMainLine)
                        # marker=self.__subplotMarkerStyle[c_markerStyle], markersize=3, markevery=4, errorevery=4,
                        # capsize=1, elinewidth=1, markeredgewidth=1)

                        df2 = self.modelConvertedDF[mode2 + fileN]
                        # da_STD = df2['da_STD'] * 2
                        # y_STD = df2[colName2_STD] * 2
                        ax1[rowCount, colCount].plot(df2['da_AVE'], df2[column2], label=column2T + ', ' + Dmlist[DmCount],
                                                     color=self.get_color(self.__A2Color, 'gray', portion=(DmCount + 1) * 0.2), linewidth=self.__A2LineWidth[c_lineStyle],
                                                     linestyle=self.__subplotLineStyle[c_lineStyle], alpha=self.__A2AlphaMainLine, )
                        # ax1[rowCount, colCount].errorbar(df2['da_AVE'], df2[column2], xerr=da_STD, yerr=y_STD, label=column2T + ', ' + Dmlist[DmCount],
                        #                                 color=self.get_color(self.__A2Color, 'gray', portion=(DmCount + 1) * 0.2), linewidth=self.__A2LineWidth[c_lineStyle],
                        #                                 linestyle=self.__subplotLineStyle[c_lineStyle], alpha=self.__A2AlphaMainLine, )
                        # errorevery=5,capsize=1, elinewidth=1, markeredgewidth=1)
                        c_lineStyle += 1
                        c_markerStyle += 1
                        c_lineWidth += 1
                        DmCount += 1

                    p8 = np.poly1d(np.polyfit(self.dfExperiment['da_AVE_EXP'], np.log(self.dfExperiment[mean].values), 8))
                    fittedEXP = []
                    for da in df2['da_AVE']:
                        fittedEXP.append(exp(p8(da)))

                    ax1[rowCount, colCount].plot(df2['da_AVE'], fittedEXP, color='black')

                    da_STD = self.dfExperiment['da_STD_EXP'] * 2
                    y_STD = self.dfExperiment[STD] * 2
                    ax1[rowCount, colCount].errorbar(self.dfExperiment['da_AVE_EXP'], self.dfExperiment[mean], xerr=da_STD, yerr=y_STD, label="Experiment, AAC-PAX-CPC",
                                                     color=self.get_color('gray', 'black'),
                                                     # linestyle='-',linewidth=3,
                                                     alpha=self.__A3AlphaMainLine, marker="D", markersize=self.__markerSize, ecolor='black',
                                                     capsize=3, elinewidth=1, markeredgewidth=1)

                    ax1[rowCount, colCount].grid(True, which='both', axis="both", alpha=0.5)
                    ax1[rowCount, colCount].set_xscale('log')
                    ax1[rowCount, colCount].set_yscale(yScale)
                    ax1[rowCount, colCount].xaxis.set_major_formatter(FormatStrFormatter("%i"))
                    ax1[rowCount, colCount].xaxis.set_major_locator(plt.MaxNLocator(6))
                    ax1[rowCount, colCount].tick_params(axis='x', which='both', labelsize=self.__xMajorTickLabelSubplotFontSize, labelrotation=self.__xTickLabelSubplotRotation)
                    ax1[rowCount, colCount].tick_params(axis='y', which='both', labelsize=self.__yMajorTickLabelSubplotFontSize)
                    ax1[rowCount, colCount].yaxis.set_major_formatter(FormatStrFormatter(yAxisFormat))
                    ax1[rowCount, colCount].set_xlim(self.__xAxisLimitsComp[0], self.__xAxisLimitsComp[1])

                    if colCount == (len(A2.unique()) - 1):
                        ax1[rowCount, colCount].set_ylabel("$\\rho_{eff,100}$= " + str(i) + " kg/m$^3$", fontsize=self.__yLabelEachSubplotFontSize, rotation=self.__yLabelEachSubplotRotation)
                        ax1[rowCount, colCount].yaxis.set_label_position("right")
                        ax1[rowCount, colCount].yaxis.labelpad = self.__yLabelEachSubplotPad

                    if rowCount == 0:
                        ax1[rowCount, colCount].set_xlabel("$\sigma_p|d_m=$" + str(j), fontsize=self.__xLabelEachSubplotFontSize)
                        ax1[rowCount, colCount].xaxis.set_label_position("top")
                        ax1[rowCount, colCount].xaxis.labelpad = self.__xLabelEachSubplotPad

                    colCount += 1
                rowCount += 1

            # rowCount = 0
            # for i in A1.unique():
            #     colCount = 0
            #     for j in A2.unique():
            #         ax1[rowCount, colCount].set_ylim(self.__xAxisLimitsComp[0], self.__xAxisLimitsComp[1])
            #         colCount += 1
            #     rowCount += 1

            F1 = f"{title}"
            T1 = title + titleAppend

            fig.subplots_adjust(top=self.__subplotGridSetup[0], wspace=self.__subplotGridSetup[1], hspace=self.__subplotGridSetup[2])
            fig.suptitle(T1, fontsize=self.__plotTitleFontSizeTotal)

            fig.text(0.5, self.__xLabelCommonSubplotAdjustment, self.__xLabelAeroDiameter, ha='center', fontsize=self.__xLabelCommonSubplotFontSize)

            if 'e' in yAxisFormat:
                fig.text(self.__yLabelCommonSubplotAdjustment[0], 0.5, yCommonLabel, va='center', rotation='vertical', fontsize=self.__yLabelCommonSubplotFontSize)
            elif 'f' in yAxisFormat:
                fig.text(self.__yLabelCommonSubplotAdjustment[1], 0.5, yCommonLabel, va='center', rotation='vertical', fontsize=self.__yLabelCommonSubplotFontSize)

            leg = plt.legend(bbox_to_anchor=(self.__legendSubplotAdjustment[0], self.__legendSubplotAdjustment[1]), markerscale=self.__legendSubplotMarkerScale, fontsize=self.__legendSubplotFontSize,
                             loc='lower left', borderaxespad=0.)
            ################# set the lineWidth of each legend object
            for legObj in leg.legendHandles:
                if hasattr(legObj, '_linewidth'):
                    w = legObj._linewidth
                    legObj.set_linewidth(w * self.__legendSubplotLineScale)
                else:
                    legObj.set_linewidth(self.__legendSubplotLineWidth)
            ########################################### Main
            buf = io.BytesIO()
            pickle.dump(fig, buf)
            buf.seek(0)
            #############################
            self.SaveAndClosePlot(folderName=folderName, F1=F1)
            ##########################
            fig2 = pickle.load(buf)
            tokeep = [4]
            for i, ax2 in enumerate(fig2.axes):
                if not i in (tokeep):
                    fig2.delaxes(ax2)
                else:
                    axestokeep = ax2
            axestokeep.change_geometry(1, 1, 1)
            ##########################
            axestokeep.tick_params(axis='x', which='both', labelbottom=True, labelsize=self.__xMajorTickLabelSubplotFontSize, labelrotation=self.__xTickLabelSubplotRotation)
            axestokeep.tick_params(axis='y', which='both', labelleft=True, labelsize=self.__yMajorTickLabelSubplotFontSize)
            leg = plt.legend(bbox_to_anchor=(self.__legendSubplotAdjustmentMain[0], self.__legendSubplotAdjustmentMain[1]),
                             markerscale=self.__legendSubplotMarkerScale,
                             fontsize=self.__legendSubplotFontSize,
                             loc='lower left', borderaxespad=0.)
            ################# set the lineWidth of each legend object
            for legObj in leg.legendHandles:
                if hasattr(legObj, '_linewidth'):
                    w = legObj._linewidth
                    legObj.set_linewidth(w * self.__legendSubplotLineScale)
                else:
                    legObj.set_linewidth(self.__legendSubplotLineWidth)
            # fig2.set_size_inches(6, 6)
            self.SaveAndClosePlotFig(fig=fig2, folderName=folderName, F1=F1 + " - Main")

        except Exception as e:
            logging.exception(e)
            raise

    def PlotSettingTotalLinesEXP(self, colName1, colName1_STD, colDetail1, mode1,
                                 colName2, colName2_STD, colDetail2, mode2,
                                 titleName, titleA, shareY, folderName,
                                 mean, STD,
                                 yScale, yLabel, yAxisFormat):
        try:

            self.PlotMultipleLinesTotalEXP(A1=self.ser_rhoEff100nmAve, A1Name=self.dfMainInfo.AGG_EFF_RHO_100NM_CENTER,
                                           A2=self.ser_SigmaMob, A2Name=self.dfMainInfo.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER,
                                           column1=colName1, colName1_STD=colName1_STD, column1T=colDetail1, mode1=mode1,
                                           column2=colName2, colName2_STD=colName2_STD, column2T=colDetail2, mode2=mode2,
                                           title=titleName, titleAppend=titleA,
                                           mean=mean, STD=STD,
                                           folderName=folderName,
                                           yScale=yScale, yCommonLabel=yLabel, yAxisFormat=yAxisFormat,
                                           shareY=shareY, )

        except Exception as e:
            logging.exception(e)
            raise

        ###################################
        ###################################
        ###################################
        ###################################

    ###################################
    ###################################
    ###################################
    ###################################

    def SaveAndClosePlot(self, folderName, F1):
        try:
            if self.__OnlySVGSaving:

                Address = GF.GetAddressTo(folderName=self.__folderNameGraph + f"\{folderName}", fileName=F1, extension="svg")
                plt.savefig(Address, format='svg', dpi=self.__figureDPI, bbox_inches='tight')

            else:
                Address = GF.GetAddressTo(folderName=self.__folderNameGraph + f"\{folderName}", fileName=F1, extension="jpg")
                plt.savefig(Address, format='jpg', dpi=self.__figureDPI, bbox_inches='tight')

                if self.__WMF_SVGSaving:
                    Address = GF.GetAddressTo(folderName=self.__folderNameGraph + f"\{folderName}", fileName=F1, extension="svg")
                    plt.savefig(Address, format='svg', dpi=self.__figureDPI, bbox_inches='tight')

                    AddressWMF = GF.GetAddressTo(folderName=self.__folderNameGraph + f"\{folderName}", fileName=F1, extension="wmf")
                    subprocess.call([self.__inkScapePath, str(Address.resolve()), '--export-wmf', str(AddressWMF.resolve())])

            plt.clf()
            plt.close()

        except Exception as e:
            logging.exception(e)
            raise

    def SaveAndClosePlotFig(self, fig, folderName, F1):
        try:
            if self.__OnlySVGSaving:

                Address = GF.GetAddressTo(folderName=self.__folderNameGraph + f"\{folderName}", fileName=F1, extension="svg")
                fig.savefig(Address, format='svg', dpi=self.__figureDPI, bbox_inches='tight')

            else:
                Address = GF.GetAddressTo(folderName=self.__folderNameGraph + f"\{folderName}", fileName=F1, extension="jpg")
                fig.savefig(Address, format='jpg', dpi=self.__figureDPI, bbox_inches='tight')

                if self.__WMF_SVGSaving:
                    Address = GF.GetAddressTo(folderName=self.__folderNameGraph + f"\{folderName}", fileName=F1, extension="svg")
                    fig.savefig(Address, format='svg', dpi=self.__figureDPI, bbox_inches='tight')

                    AddressWMF = GF.GetAddressTo(folderName=self.__folderNameGraph + f"\{folderName}", fileName=F1, extension="wmf")
                    subprocess.call([self.__inkScapePath, str(Address.resolve()), '--export-wmf', str(AddressWMF.resolve())])

            fig.clf()

        except Exception as e:
            logging.exception(e)
            raise

    def get_color(self, color1, color2, portion=1):
        try:
            colorRGBA1 = colors.to_rgba(color1)
            colorRGBA2 = colors.to_rgba(color2)
            alpha = (colorRGBA1[3] + colorRGBA2[3] * portion) / 2
            red = ((colorRGBA1[0] * colorRGBA1[3]) + (colorRGBA2[0] * colorRGBA2[3] * portion)) / (colorRGBA1[3] + colorRGBA2[3] * portion)
            green = ((colorRGBA1[1] * colorRGBA1[3]) + (colorRGBA2[1] * colorRGBA2[3] * portion)) / (colorRGBA1[3] + colorRGBA2[3] * portion)
            blue = ((colorRGBA1[2] * colorRGBA1[3]) + (colorRGBA2[2] * colorRGBA2[3] * portion)) / (colorRGBA1[3] + colorRGBA2[3] * portion)
            return (red, green, blue, alpha)
        except Exception as e:
            logging.exception(e)
            raise

    ######################################################################
    ######################################################################
    ######################################################################
    ######################################################################
