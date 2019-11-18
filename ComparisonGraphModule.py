from ConfigReaderModule import logging
import GeneralFunctions as GF
##############
import numpy as np
from matplotlib import rcParams
from matplotlib import pyplot as plt
import pandas as pd
from matplotlib.ticker import FormatStrFormatter
from math import pi
from math import exp
from math import isnan
from math import log
import subprocess

from matplotlib import colors

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
            self.__dmSigmaLogN = [1.867]
            self.__dmMedianLogN = [133.3]
            #################################################
            #################################################
            ################################################# Comparison Graphs
            self.__CRSGraphs = True
            self.__ErrorAndRatioGraphs = True
            self.__MACMSCGraphs = True
            self.__EfficiencyGraphs = True
            self.__SSAGraphs = True
            self.__dpInfoGraphs = True
            self.__dmDistributionGraphsBar = True
            self.__dmDistributionGraphsLine = True
            self.__dmDistributionRatioGraphsBar = True
            ###############
            self.__DmCTE = True
            self.__rhoEff100nmCTE = True
            self.__sigmaEachMobCTE = True
            #################################################
            self.__TotalDFLabels = ['sigma', 'median', 'TotalMassGram',
                                    'ABS_RDG_Total', 'SCA_RDG_Total',
                                    'ABS_TMatrix_Total', 'SCA_TMatrix_Total',
                                    'MAC_RDG_Total', 'MSC_RDG_Total',
                                    'MAC_TMatrix_Total', 'MSC_TMatrix_Total',
                                    'SSA_RDG_Total', 'SSA_TMatrix_Total', 'sum']
            ################################################# General Variables
            self.__folderNameGraph = folderInfo['FOLDER_NAME_GRAPH']
            self.__folderNameData = folderInfo['FOLDER_NAME_DATA']
            self.__WMF_SVGSaving = False
            self.__OnlySVGSaving = True
            self.__inkScapePath = "C://Program Files//inkscape//inkscape.exe"
            ################################################# Graphs Labels
            self.__xLabelMobDiameter = "Mobility-equivalent Diameter (nm)"
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
            self.__A1LineWidth = [3, 1.9, 1.3]
            self.__A1AlphaMainLine = 0.65
            self.__A2Color = 'blue'
            self.__A2LineWidth = [3, 1.9, 1.3]
            self.__A2AlphaMainLine = 0.7
            self.__subplotLineStyle = ['-', '--', '-.']
            self.__subplotMarkerStyle = ["o", "X", "^"]
            self.__subplotMarkerSize = 8
            self.__xMajorTickLabelSubplotFontSize = 18
            self.__xTickLabelSubplotRotation = 45
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
            self.__yLabelCommonSubplotAdjustment = [0.018, 0.044]  # lower to get left
            self.__legendSubplotAdjustment = [1.25, 1.7]
            self.__legendSubplotMarkerScale = 2.25
            self.__legendSubplotFontSize = 20
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
            self.__xAxisLimitsComp = [49, 960]
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

    def GetDataDefinition(self):
        try:
            self.serMainFileName = self.dfMainInfo.loc[self.dfMainInfo['AA_Plot'] == 1, 'AA_FileName']
            self.ser_Dm = self.dfMainInfo.loc[self.dfMainInfo['AA_Plot'] == 1, 'AGG_EFF_DM_CENTER']
            self.ser_rhoEff100nm = self.dfMainInfo.loc[self.dfMainInfo['AA_Plot'] == 1, 'AGG_EFF_RHO_100NM_CENTER']
            self.ser_SigmaMob = self.dfMainInfo.loc[self.dfMainInfo['AA_Plot'] == 1, 'AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER']
            self.ser_D_TEM = self.dfMainInfo.loc[self.dfMainInfo['AA_Plot'] == 1, 'D_TEM']
            self.ser_dp100_nano = self.dfMainInfo.loc[self.dfMainInfo['AA_Plot'] == 1, 'dp100_nano']

        except Exception as e:
            logging.exception(e)
            raise

    def GetDataLabel(self, fileName, index):
        try:
            self.dictDataLabel[fileName] = "$D_m=$" + str(round(self.ser_Dm.loc[index], 2)) + ", " + \
                                           "$\\rho_{eff,100}$=" + str(round(self.ser_rhoEff100nm.loc[index], 1)) + ", " + \
                                           "$\sigma_p|d_m=$" + str(round(self.ser_SigmaMob.loc[index], 2)) + ", " + \
                                           "$D_{TEM}=$" + str(abs(round(self.ser_D_TEM.loc[index], 2))) + ", " + \
                                           "$d_{p,100}=$" + str(round(self.ser_dp100_nano.loc[index], 1))

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

    def CalcTotalMACMSC(self, dataFrame, fileName):
        try:

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

                    else:
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

                    df.loc[len(df)] = result

            self.dictTotalOpticalProp[fileName] = df

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

    def RDG_TMatrixComparisonGraphs(self):
        try:

            self.dfMainInfo = pd.read_csv(GF.GetAddressTo(folderName=self.__folderNameData, fileName='beacon', extension='csv'))
            self.GetDataDefinition()
            self.dictData = {}
            self.dictDataLabel = {}
            self.dictTotalOpticalProp = {}

            for index, fileN in self.serMainFileName.iteritems():
                self.CoreReader("TR_" + fileN, index)
                self.CoreReader("RE_" + fileN, index)

            #### compring RE and TR
            self.CoreComparer("TR_", "RDG", "TR_", "Tmatrix")
            self.CoreComparer("RE_", "RDG", "RE_", "Tmatrix")

            self.CoreComparer("TR_", "RDG", "RE_", "RDG")
            self.CoreComparer("TR_", "Tmatrix", "RE_", "Tmatrix")

        except Exception as e:
            logging.exception(e)
            raise

    def CoreReader(self, fileN, index):
        try:
            self.dictData[fileN] = pd.read_csv(GF.GetAddressTo(folderName=self.__folderNameData, fileName=fileN))
            self.dictData[fileN] = self.dictData[fileN].replace(0, np.nan)

            if self.__dpDeviation:
                self.Check_dpRatio(fileName=fileN)
            self.CalcTotalMACMSC(self.dictData[fileN], fileN)
            self.CalcEfficiency(fileName=fileN)
            self.GetDataLabel(fileName=fileN, index=index)

        except Exception as e:
            logging.exception(e)
            raise

    def CoreComparer(self, mode1, model1, mode2, model2):
        try:
            # if self.__CRSGraphs:
            #   self.PlotCrossSections()
            # if self.__ErrorAndRatioGraphs:
            #    self.PlotErrorAndRatio()
            # if self.__SSAGraphs:
            #   self.PlotSSA()
            # if self.__MACMSCGraphs:
            #    self.PlotMACMSC()
            # if self.__EfficiencyGraphs:
            #   self.PlotEff()
            # if self.__dpInfoGraphs:
            #    self.PlotdpInfo()
            if self.__dmDistributionGraphsLine:
                self.PlotTotalGraphsLine(mode1, model1, mode2, model2)
            if self.__dmDistributionGraphsBar:
                self.PlotTotalGraphsBar(mode1, model1, mode2, model2)

            # if self.__dmDistributionRatioGraphsBar:
            #   self.PlotTotalRatioGraphsBar()

        except Exception as e:
            logging.exception(e)
            raise

    ###################################
    ###################################
    ###################################
    ###################################

    def PlotTotalRatioGraphsBar(self):
        try:
            folderName = '_RatioBarGraphs'
            ######################## ABS Cross Section
            columnName1 = 'ABS_RDG_Total'
            columnDetail1 = 'RDG-FA'
            columnName2 = 'ABS_TMatrix_Total'
            columnDetail2 = 'T-matrix'
            titleName = "T-matrix and RDG-FA Absorption Cross-section Ratio"
            titleAppend = ""
            self.PlotSettingTotalRatio(colName1=columnName1, colDetail1=columnDetail1,
                                       colName2=columnName2, colDetail2=columnDetail2,
                                       titleName=titleName, titleA=titleAppend,
                                       shareY='all', folderName=folderName,
                                       yScale='linear', yLabel=self.__yLabelRatioCRS,
                                       showValue=True, valueFormat='{:.2f}', yAxisFormat='%1.2f')
            ######################## SCA Cross Section
            columnName1 = 'SCA_RDG_Total'
            columnDetail1 = 'RDG-FA'
            columnName2 = 'SCA_TMatrix_Total'
            columnDetail2 = 'T-matrix'
            titleName = "T-matrix and RDG-FA Scattering Cross-section Ratio"
            titleAppend = ""
            self.PlotSettingTotalRatio(colName1=columnName1, colDetail1=columnDetail1,
                                       colName2=columnName2, colDetail2=columnDetail2,
                                       titleName=titleName, titleA=titleAppend,
                                       shareY='all', folderName=folderName,
                                       yScale='linear', yLabel=self.__yLabelRatioCRS,
                                       showValue=True, valueFormat='{:.2f}', yAxisFormat='%1.2f')
            ######################## MAC
            columnName1 = 'MAC_RDG_Total'
            columnDetail1 = 'RDG-FA'
            columnName2 = 'MAC_TMatrix_Total'
            columnDetail2 = 'T-matrix'
            titleName = "T-matrix and RDG-FA Mass-absorption Coefficient (MAC) Ratio"
            titleAppend = ""
            self.PlotSettingTotalRatio(colName1=columnName1, colDetail1=columnDetail1,
                                       colName2=columnName2, colDetail2=columnDetail2,
                                       titleName=titleName, titleA=titleAppend,
                                       shareY='all', folderName=folderName,
                                       yScale='linear', yLabel=self.__yLabelRatioMAC,
                                       showValue=True, valueFormat='{:.2f}', yAxisFormat='%1.2f')
            ######################## MSC
            columnName1 = 'MSC_RDG_Total'
            columnDetail1 = 'RDG-FA'
            columnName2 = 'MSC_TMatrix_Total'
            columnDetail2 = 'T-matrix'
            titleName = "T-matrix and RDG-FA Mass-scattering Coefficient (MSC) Ratio"
            titleAppend = ""
            self.PlotSettingTotalRatio(colName1=columnName1, colDetail1=columnDetail1,
                                       colName2=columnName2, colDetail2=columnDetail2,
                                       titleName=titleName, titleA=titleAppend,
                                       shareY='all', folderName=folderName,
                                       yScale='linear', yLabel=self.__yLabelRatioMSC,
                                       showValue=True, valueFormat='{:.2f}', yAxisFormat='%1.2f')
            ######################## SSA
            columnName1 = 'SSA_RDG_Total'
            columnDetail1 = 'RDG-FA'
            columnName2 = 'SSA_TMatrix_Total'
            columnDetail2 = 'T-matrix'
            titleName = "T-matrix and RDG-FA Single-scattering Albedo Ratio"
            titleAppend = ""
            self.PlotSettingTotalRatio(colName1=columnName1, colDetail1=columnDetail1,
                                       colName2=columnName2, colDetail2=columnDetail2,
                                       titleName=titleName, titleA=titleAppend,
                                       shareY='all', folderName=folderName,
                                       yScale='linear', yLabel=self.__yLabelRatioSSA,
                                       showValue=True, valueFormat='{:.2f}', yAxisFormat='%1.2f')


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

    def PlotTotalGraphsLine(self, mode1, model1, mode2, model2):
        try:
            columnNameA, columnDetailA, columnNameB, columnDetailB = self.ModeModelSelector(mode1, model1, mode2, model2)
            #########################################
            #########################################
            #########################################
            folderName = f'{mode1}_{model1}--{mode2}_{model2}/_TotalLineGraphs'
            ######################## ABS Cross Section
            columnName1 = 'ABS' + columnNameA
            columnDetail1 = columnDetailA
            columnName2 = 'ABS' + columnNameB
            columnDetail2 = columnDetailB
            titleName = "Absorption Cross-section"
            titleAppend = " Vs. Mobility-equivalent Diameter"
            self.PlotSettingTotalLines(colName1=columnName1, colDetail1=columnDetail1, mode1=mode1,
                                       colName2=columnName2, colDetail2=columnDetail2, mode2=mode2,
                                       titleName=titleName, titleA=titleAppend,
                                       shareY='all', folderName=folderName,
                                       yScale='log', yLabel=self.__yLabelPlotCrossSections,
                                       yAxisFormat='%1.1e')
            ######################## SCA Cross Section
            columnName1 = 'SCA' + columnNameA
            columnDetail1 = columnDetailA
            columnName2 = 'SCA' + columnNameB
            columnDetail2 = columnDetailB
            titleName = "Scattering Cross-section"
            titleAppend = " Vs. Mobility-equivalent Diameter"
            self.PlotSettingTotalLines(colName1=columnName1, colDetail1=columnDetail1, mode1=mode1,
                                       colName2=columnName2, colDetail2=columnDetail2, mode2=mode2,
                                       titleName=titleName, titleA=titleAppend,
                                       shareY='all', folderName=folderName,
                                       yScale='log', yLabel=self.__yLabelPlotCrossSections,
                                       yAxisFormat='%1.1e')
            ######################## MAC Cross Section
            columnName1 = 'MAC' + columnNameA
            columnDetail1 = columnDetailA
            columnName2 = 'MAC' + columnNameB
            columnDetail2 = columnDetailB
            titleName = "Mass-absorption Coefficient (MAC)"
            titleAppend = " Vs. Mobility-equivalent Diameter"
            self.PlotSettingTotalLines(colName1=columnName1, colDetail1=columnDetail1, mode1=mode1,
                                       colName2=columnName2, colDetail2=columnDetail2, mode2=mode2,
                                       titleName=titleName, titleA=titleAppend,
                                       shareY='all', folderName=folderName,
                                       yScale='linear', yLabel=self.__yLabelMAC,
                                       yAxisFormat='%1.2f')
            ######################## MSC Cross Section
            columnName1 = 'MSC' + columnNameA
            columnDetail1 = columnDetailA
            columnName2 = 'MSC' + columnNameB
            columnDetail2 = columnDetailB
            titleName = "Mass-scattering Coefficient (MSC)"
            titleAppend = " Vs. Mobility-equivalent Diameter"
            self.PlotSettingTotalLines(colName1=columnName1, colDetail1=columnDetail1, mode1=mode1,
                                       colName2=columnName2, colDetail2=columnDetail2, mode2=mode2,
                                       titleName=titleName, titleA=titleAppend,
                                       shareY='all', folderName=folderName,
                                       yScale='log', yLabel=self.__yLabelMSC,
                                       yAxisFormat='%1.2f')
            ######################## SSA
            columnName1 = 'SSA' + columnNameA
            columnDetail1 = columnDetailA
            columnName2 = 'SSA' + columnNameB
            columnDetail2 = columnDetailB
            titleName = "Single-scattering Albedo (SSA)"
            titleAppend = " Vs. Mobility-equivalent Diameter"
            self.PlotSettingTotalLines(colName1=columnName1, colDetail1=columnDetail1, mode1=mode1,
                                       colName2=columnName2, colDetail2=columnDetail2, mode2=mode2,
                                       titleName=titleName, titleA=titleAppend,
                                       shareY='all', folderName=folderName,
                                       yScale='log', yLabel=self.__yLabelSSA,
                                       yAxisFormat='%1.2f')
            ######################## Absorption Eff
            columnName1 = 'ABS' + columnNameA + '_Eff'
            columnDetail1 = columnDetailA
            columnName2 = 'ABS' + columnNameB + '_Eff'
            columnDetail2 = columnDetailB
            titleName = "Absorption Efficiency"
            titleAppend = " Vs. Mobility-equivalent Diameter"
            self.PlotSettingTotalLines(colName1=columnName1, colDetail1=columnDetail1, mode1=mode1,
                                       colName2=columnName2, colDetail2=columnDetail2, mode2=mode2,
                                       titleName=titleName, titleA=titleAppend,
                                       shareY='all', folderName=folderName,
                                       yScale='log', yLabel=self.__yLabelEFF,
                                       yAxisFormat='%1.2f')
            ######################## Scattering Eff
            columnName1 = 'SCA' + columnNameA + '_Eff'
            columnDetail1 = columnDetailA
            columnName2 = 'SCA' + columnNameB + '_Eff'
            columnDetail2 = columnDetailB
            titleName = "Scattering Efficiency"
            titleAppend = " Vs. Mobility-equivalent Diameter"
            self.PlotSettingTotalLines(colName1=columnName1, colDetail1=columnDetail1, mode1=mode1,
                                       colName2=columnName2, colDetail2=columnDetail2, mode2=mode2,
                                       titleName=titleName, titleA=titleAppend,
                                       shareY='all', folderName=folderName,
                                       yScale='log', yLabel=self.__yLabelEFF,
                                       yAxisFormat='%1.2f')
            ######################## Ratio ABS and SCA
            columnName1 = 'RatioABS'
            columnDetail1 = 'ABS Ratio'
            columnName2 = 'RatioSCA'
            columnDetail2 = 'SCA Ratio'
            titleName = "T-matrix and RDG-FA Cross Section Ratio"
            titleAppend = " Vs. Mobility-equivalent Diameter"
            self.PlotSettingTotalLines(colName1=columnName1, colDetail1=columnDetail1, mode1=mode1,
                                       colName2=columnName2, colDetail2=columnDetail2, mode2=mode2,
                                       titleName=titleName, titleA=titleAppend,
                                       shareY='all', folderName=folderName,
                                       yScale='linear', yLabel=self.__yLabelRatioCRS,
                                       yAxisFormat='%1.2f')
            ######################## dp median and dp Avg
            columnName1 = 'dp_median'
            columnDetail1 = 'Primary Particle Median Diameter'
            columnName2 = 'dp_Ave'
            columnDetail2 = 'Primary Particle Average Diameter'
            titleName = "Primary Particle Median and Average Diameter"
            titleAppend = " Vs. Mobility-equivalent Diameter"
            self.PlotSettingTotalLines(colName1=columnName1, colDetail1=columnDetail1, mode1=mode1,
                                       colName2=columnName2, colDetail2=columnDetail2, mode2=mode2,
                                       titleName=titleName, titleA=titleAppend,
                                       shareY='all', folderName=folderName,
                                       yScale='linear', yLabel=self.__yLabel_dp,
                                       yAxisFormat='%1.2f')

            A = 33
        except Exception as e:
            logging.exception(e)
            raise

    ###################################
    ###################################
    ###################################
    ###################################

    def PlotTotalGraphsBar(self, mode1, model1, mode2, model2):
        try:
            columnNameA, columnDetailA, columnNameB, columnDetailB = self.ModeModelSelector(mode1, model1, mode2, model2)
            #########################################
            #########################################
            #########################################
            folderName = f'{mode1}_{model1}--{mode2}_{model2}/__TotalBarGraphs'
            ######################## ABS Cross Section
            columnName1 = 'ABS' + columnNameA + '_Total'
            columnDetail1 = columnDetailA
            columnName2 = 'ABS' + columnNameB + '_Total'
            columnDetail2 = columnDetailB
            titleName = "Total Absorption Cross-section"
            titleAppend = ""
            self.PlotSettingTotalBars(colName1=columnName1, colDetail1=columnDetail1, mode1=mode1,
                                      colName2=columnName2, colDetail2=columnDetail2, mode2=mode2,
                                      titleName=titleName, titleA=titleAppend,
                                      shareY='all', folderName=folderName,
                                      yScale='linear', yLabel=self.__yLabelPlotCrossSections,
                                      showValue=True, valueFormat='{:.2E}', yAxisFormat='%1.1e')
            ######################## SCA Cross Section
            columnName1 = 'SCA' + columnNameA + '_Total'
            columnDetail1 = columnDetailA
            columnName2 = 'SCA' + columnNameB + '_Total'
            columnDetail2 = columnDetailB
            titleName = "Total Scattering Cross-section"
            titleAppend = ""
            self.PlotSettingTotalBars(colName1=columnName1, colDetail1=columnDetail1, mode1=mode1,
                                      colName2=columnName2, colDetail2=columnDetail2, mode2=mode2,
                                      titleName=titleName, titleA=titleAppend,
                                      shareY='all', folderName=folderName,
                                      yScale='linear', yLabel=self.__yLabelPlotCrossSections,
                                      showValue=True, valueFormat='{:.2E}', yAxisFormat='%1.1e')
            ######################## MAC
            columnName1 = 'MAC' + columnNameA + '_Total'
            columnDetail1 = columnDetailA
            columnName2 = 'MAC' + columnNameB + '_Total'
            columnDetail2 = columnDetailB
            titleName = "Total Mass-absorption Coefficient (MAC)"
            titleAppend = ""
            self.PlotSettingTotalBars(colName1=columnName1, colDetail1=columnDetail1, mode1=mode1,
                                      colName2=columnName2, colDetail2=columnDetail2, mode2=mode2,
                                      titleName=titleName, titleA=titleAppend,
                                      shareY='all', folderName=folderName,
                                      yScale='linear', yLabel=self.__yLabelMAC,
                                      showValue=True, valueFormat='{:.2f}', yAxisFormat='%1.2f')
            ######################## MSC
            columnName1 = 'MSC' + columnNameA + '_Total'
            columnDetail1 = columnDetailA
            columnName2 = 'MSC' + columnNameB + '_Total'
            columnDetail2 = columnDetailB
            titleName = "Total Mass-scattering Coefficient (MSC)"
            titleAppend = ""
            self.PlotSettingTotalBars(colName1=columnName1, colDetail1=columnDetail1, mode1=mode1,
                                      colName2=columnName2, colDetail2=columnDetail2, mode2=mode2,
                                      titleName=titleName, titleA=titleAppend,
                                      shareY='all', folderName=folderName,
                                      yScale='linear', yLabel=self.__yLabelMSC,
                                      showValue=True, valueFormat='{:.2f}', yAxisFormat='%1.2f')
            ######################## SSA
            columnName1 = 'SSA' + columnNameA + '_Total'
            columnDetail1 = columnDetailA
            columnName2 = 'SSA' + columnNameB + '_Total'
            columnDetail2 = columnDetailB
            titleName = "Total Single-scattering Albedo"
            titleAppend = ""
            self.PlotSettingTotalBars(colName1=columnName1, colDetail1=columnDetail1, mode1=mode1,
                                      colName2=columnName2, colDetail2=columnDetail2, mode2=mode2,
                                      titleName=titleName, titleA=titleAppend,
                                      shareY='all', folderName=folderName,
                                      yScale='linear', yLabel=self.__yLabelSSA,
                                      showValue=True, valueFormat='{:.2f}', yAxisFormat='%1.2f')


        except Exception as e:
            logging.exception(e)
            raise

    ###################################
    ###################################
    ###################################
    ###################################

    def PlotdpInfo(self):
        try:
            folderName = '_PrimaryParticle'
            ######################## dp median
            columnName = 'dp_median'
            titleName = "Primary Particle Median Diameter for each Mobility Diameter"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='linear', yLabel=self.__yLabel_dpMedian)
            ######################## dp ave
            columnName = 'dp_Ave'
            titleName = "Primary Particle Average Diameter for each Mobility Diameter"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='linear', yLabel=self.__yLabel_dpAve)
            ######################## Number of calcs
            columnName = 'NumberOfCalcs'
            titleName = "Number of Calculations for each Mobility Diameter"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='linear', yLabel=self.__yLabelCalcNumber)

        except Exception as e:
            logging.exception(e)
            raise

    ###################################
    ###################################
    ###################################
    ###################################

    def PlotEff(self):
        try:
            folderName = '_OpticalEfficiency'
            ######################## ABS RDG
            columnName = 'ABS_RDG_Eff'
            titleName = "Absorption Efficiency for RDG-FA"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='linear', yLabel=self.__yLabelEFF)
            ######################## SCA RDG
            columnName = 'SCA_RDG_Eff'
            titleName = "Scattering Efficiency for RDG-FA"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='linear', yLabel=self.__yLabelEFF)
            ######################## ABS TMatrix
            columnName = 'ABS_TMatrix_Eff'
            titleName = "Absorption Efficiency for T-matrix"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='linear', yLabel=self.__yLabelEFF)
            ######################## SCA TMatrix
            columnName = 'SCA_TMatrix_Eff'
            titleName = "Scattering Efficiency for T-matrix"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='linear', yLabel=self.__yLabelEFF)



        except Exception as e:
            logging.exception(e)
            raise

    ###################################
    ###################################
    ###################################
    ###################################

    def PlotMACMSC(self):
        try:
            folderName = '_MACs_MSCs'
            ######################## MAC RDG
            columnName = 'MAC_RDG'
            titleName = "MAC for RDG-FA"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='linear', yLabel=self.__yLabelMAC)
            ######################## MSC RDG
            columnName = 'MSC_RDG'
            titleName = "MSC for RDG-FA"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='linear', yLabel=self.__yLabelMSC)
            ######################## MAC TMatrix
            columnName = 'MAC_TMatrix'
            titleName = "MAC for T-matrix"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='linear', yLabel=self.__yLabelMAC)
            ######################## MSC TMatrix
            columnName = 'MSC_TMatrix'
            titleName = "MSC for T-matrix"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='linear', yLabel=self.__yLabelMSC)



        except Exception as e:
            logging.exception(e)
            raise

    ###################################
    ###################################
    ###################################
    ###################################

    def PlotSSA(self):
        try:
            folderName = '_SSAs'
            ######################## SSA TMatrix
            columnName = 'SSA_TMatrix'
            titleName = "SSA for T-matrix"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='linear', yLabel=self.__yLabelSSA)

            ######################## SSA RDG
            columnName = 'SSA_RDG'
            titleName = "SSA for RDG-FA"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='linear', yLabel=self.__yLabelSSA)



        except Exception as e:
            logging.exception(e)
            raise

    ###################################
    ###################################
    ###################################
    ###################################

    def PlotErrorAndRatio(self):
        try:
            folderName = '_ErrorsAndRatios'
            ######################## Absorption Percentage Error
            columnName = 'RealPercentErrorABS'
            titleName = "Percentage Error Between T-matrix and RDG-FA Absorption Cross-section"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='linear', yLabel=self.__yLabelErrorPercent)
            ######################## Absorption Ratio
            columnName = 'RatioABS'
            titleName = "T-matrix and RDG-FA Absorption Cross-section Ratio"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='linear', yLabel=self.__yLabelRatioCRS)
            ######################## Scattering Percentage Error
            columnName = 'RealPercentErrorSCA'
            titleName = "Percentage Error Between T-matrix and RDG-FA Scattering Cross-section"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='linear', yLabel=self.__yLabelErrorPercent)
            ######################## Scattering Ratio
            columnName = 'RatioSCA'
            titleName = "T-matrix and RDG-FA Scattering Cross-section Ratio"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='linear', yLabel=self.__yLabelRatioCRS)


        except Exception as e:
            logging.exception(e)
            raise

    ###################################
    ###################################
    ###################################
    ###################################

    def PlotCrossSections(self):
        try:
            folderName = '_OpticalCrossSections'
            ######################## Absorption RDG
            columnName = 'ABS_RDG'
            titleName = "RDG-FA Absorption Cross-section"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='log', yLabel=self.__yLabelPlotCrossSections)
            ######################## Scattering RDG
            columnName = 'SCA_RDG'
            titleName = "RDG-FA Scattering Cross-section"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='log', yLabel=self.__yLabelPlotCrossSections)
            ######################## Absorption TMatrix
            columnName = 'ABS_TMatrix'
            titleName = "T-matrix Absorption Cross-section"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='log', yLabel=self.__yLabelPlotCrossSections)
            ######################## Scattering TMatrix
            columnName = 'SCA_TMatrix'
            titleName = "T-matrix Scattering Cross-section"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='log', yLabel=self.__yLabelPlotCrossSections)

        except Exception as e:
            logging.exception(e)
            raise

    ##################################################################################
    ##################################################################################
    ##################################################################################
    ##################################################################################

    def PlotSimpleLines(self, A1, A1Name, A2, A2Name,
                        title, column, folderName,
                        titleAppend1, roundNTitle,
                        fileAppend1, roundNFile,
                        yLabel,
                        titleAppend2=None, fileAppend2=None,
                        xScale='log', yScale='log'):
        try:
            alphaMainLine = self.__alphaMainLine

            for i in A1.unique():
                fig, ax1 = plt.subplots()
                c_lineColor = 0
                c_lineWidth = 0
                for j in A2.unique():
                    fileNames = self.dfMainInfo[(A1Name == i) & (A2Name == j)].AA_FileName
                    c_lineStyle = 0
                    c_markerStyle = 0

                    for index, fileN in fileNames.iteritems():
                        if fileN in self.dictData:
                            df = self.dictData[fileN]
                            ax1.plot(df['dm'], df[column], label=self.dictDataLabel[fileN],
                                     color=self.__lineColor[c_lineColor], linewidth=self.__lineWidth[c_lineWidth],
                                     linestyle=self.__lineStyle[c_lineStyle], alpha=alphaMainLine,
                                     marker=self.__markerStyle[c_markerStyle], markersize=self.__markerSize)
                        c_lineStyle += 1
                        c_markerStyle += 1

                    c_lineColor += 1
                    c_lineWidth += 1

                if titleAppend2:
                    T1 = title + titleAppend1 + str(round(i, roundNTitle)) + titleAppend2
                else:
                    T1 = title + titleAppend1 + str(round(i, roundNTitle))
                if fileAppend2:
                    F1 = title + fileAppend1 + str(round(i, roundNFile)) + fileAppend2
                else:
                    F1 = title + fileAppend1 + str(round(i, roundNFile))

                ax1.grid(True, which='both', axis="both", alpha=0.5)
                ax1.set_xscale(xScale)
                ax1.set_yscale(yScale)
                ax1.set_xlabel(self.__xLabelMobDiameter, fontsize=self.__xLabelGeneralFontSize)
                ax1.set_ylabel(yLabel, fontsize=self.__yLabelGeneralFontSize)
                ax1.set_xlim(self.__xAxisLimitsComp[0], self.__xAxisLimitsComp[1])

                ax1.set_title(T1, fontsize=self.__plotTitleGeneralFontSize)
                ax1.xaxis.set_major_formatter(FormatStrFormatter("%i"))
                ax1.xaxis.set_minor_formatter(FormatStrFormatter("%i"))
                ax1.tick_params(axis='x', which='major', labelsize=self.__xMajorTickLabelGeneralFontSize)
                ax1.tick_params(axis='x', which='minor', labelsize=self.__xMinorTickLabelGeneralFontSize)
                ax1.tick_params(axis='y', which='major', labelsize=self.__yMajorTickLabelGeneralFontSize)
                ax1.tick_params(axis='y', which='minor', labelsize=self.__yMinorTickLabelGeneralFontSize)
                ax1.legend(bbox_to_anchor=(1.00, 0.5), markerscale=self.__legendMarkerScale, loc='center left', fontsize='large')
                self.SaveAndClosePlot(folderName=folderName, F1=F1)

        except Exception as e:
            logging.exception(e)
            raise

    def PlotSetting_3CatLine(self, columnName, titleName, folderName,
                             yScale, yLabel):
        try:
            if self.__rhoEff100nmCTE == True:
                self.PlotSimpleLines(A1=self.ser_rhoEff100nm, A1Name=self.dfMainInfo.AGG_EFF_RHO_100NM_CENTER,
                                     A2=self.ser_Dm, A2Name=self.dfMainInfo.AGG_EFF_DM_CENTER,
                                     column=columnName, title=titleName, folderName=folderName,
                                     yScale=yScale, yLabel=yLabel,
                                     titleAppend1=" for " + "$\\rho_{eff,100}$= ", roundNTitle=1, titleAppend2=" kg/m$^3$",
                                     fileAppend1=" for " + "rho_eff=", roundNFile=0)

            if self.__DmCTE == True:
                self.PlotSimpleLines(A1=self.ser_Dm, A1Name=self.dfMainInfo.AGG_EFF_DM_CENTER,
                                     A2=self.ser_SigmaMob, A2Name=self.dfMainInfo.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER,
                                     column=columnName, title=titleName, folderName=folderName,
                                     yScale=yScale, yLabel=yLabel,
                                     titleAppend1=" for " + "$D_m=$", roundNTitle=2,
                                     fileAppend1=" for " + "Dm=", roundNFile=2)

            if self.__sigmaEachMobCTE == True:
                self.PlotSimpleLines(A1=self.ser_SigmaMob, A1Name=self.dfMainInfo.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER,
                                     A2=self.ser_rhoEff100nm, A2Name=self.dfMainInfo.AGG_EFF_RHO_100NM_CENTER,
                                     column=columnName, title=titleName, folderName=folderName,
                                     yScale=yScale, yLabel=yLabel,
                                     titleAppend1=" for " + "$\sigma_p|d_m=$", roundNTitle=2,
                                     fileAppend1=" for " + "sigma_p=", roundNFile=2)

        except Exception as e:
            logging.exception(e)
            raise

    ###################################
    ###################################
    ###################################
    ###################################

    def PlotMultipleBarsTotal(self, A1, A1Name, A2, A2Name,
                              column1, column1T, mode1,
                              column2, column2T, mode2,
                              title, titleAppend,
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
                # DmCount=0
                for j in A2.unique():
                    fileNames = self.dfMainInfo[(A1Name == i) & (A2Name == j)].AA_FileName

                    indexBarPlot = np.arange(3)

                    C1Height = []
                    C2Height = []

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

                    AXES_A = ax1[rowCount, colCount].bar(indexBarPlot - self.__barWidth / 2, C1Height, color=self.__A1Color, width=self.__barWidth, edgecolor='white', hatch=self.__patterns[0],
                                                         label=column1T)
                    AXES_B = ax1[rowCount, colCount].bar(indexBarPlot + self.__barWidth / 2, C2Height, color=self.__A2Color, width=self.__barWidth, edgecolor='white', hatch=self.__patterns[4],
                                                         label=column2T)
                    # DmCount+=1

                    if showValue:
                        autoLabel(self, rects=AXES_A, row=rowCount, col=colCount, format=valueFormat)
                        autoLabel(self, rects=AXES_B, row=rowCount, col=colCount, format=valueFormat)

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

            F1 = f"{title} for median={median}-sigma={sigma}"
            T1 = title + titleAppend + " for " + "$d_{m,g}=$" + str(round(median)) + " nm " + "and " + "$\sigma_g=$" + str(round(sigma, 2))

            fig.subplots_adjust(top=self.__subplotGridSetup[0], wspace=self.__subplotGridSetup[1], hspace=self.__subplotGridSetup[2])
            fig.suptitle(T1, fontsize=self.__plotTitleFontSizeTotal)

            fig.text(0.5, self.__xLabelCommonSubplotAdjustment, self.__xLabelMassMob, ha='center', fontsize=self.__xLabelCommonSubplotFontSize)

            if 'e' in yAxisFormat:
                fig.text(self.__yLabelCommonSubplotAdjustment[0], 0.5, yCommonLabel, va='center', rotation='vertical', fontsize=self.__yLabelCommonSubplotFontSize)
            elif 'f' in yAxisFormat:
                fig.text(self.__yLabelCommonSubplotAdjustment[1], 0.5, yCommonLabel, va='center', rotation='vertical', fontsize=self.__yLabelCommonSubplotFontSize)

            leg = plt.legend(bbox_to_anchor=(self.__legendSubplotAdjustment[0], self.__legendSubplotAdjustment[1]), markerscale=self.__legendSubplotMarkerScale, fontsize=self.__legendSubplotFontSize,
                             loc='lower left', borderaxespad=0.)
            ################# set the lineWidth of each legend object
            for legObj in leg.legendHandles:
                legObj.set_linewidth(self.__legendSubplotLineWidth)
            #################
            self.SaveAndClosePlot(folderName=folderName, F1=F1)


        except Exception as e:
            logging.exception(e)
            raise

    def PlotSettingTotalBars(self, colName1, colDetail1, mode1,
                             colName2, colDetail2, mode2,
                             titleName, titleA, shareY, folderName,
                             yScale, yLabel, showValue, valueFormat, yAxisFormat):
        try:
            for sigma in self.__dmSigmaLogN:
                for median in self.__dmMedianLogN:
                    self.PlotMultipleBarsTotal(A1=self.ser_rhoEff100nm, A1Name=self.dfMainInfo.AGG_EFF_RHO_100NM_CENTER,
                                               A2=self.ser_SigmaMob, A2Name=self.dfMainInfo.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER,
                                               column1=colName1, column1T=colDetail1, mode1=mode1,
                                               column2=colName2, column2T=colDetail2, mode2=mode2,
                                               title=titleName, titleAppend=titleA,
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

    def PlotMultipleLinesTotal(self, A1, A1Name, A2, A2Name,
                               column1, column1T, mode1,
                               column2, column2T, mode2,
                               title, titleAppend,
                               folderName, yCommonLabel,
                               shareY, yAxisFormat, yScale='log'):
        try:

            fig, ax1 = plt.subplots(nrows=len(A1.unique()), ncols=len(A2.unique()), sharex=True, sharey=shareY, figsize=(12, 12), constrained_layout=True)

            alphaMainLine = self.__alphaMainLine
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
                        df1 = self.dictData[mode1 + fileN]
                        ax1[rowCount, colCount].plot(df1['dm'], df1[column1], label=column1T + ', ' + Dmlist[DmCount],
                                                     color=self.get_color(self.__A1Color, 'yellow', portion=(DmCount + 1) * 0.1), linewidth=self.__A1LineWidth[c_lineWidth],
                                                     linestyle=self.__subplotLineStyle[c_lineStyle], alpha=self.__A1AlphaMainLine)

                        df2 = self.dictData[mode2 + fileN]
                        ax1[rowCount, colCount].plot(df2['dm'], df2[column2], label=column2T + ', ' + Dmlist[DmCount],
                                                     color=self.get_color(self.__A2Color, 'gray', portion=(DmCount + 1) * 0.2), linewidth=self.__A2LineWidth[c_lineWidth],
                                                     linestyle=self.__subplotLineStyle[c_lineStyle], alpha=self.__A2AlphaMainLine)

                        c_lineStyle += 1
                        c_markerStyle += 1
                        c_lineWidth += 1
                        DmCount += 1

                    ax1[rowCount, colCount].grid(True, which='both', axis="both", alpha=0.5)
                    ax1[rowCount, colCount].set_xscale('log')
                    ax1[rowCount, colCount].set_yscale(yScale)
                    ax1[rowCount, colCount].xaxis.set_major_formatter(FormatStrFormatter("%i"))
                    ax1[rowCount, colCount].xaxis.set_major_locator(plt.MaxNLocator(6))
                    ax1[rowCount, colCount].tick_params(axis='x', which='both', labelsize=self.__xMajorTickLabelSubplotFontSize, labelrotation=self.__xTickLabelSubplotRotation)
                    ax1[rowCount, colCount].tick_params(axis='y', which='both', labelsize=self.__yMajorTickLabelSubplotFontSize)
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

            F1 = f"{title}"
            T1 = title + titleAppend

            fig.subplots_adjust(top=self.__subplotGridSetup[0], wspace=self.__subplotGridSetup[1], hspace=self.__subplotGridSetup[2])
            fig.suptitle(T1, fontsize=self.__plotTitleFontSizeTotal)

            fig.text(0.5, self.__xLabelCommonSubplotAdjustment, self.__xLabelMobDiameter, ha='center', fontsize=self.__xLabelCommonSubplotFontSize)

            if 'e' in yAxisFormat:
                fig.text(self.__yLabelCommonSubplotAdjustment[0], 0.5, yCommonLabel, va='center', rotation='vertical', fontsize=self.__yLabelCommonSubplotFontSize)
            elif 'f' in yAxisFormat:
                fig.text(self.__yLabelCommonSubplotAdjustment[1], 0.5, yCommonLabel, va='center', rotation='vertical', fontsize=self.__yLabelCommonSubplotFontSize)

            leg = plt.legend(bbox_to_anchor=(self.__legendSubplotAdjustment[0], self.__legendSubplotAdjustment[1]), markerscale=self.__legendSubplotMarkerScale, fontsize=self.__legendSubplotFontSize,
                             loc='lower left', borderaxespad=0.)
            ################# set the lineWidth of each legend object
            for legObj in leg.legendHandles:
                legObj.set_linewidth(self.__legendSubplotLineWidth)
            #################
            self.SaveAndClosePlot(folderName=folderName, F1=F1)

        except Exception as e:
            logging.exception(e)
            raise

    def PlotSettingTotalLines(self, colName1, colDetail1, mode1,
                              colName2, colDetail2, mode2,
                              titleName, titleA, shareY, folderName,
                              yScale, yLabel, yAxisFormat):
        try:
            self.PlotMultipleLinesTotal(A1=self.ser_rhoEff100nm, A1Name=self.dfMainInfo.AGG_EFF_RHO_100NM_CENTER,
                                        A2=self.ser_SigmaMob, A2Name=self.dfMainInfo.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER,
                                        column1=colName1, column1T=colDetail1, mode1=mode1,
                                        column2=colName2, column2T=colDetail2, mode2=mode2,
                                        title=titleName, titleAppend=titleA,
                                        folderName=folderName,
                                        yScale=yScale, yCommonLabel=yLabel, yAxisFormat=yAxisFormat,
                                        shareY=shareY,
                                        )

        except Exception as e:
            logging.exception(e)
            raise

    ###################################
    ###################################
    ###################################
    ###################################

    def PlotMultipleBarsRatio(self, A1, A1Name, A2, A2Name,
                              column1, column1T,
                              column2, column2T,
                              title, titleAppend,
                              folderName, yCommonLabel,
                              shareY, yAxisFormat,
                              sigma, median,
                              showValue=True, valueFormat='{:.2E}', yScale='log'):
        try:

            fig, ax1 = plt.subplots(nrows=len(A1.unique()), ncols=len(A2.unique()), sharex=True, sharey=shareY, figsize=(12, 12))

            def autoLabel(self, rects, row, col, format='{:.2E}', xPos='center'):

                if 'E' in format:
                    fontS = self.__valueFontSizeTotal
                if 'f' in format:
                    fontS = self.__valueFontSizeTotal * 1.4
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

                    if sigma:
                        Ratio = {}
                        for median1 in self.__dmMedianLogN:
                            Ratio[median1] = []
                    if median:
                        Ratio = {}
                        for sigma1 in self.__dmSigmaLogN:
                            Ratio[sigma1] = []

                    for index, fileN in fileNames.iteritems():
                        if fileN in self.dictData:

                            df = self.dictTotalOpticalProp[fileN]
                            if sigma:
                                for median1 in self.__dmMedianLogN:

                                    selectedRow = df.loc[(df['sigma'] == sigma) & (df['median'] == median1)]
                                    if (len(selectedRow[column1]) != 1) or (len(selectedRow[column2]) != 1):
                                        raise
                                    else:
                                        for c1 in selectedRow[column1]:
                                            for c2 in selectedRow[column2]:
                                                r1 = c2 / c1
                                    Ratio[median1].append(r1)

                            if median:
                                for sigma1 in self.__dmSigmaLogN:

                                    selectedRow = df.loc[(df['sigma'] == sigma1) & (df['median'] == median)]
                                    if (len(selectedRow[column1]) != 1) or (len(selectedRow[column2]) != 1):
                                        raise
                                    else:
                                        for c1 in selectedRow[column1]:
                                            for c2 in selectedRow[column2]:
                                                r1 = c2 / c1
                                    Ratio[sigma1].append(r1)

                    A, B = [], []
                    for key in Ratio:
                        if sigma:
                            A.append("$d_{m,g}=$" + str(round(key)) + " nm ")
                        if median:
                            A.append("$\sigma_g=$" + str(round(key, 2)))
                        B.append(Ratio[key])

                    AXES_A = ax1[rowCount, colCount].bar(indexBarPlot - self.__barRatioWidth, B[0], color=self.__barColor[0], width=self.__barRatioWidth, edgecolor='white', label=A[0],
                                                         hatch=self.__patterns[0])
                    AXES_B = ax1[rowCount, colCount].bar(indexBarPlot, B[1], color=self.__barColor[1], width=self.__barRatioWidth, edgecolor='white', label=A[1], hatch=self.__patterns[4])
                    AXES_C = ax1[rowCount, colCount].bar(indexBarPlot + self.__barRatioWidth, B[2], color=self.__barColor[2], width=self.__barRatioWidth, edgecolor='white', label=A[2],
                                                         hatch=self.__patterns[5])
                    if showValue:
                        autoLabel(self, rects=AXES_A, row=rowCount, col=colCount, format=valueFormat)
                        autoLabel(self, rects=AXES_B, row=rowCount, col=colCount, format=valueFormat)
                        autoLabel(self, rects=AXES_C, row=rowCount, col=colCount, format=valueFormat)

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

            if median:
                F1 = f"{title} for median={median}"
                T1 = title + titleAppend + " for " + "$d_{m,g}=$" + str(round(median)) + " nm "
            if sigma:
                F1 = f"{title} for sigma={sigma}"
                T1 = title + titleAppend + " for " + "$\sigma_g=$" + str(round(sigma, 2))

            fig.subplots_adjust(top=self.__subplotGridSetup[0], wspace=self.__subplotGridSetup[1], hspace=self.__subplotGridSetup[2])
            fig.suptitle(T1, fontsize=self.__plotTitleFontSizeTotal)

            fig.text(0.5, self.__xLabelCommonSubplotAdjustment, self.__xLabelMassMob, ha='center', fontsize=self.__xLabelCommonSubplotFontSize)

            if 'e' in yAxisFormat:
                fig.text(self.__yLabelCommonSubplotAdjustment[0], 0.5, yCommonLabel, va='center', rotation='vertical', fontsize=self.__yLabelCommonSubplotFontSize)
            elif 'f' in yAxisFormat:
                fig.text(self.__yLabelCommonSubplotAdjustment[1], 0.5, yCommonLabel, va='center', rotation='vertical', fontsize=self.__yLabelCommonSubplotFontSize)

            leg = plt.legend(bbox_to_anchor=(self.__legendSubplotAdjustment[0], self.__legendSubplotAdjustment[1]), markerscale=self.__legendSubplotMarkerScale, fontsize=self.__legendSubplotFontSize,
                             loc='lower left', borderaxespad=0.)
            ################# set the lineWidth of each legend object
            for legObj in leg.legendHandles:
                legObj.set_linewidth(self.__legendSubplotLineWidth)
            #################
            self.SaveAndClosePlot(folderName=folderName, F1=F1)

        except Exception as e:
            logging.exception(e)
            raise

    def PlotSettingTotalRatio(self, colName1, colDetail1, colName2, colDetail2,
                              titleName, titleA, shareY, folderName,
                              yScale, yLabel, showValue, valueFormat, yAxisFormat):
        try:
            for sigma in self.__dmSigmaLogN:
                self.PlotMultipleBarsRatio(A1=self.ser_rhoEff100nm, A1Name=self.dfMainInfo.AGG_EFF_RHO_100NM_CENTER,
                                           A2=self.ser_SigmaMob, A2Name=self.dfMainInfo.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER,
                                           column1=colName1, column1T=colDetail1,
                                           column2=colName2, column2T=colDetail2,
                                           title=titleName, titleAppend=titleA,
                                           folderName=folderName,
                                           yScale=yScale, yCommonLabel=yLabel, yAxisFormat=yAxisFormat,
                                           shareY=shareY,
                                           sigma=sigma, median=None,
                                           showValue=showValue, valueFormat=valueFormat)
            for median in self.__dmMedianLogN:
                self.PlotMultipleBarsRatio(A1=self.ser_rhoEff100nm, A1Name=self.dfMainInfo.AGG_EFF_RHO_100NM_CENTER,
                                           A2=self.ser_SigmaMob, A2Name=self.dfMainInfo.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER,
                                           column1=colName1, column1T=colDetail1,
                                           column2=colName2, column2T=colDetail2,
                                           title=titleName, titleAppend=titleA,
                                           folderName=folderName,
                                           yScale=yScale, yCommonLabel=yLabel, yAxisFormat=yAxisFormat,
                                           shareY=shareY,
                                           sigma=None, median=median,
                                           showValue=showValue, valueFormat=valueFormat)


        except Exception as e:
            logging.exception(e)
            raise

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
