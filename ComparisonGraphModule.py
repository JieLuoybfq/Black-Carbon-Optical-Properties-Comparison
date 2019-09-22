from ConfigParserModule import logging
import GeneralFunctions as GF
import numpy as np
from matplotlib import rcParams
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.ticker import FormatStrFormatter
from math import pi
from math import exp
from math import isnan
from math import log
from scipy.optimize import fsolve
import subprocess

####### Plotting Parameters
rcParams['mathtext.fontset'] = 'stix'
rcParams['font.family'] = 'STIXGeneral'


##############

class GraphTools:
    def __init__(self, FolderInfo):
        try:
            #################################################
            #################################################
            ################################################# dm and dp distribution variables
            self.__dpDeviation = True
            self.__dpDeviationPercent = 25
            self.__dmSigmaLogN = [1.2, 1.4, 1.6]
            self.__dmMedianLogN = [100, 200, 300]
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
            #################################################
            ################################################# Experiment Graphs
            self.__CRSGraphs_EXP = True
            self.__ErrorAndRatioGraphs_EXP = True
            self.__MACMSCGraphs_EXP = True
            self.__SSAGraphs_EXP = True
            self.__dpMedianGraphs_EXP = True
            self.__DmCTE_EXP = True
            self.__rhoEff100nmCTE_EXP = True
            self.__sigmaMobCTE_EXP = True
            #################################################
            self.__TotalDFLabels = ['sigma', 'median', 'TotalMassGram',
                                    'ABS_RDG_Total', 'SCA_RDG_Total',
                                    'ABS_TMatrix_Total', 'SCA_TMatrix_Total',
                                    'MAC_RDG_Total', 'MSC_RDG_Total',
                                    'MAC_TMatrix_Total', 'MSC_TMatrix_Total',
                                    'SSA_RDG_Total', 'SSA_TMatrix_Total', 'sum']
            ################################################# General Variables
            self.__folderNameGraph = FolderInfo['FOLDER_NAME_GRAPH']
            self.__folderNameData = FolderInfo['FOLDER_NAME_DATA']
            self.__WMF_SVGSaving = False
            self.__inkScapePath = "C://Program Files//inkscape//inkscape.exe"
            self.__xLabelMobDiameter = "Mobility-equivalent Diameter (nm)"
            self.__yLabelPlotCrossSections = "Cross Section (µm" + "$^{}$".format(2) + ")"
            self.__xLabelMassMob = "Mass-Mobility Exponent"
            self.__yLabelErrorPercent = "Error Percent (%)"
            self.__yLabelRatioCRS = "Cross Section Ratio (TMatrix/RDG)"
            self.__yLabelRatioMAC = "MAC Ratio (TMatrix/RDG)"
            self.__yLabelRatioMSC = "MSC Ratio(TMatrix/RDG)"
            self.__yLabelRatioSSA = "SSA Ratio (TMatrix/RDG)"
            self.__yLabelSSA = "SSA"
            self.__yLabelMAC = "MAC (m$^2$/g)"
            self.__yLabelMSC = "MSC (m$^2$/g)"
            self.__yLabelEFF = "Optical Efficiency"
            self.__yLabel_dpMedian = "Primary Particle Median Diameter (nm)"
            self.__yLabel_dpAve = "Primary Particle Average Diameter (nm)"
            self.__yLabel_dp = "Primary Particle Diameter (nm)"
            self.__yLabelCalcNumber = "#"
            ################################################# Plot Setting
            self.__plotTitleGeneralFontSize = 14
            self.__xLabelGeneralFontSize = 12
            self.__yLabelGeneralFontSize = 12
            self.__xMajorTickLabelGeneralFontSize = 9
            self.__xMinorTickLabelGeneralFontSize = 8
            self.__yMajorTickLabelGeneralFontSize = 9
            self.__yMinorTickLabelGeneralFontSize = 8
            self.__legendMarkerScale = 2
            self.__figureDPI = 500
            self.__markerSize = 3
            self.__alphaMainLine = 0.55
            self.__lineColor = ['red', 'blue', 'green']
            self.__barColor = ['red', 'blue', 'green']
            self.__lineStyle = ['-', '--', ':']
            self.__markerStyle = ["o", "X", "^"]
            self.__lineWidth = [1.5, 1, 0.5]
            self.__xAxisLimitsComp = [49, 960]
            self.__xAxisLimits2 = [49, 960]
            #####
            self.__barWidth = 0.34
            self.__barRatioWidth = 0.25
            self.__valueFontSizeTotal = 7
            self.__plotTitleFontSizeTotal = 30
            self.__xLabelGroupFontSizeTotal = 18
            self.__xLabelMultipleSubplotFontSize = 27
            #####
            self.__yTickFontSizeTotal = 17
            self.__yLabelMultipleSubplotFontSize = 25
            #####
            self.__xLabelCommonFontSizeTotal = 12
            self.__yLabelCommonFontSizeTotal = 12

        except Exception as e:
            logging.exception(e)
            raise

    def ExeecuteMobility(self, AGG_Info):
        try:
            def getAverageAndSTD(Array):
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

            def getRandomFromArr(Array, Number):
                try:
                    # uniform choice
                    A = np.random.choice(Array, size=int(Number), replace=False)
                    return A

                except Exception as e:
                    logging.exception(e)
                    raise

            def createRandomNormalArr(Center, Width, Number):
                try:
                    A = np.random.normal(Center, Width, int(Number))
                    return A

                except Exception as e:
                    logging.exception(e)
                    raise

            def calcViscosity(Temperature, Pressure):  # ASSUME INDEPENDENT OF P HERE, Sutherland Equation
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

            def calcMonteCarlo_Viscosity(Length, Temperature, Pressure):
                try:
                    results = np.empty(shape=int(Length))
                    for i in range(int(Length)):
                        results[i] = calcViscosity(Temperature=Temperature[i], Pressure=Pressure[i])
                    return results

                except Exception as e:
                    logging.exception(e)
                    raise

            def calcMeanFreePath(Viscosity, Temperature, Pressure):
                try:
                    # For Air
                    M = .029  # kg/mol/K
                    R = 8.314  # J/mol/K
                    l = 2 * Viscosity / (Pressure * (8 * M / pi / R / Temperature) ** 0.5)
                    return l

                except Exception as e:
                    logging.exception(e)
                    raise

            def calcMonteCarlo_MeanFreePath(Length, Viscosity, Temperature, Pressure):
                try:
                    results = np.empty(shape=int(Length))
                    for i in range(int(Length)):
                        results[i] = calcMeanFreePath(Viscosity=Viscosity[i], Temperature=Temperature[i], Pressure=Pressure[i])
                    return results

                except Exception as e:
                    logging.exception(e)
                    raise

            def calcMobilityDiameter(AerodynamicDiameter_nm, MeanFreePath_m, Eff_Dm, Eff_Rho_100nm):
                try:
                    da = AerodynamicDiameter_nm * (10 ** (-9))
                    func = lambda mobilityDiameterIn_nm: ((mobilityDiameterIn_nm ** 2) * (
                            1 + (2 * MeanFreePath_m / mobilityDiameterIn_nm) * (1.257 + .4 * exp(-1.1 / (2 * MeanFreePath_m / mobilityDiameterIn_nm)))) * (
                                                                  (Eff_Rho_100nm / ((100 * (10 ** -9)) ** (Eff_Dm - 3))) * (mobilityDiameterIn_nm) ** (Eff_Dm - 3))) - (
                                                                 (da ** 2) * (1 + (2 * MeanFreePath_m / da) * (1.257 + .4 * exp(-1.1 / (2 * MeanFreePath_m / da)))) * 1000)
                    dm_initial_guess = da
                    dm_solution = fsolve(func, dm_initial_guess)
                    dm_Nano = dm_solution[0] * (10 ** 9)
                    return dm_Nano

                except Exception as e:
                    logging.exception(e)
                    raise

            def calcMonteCarlo_MobilityDiameter(Length, AerodynamicDiameter_nm, MeanFreePath_m, Eff_Dm, Eff_Rho_100nm):
                try:
                    results = np.empty(shape=int(Length))
                    for i in range(int(Length)):
                        results[i] = calcMobilityDiameter(AerodynamicDiameter_nm=AerodynamicDiameter_nm[i], MeanFreePath_m=MeanFreePath_m[i], Eff_Dm=Eff_Dm[i], Eff_Rho_100nm=Eff_Rho_100nm[i])
                    return results

                except Exception as e:
                    logging.exception(e)
                    raise

            def calcAbsorptionCS(AbsorptionCoefficient, NumberDensity):
                try:
                    # in um^2 if AbsorptionCoefficient in(1/Mm) and NumberDensity in (#/cc)
                    A = AbsorptionCoefficient / NumberDensity
                    return A

                except Exception as e:
                    logging.exception(e)
                    raise

            def calcMonteCarlo_AbsorptionCS(Length, AbsorptionCoefficient, NumberDensity):
                try:
                    results = np.empty(shape=int(Length))
                    for i in range(int(Length)):
                        results[i] = calcAbsorptionCS(AbsorptionCoefficient=AbsorptionCoefficient[i], NumberDensity=NumberDensity[i])
                    return results

                except Exception as e:
                    logging.exception(e)
                    raise

            def calcScatteringCS(ScatteringCoefficient, NumberDensity):
                try:
                    # in um^2 if ScatteringCoefficient in(1/Mm) and NumberDensity in (#/cc)
                    A = ScatteringCoefficient / NumberDensity
                    return A

                except Exception as e:
                    logging.exception(e)
                    raise

            def calcMonteCarlo_ScatteringCS(Length, ScatteringCoefficient, NumberDensity):
                try:
                    results = np.empty(shape=int(Length))
                    for i in range(int(Length)):
                        results[i] = calcScatteringCS(ScatteringCoefficient=ScatteringCoefficient[i], NumberDensity=NumberDensity[i])
                    return results

                except Exception as e:
                    logging.exception(e)
                    raise

            def calcSSA(Scattering, Absorption):
                try:
                    # in um^2 if ScatteringCoefficient in(1/Mm) and NumberDensity in (#/cc)
                    A = Scattering / (Scattering + Absorption)
                    return A

                except Exception as e:
                    logging.exception(e)
                    raise

            def calcMonteCarlo_SSA(Length, Scattering, Absorption):
                try:
                    results = np.empty(shape=int(Length))
                    for i in range(int(Length)):
                        results[i] = calcSSA(Scattering=Scattering[i], Absorption=Absorption[i])
                    return results

                except Exception as e:
                    logging.exception(e)
                    raise

            ###########################################
            ###########################################
            ###########################################
            ###########################################
            ###########################################
            ###########################################
            self.ExpDF = pd.read_csv(GF.getAddressTo(FolderName=self.__folderNameData, FileName='Experiment', Extension='csv'))
            AeroDiam = self.ExpDF['AAC_Diam_Ave']
            AAC_Diam_STD = self.ExpDF['AAC_Diam_STD']
            PAX_ABS_Coef_Ave = self.ExpDF['PAX_ABS_Ave']
            PAX_ABS_Coef_STD = self.ExpDF['PAX_ABS_STD']
            PAX_SCA_Coef_Ave = self.ExpDF['PAX_SCA_Ave']
            PAX_SCA_Coef_STD = self.ExpDF['PAX_SCA_STD']
            CPC_Ave = self.ExpDF['CPC_Conc_Ave']
            CPC_STD = self.ExpDF['CPC_Conc_STD']
            ##################
            ##################
            dm = []
            dm_STD = []
            ABS_CRS = []
            ABS_CRS_STD = []
            SCA_CRS = []
            SCA_CRS_STD = []
            SSA = []
            SSA_STD = []
            ##################
            ##################
            for index, diam_Aero in AeroDiam.iteritems():
                arrAgg_Aerodynamic_Diameter = createRandomNormalArr(Center=AeroDiam[index],
                                                                    Width=AAC_Diam_STD[index],
                                                                    Number=AGG_Info['MONTECARLO_ARRAY_SIZE'])
                arrAgg_Eff_Dm = createRandomNormalArr(Center=AGG_Info['AGG_EFF_DM_CENTER'],
                                                      Width=AGG_Info['AGG_EFF_DM_STANDARD_DEVIATION'],
                                                      Number=AGG_Info['MONTECARLO_ARRAY_SIZE'])
                arrAgg_Eff_Rho_100nm = createRandomNormalArr(Center=AGG_Info['AGG_EFF_RHO_100NM_CENTER'],
                                                             Width=AGG_Info['AGG_EFF_RHO_100NM_STANDARD_DEVIATION'],
                                                             Number=AGG_Info['MONTECARLO_ARRAY_SIZE'])
                arrAir_Temperature = createRandomNormalArr(Center=AGG_Info['AIR_TEMPERATURE_CENTER'],
                                                           Width=AGG_Info['AIR_TEMPERATURE_STANDARD_DEVIATION'],
                                                           Number=AGG_Info['MONTECARLO_ARRAY_SIZE'])
                arrAir_Pressure = createRandomNormalArr(Center=AGG_Info['AIR_PRESSURE_CENTER'],
                                                        Width=AGG_Info['AIR_PRESSURE_STANDARD_DEVIATION'],
                                                        Number=AGG_Info['MONTECARLO_ARRAY_SIZE'])
                arrPAX_Absorption_Coefficient = createRandomNormalArr(Center=PAX_ABS_Coef_Ave[index],
                                                                      Width=PAX_ABS_Coef_STD[index],
                                                                      Number=AGG_Info['MONTECARLO_ARRAY_SIZE'])
                arrPAX_Scattering_Coefficient = createRandomNormalArr(Center=PAX_SCA_Coef_Ave[index],
                                                                      Width=PAX_SCA_Coef_STD[index],
                                                                      Number=AGG_Info['MONTECARLO_ARRAY_SIZE'])
                arrCPC_Number_Density = createRandomNormalArr(Center=CPC_Ave[index],
                                                              Width=CPC_STD[index],
                                                              Number=AGG_Info['MONTECARLO_ARRAY_SIZE'])

                ######################################################

                arrAgg_Aerodynamic_Diameter_Random = getRandomFromArr(Array=arrAgg_Aerodynamic_Diameter, Number=AGG_Info['MONTECARLO_RANDOM_SIZE'])
                arrAgg_Eff_Dm_Random = getRandomFromArr(Array=arrAgg_Eff_Dm, Number=AGG_Info['MONTECARLO_RANDOM_SIZE'])
                arrAgg_Eff_Rho_100nm_Random = getRandomFromArr(Array=arrAgg_Eff_Rho_100nm, Number=AGG_Info['MONTECARLO_RANDOM_SIZE'])
                arrAir_Temperature_Random = getRandomFromArr(Array=arrAir_Temperature, Number=AGG_Info['MONTECARLO_RANDOM_SIZE'])
                arrAir_Pressure_Random = getRandomFromArr(Array=arrAir_Pressure, Number=AGG_Info['MONTECARLO_RANDOM_SIZE'])
                arrPAX_Absorption_Coefficient_Random = getRandomFromArr(Array=arrPAX_Absorption_Coefficient, Number=AGG_Info['MONTECARLO_RANDOM_SIZE'])
                arrPAX_Scattering_Coefficient_Random = getRandomFromArr(Array=arrPAX_Scattering_Coefficient, Number=AGG_Info['MONTECARLO_RANDOM_SIZE'])
                arrCPC_Number_Density_Random = getRandomFromArr(Array=arrCPC_Number_Density, Number=AGG_Info['MONTECARLO_RANDOM_SIZE'])

                ####################### Calculation
                arrMonteCarloViscosity = calcMonteCarlo_Viscosity(Length=AGG_Info['MONTECARLO_RANDOM_SIZE'], Temperature=arrAir_Temperature_Random,
                                                                  Pressure=arrAir_Pressure_Random)
                arrMonteCarloMeanFreePath = calcMonteCarlo_MeanFreePath(Length=AGG_Info['MONTECARLO_RANDOM_SIZE'],
                                                                        Viscosity=arrMonteCarloViscosity,
                                                                        Temperature=arrAir_Temperature_Random,
                                                                        Pressure=arrAir_Pressure_Random)
                # arrViscosityAverageSTD = getAverageAndSTD(arrMonteCarloViscosity)
                # arrMeanFreePathAverageSTD = getAverageAndSTD(arrMonteCarloMeanFreePath)

                # Mobility diameter in nm
                arrMonteCarloMobilityDiameter = calcMonteCarlo_MobilityDiameter(Length=AGG_Info['MONTECARLO_RANDOM_SIZE'],
                                                                                AerodynamicDiameter_nm=arrAgg_Aerodynamic_Diameter_Random,
                                                                                MeanFreePath_m=arrMonteCarloMeanFreePath,
                                                                                Eff_Dm=arrAgg_Eff_Dm_Random,
                                                                                Eff_Rho_100nm=arrAgg_Eff_Rho_100nm_Random)
                arrMobilityDiameterAverageSTD = getAverageAndSTD(arrMonteCarloMobilityDiameter)

                # Cross section in um^2
                arrMonteCarloAbsorptionCS = calcMonteCarlo_AbsorptionCS(Length=AGG_Info['MONTECARLO_RANDOM_SIZE'],
                                                                        AbsorptionCoefficient=arrPAX_Absorption_Coefficient_Random,
                                                                        NumberDensity=arrCPC_Number_Density_Random)
                arrMonteCarloScatteringCS = calcMonteCarlo_ScatteringCS(Length=AGG_Info['MONTECARLO_RANDOM_SIZE'],
                                                                        ScatteringCoefficient=arrPAX_Scattering_Coefficient_Random,
                                                                        NumberDensity=arrCPC_Number_Density_Random)
                arrAbsorptionCSAverageSTD = getAverageAndSTD(arrMonteCarloAbsorptionCS)
                arrScatteringCSAverageSTD = getAverageAndSTD(arrMonteCarloScatteringCS)
                arrMonteCarloSSA = calcMonteCarlo_SSA(Length=AGG_Info['MONTECARLO_RANDOM_SIZE'], Scattering=arrMonteCarloScatteringCS, Absorption=arrMonteCarloAbsorptionCS)
                arrSSAAverageSTD = getAverageAndSTD(arrMonteCarloSSA)
                #######################

                dm.append(arrMobilityDiameterAverageSTD[0])
                dm_STD.append(arrMobilityDiameterAverageSTD[1])
                ABS_CRS.append(arrAbsorptionCSAverageSTD[0])
                ABS_CRS_STD.append(arrAbsorptionCSAverageSTD[1])
                SCA_CRS.append(arrScatteringCSAverageSTD[0])
                SCA_CRS_STD.append(arrScatteringCSAverageSTD[1])
                SSA.append(arrSSAAverageSTD[0])
                SSA_STD.append(arrSSAAverageSTD[1])
            ExpDict = {'dm': dm,
                       'dm_STD': dm_STD,
                       'ABS_CRS': ABS_CRS,
                       'ABS_CRS_STD': ABS_CRS_STD,
                       'SCA_CRS': SCA_CRS,
                       'SCA_CRS_STD': SCA_CRS_STD,
                       'SSA': SSA,
                       'SSA_STD': SSA_STD
                       }
            self.ExpDataFrame = pd.DataFrame(ExpDict)
        except Exception as e:
            logging.exception(e)
            raise

    def Experiment_RDG_TMatrixComparisonGraphs(self, AGG_Info):
        try:
            self.ExeecuteMobility(AGG_Info)
            ###########################################
            ###########################################
            ###########################################

            self.dfMainInfo = pd.read_csv(GF.getAddressTo(FolderName=self.__folderNameData, FileName='beacon', Extension='csv'))
            self.dictData = {}
            self.dictDataLabel = {}
            self.serMainFileName = self.dfMainInfo.loc[self.dfMainInfo['AA_Plot'] == 1, 'AA_FileName']
            self.ser_Dm = self.dfMainInfo.loc[self.dfMainInfo['AA_Plot'] == 1, 'AGG_EFF_DM_CENTER']
            self.ser_rhoEff100nm = self.dfMainInfo.loc[self.dfMainInfo['AA_Plot'] == 1, 'AGG_EFF_RHO_100NM_CENTER']
            self.ser_SigmaMob = self.dfMainInfo.loc[self.dfMainInfo['AA_Plot'] == 1, 'AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER']
            self.ser_D_TEM = self.dfMainInfo.loc[self.dfMainInfo['AA_Plot'] == 1, 'D_TEM']
            self.ser_dp100_nano = self.dfMainInfo.loc[self.dfMainInfo['AA_Plot'] == 1, 'dp100_nano']

            for index, item in self.serMainFileName.iteritems():
                self.dictData[item] = pd.read_csv(GF.getAddressTo(FolderName=self.__folderNameData, FileName=item))
                self.dictData[item] = self.dictData[item].replace(0, np.nan)
                self.dictDataLabel[item] = "$D_m=$" + str(round(self.ser_Dm.loc[index], 2)) + ", " + "$\\rho_{eff,100}$=" + str(round(self.ser_rhoEff100nm.loc[index], 1)) \
                                           + ", " + "$\sigma_p|d_m=$" + str(round(self.ser_SigmaMob.loc[index], 2)) + ", " + "$D_{TEM}=$" + str(round(self.ser_D_TEM.loc[index], 2)) \
                                           + ", " + "$d_{p,100}=$" + str(round(self.ser_dp100_nano.loc[index], 1))

            if self.__CRSGraphs_EXP:
                self.PlotCrossSections_EXP()
            if self.__ErrorAndRatioGraphs_EXP:
                pass
                # self.PlotErrorAndRatio_EXP()
            if self.__SSAGraphs_EXP:
                self.PlotSSA_EXP()
            if self.__MACMSCGraphs_EXP:
                pass
                # self.PlotMACMSC_EXP()
            if self.__dpMedianGraphs_EXP:
                pass
                # self.PlotdpMedian_EXP()

        except Exception as e:
            logging.exception(e)
            raise
        ###################################
        ###################################
        ###################################
        ###################################

    def PlotSSACore_EXP(self, A1, A1Name, A2, A2Name, T1A, T1I, F1A, F1I, title, column, yLabel, mean, STD, T1B=None, F1B=None):
        try:
            alphaMainLine = self.__alphaMainLine

            for i in A1.unique():

                fig, ax1 = plt.subplots()

                c_lineWidth = 0

                for j in A2.unique():

                    fileNames = self.dfMainInfo[(A1Name == i) & (A2Name == j)].AA_FileName

                    c_lineStyle = 0
                    c_markerStyle = 0
                    c_lineColor = 0

                    for index, item in fileNames.iteritems():
                        if item in self.dictData:
                            df = self.dictData[item]
                            ax1.plot(df['dm'], df[column], label=self.dictDataLabel[item], color=self.__lineColor[c_lineColor], linewidth=self.__lineWidth[c_lineWidth],
                                     linestyle=self.__lineStyle[c_lineStyle], alpha=alphaMainLine, marker=self.__markerStyle[c_markerStyle], markersize=self.__markerSize)

                        c_lineStyle += 1
                        c_markerStyle += 1
                        c_lineColor += 1
                    # ax1.plot(self.ExpDataFrame['dm'], self.ExpDataFrame[mean],)
                    dm_STD = self.ExpDataFrame['dm_STD'] * 2
                    y_STD = self.ExpDataFrame[STD] * 2
                    ax1.errorbar(self.ExpDataFrame['dm'], self.ExpDataFrame[mean], xerr=dm_STD, yerr=y_STD, label="Experiment", color='black',
                                 marker="D", markersize=self.__markerSize, ecolor='black', capsize=4, elinewidth=1.5, markeredgewidth=1.5)

                    c_lineWidth += 1
                if T1B:
                    T1 = title + T1A + str(round(i, T1I)) + T1B
                else:
                    T1 = title + T1A + str(round(i, T1I))
                if F1B:
                    F1 = title + F1A + str(round(i, F1I)) + F1B
                else:
                    F1 = title + F1A + str(round(i, F1I))

                ax1.grid(True, which='both', axis="both", alpha=0.5)
                ax1.set_xscale("log")
                # ax1.set_yscale("log")
                ax1.set_xlabel(self.__xLabelMobDiameter, fontsize=self.__xLabelGeneralFontSize)
                ax1.set_ylabel(yLabel, fontsize=self.__yLabelGeneralFontSize)
                ax1.set_xlim(self.__xAxisLimitsComp[0], self.__xAxisLimitsComp[1])

                ax1.set_title(T1, fontsize=self.__plotTitleGeneralFontSize)
                ax1.xaxis.set_major_formatter(FormatStrFormatter("%i"))
                ax1.xaxis.set_minor_formatter(FormatStrFormatter("%i"))
                ax1.legend(bbox_to_anchor=(1.00, 0.5), loc='center left', fontsize='large')
                Address = GF.getAddressTo(FolderName=self.__folderNameGraph + "\SSAs_EXP", FileName=F1, Extension="jpg")
                plt.savefig(Address, format='jpg', dpi=self.__figureDPI, bbox_inches='tight')
                plt.clf()
                plt.close()

        except Exception as e:
            logging.exception(e)
            raise

    def PlotSSAIntermediate_EXP(self, clName, tlName, yLabel, mean, STD):
        try:
            if self.__rhoEff100nmCTE == True:
                self.PlotSSACore_EXP(A1=self.ser_rhoEff100nm, A1Name=self.dfMainInfo.AGG_EFF_RHO_100NM_CENTER, A2=self.ser_Dm, A2Name=self.dfMainInfo.AGG_EFF_DM_CENTER,
                                     column=clName, mean=mean, STD=STD, yLabel=yLabel, title=tlName, T1A=" for " + "$\\rho_{eff,100}$= ", T1I=1, T1B=" kg/m$^3$", F1A=" for " + "rho_eff=", F1I=0)

            # if self.__DmCTE == True:
            #     self.PlotSSACore(A1=self.DmSeries, A1Name=self.MainDF.AGG_EFF_DM_CENTER, A2=self.sigmaMobSeries, A2Name=self.MainDF.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER,
            #                      column=clName, yLabel=yLabel, title=tlName, T1A=" for " + "$D_m=$", T1I=2, F1A=" for " + "Dm=", F1I=2)
            #
            # if self.__sigmaMobCTE == True:
            #     self.PlotSSACore(A1=self.sigmaMobSeries, A1Name=self.MainDF.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER, A2=self.rhoEff100nmSeries,
            #                      A2Name=self.MainDF.AGG_EFF_RHO_100NM_CENTER,
            #                      column=clName, yLabel=yLabel, title=tlName, T1A=" for " + "$\sigma_p|d_m=$", T1I=2, F1A=" for " + "sigma_p=", F1I=2)

        except Exception as e:
            logging.exception(e)
            raise

    def PlotSSA_EXP(self):
        try:
            ######################## SSA TMatrix
            clName = 'SSA_TMatrix'
            mean = 'SSA'
            STD = 'SSA_STD'
            tlName = "SSA for TMatrix"
            self.PlotSSAIntermediate_EXP(clName=clName, tlName=tlName, yLabel=self.__yLabelSSA, mean=mean, STD=STD)

            ######################## SSA RDG
            clName = 'SSA_RDG'
            mean = 'SSA'
            STD = 'SSA_STD'
            tlName = "SSA for RDG-FA"
            self.PlotSSAIntermediate_EXP(clName=clName, tlName=tlName, yLabel=self.__yLabelSSA, mean=mean, STD=STD)



        except Exception as e:
            logging.exception(e)
            raise

    ###################################
    ###################################
    ###################################
    ###################################
    def PlotCrossSectionsCore_EXP(self, A1, A1Name, A2, A2Name, T1A, T1I, F1A, F1I, title, column, mean, STD, T1B=None, F1B=None):
        try:
            alphaMainLine = self.__alphaMainLine

            for i in A1.unique():

                fig, ax1 = plt.subplots()

                c_lineWidth = 0

                for j in A2.unique():

                    fileNames = self.dfMainInfo[(A1Name == i) & (A2Name == j)].AA_FileName

                    c_lineStyle = 0
                    c_markerStyle = 0
                    c_lineColor = 0

                    for index, item in fileNames.iteritems():
                        if item in self.dictData:
                            df = self.dictData[item]
                            ax1.plot(df['dm'], df[column], label=self.dictDataLabel[item], color=self.__lineColor[c_lineColor], linewidth=self.__lineWidth[c_lineWidth],
                                     linestyle=self.__lineStyle[c_lineStyle], alpha=alphaMainLine, marker=self.__markerStyle[c_markerStyle], markersize=self.__markerSize)

                        c_lineStyle += 1
                        c_markerStyle += 1
                        c_lineColor += 1
                    # ax1.plot(self.ExpDataFrame['dm'], self.ExpDataFrame[mean],)
                    dm_STD = self.ExpDataFrame['dm_STD'] * 2
                    y_STD = self.ExpDataFrame[STD] * 2
                    ax1.errorbar(self.ExpDataFrame['dm'], self.ExpDataFrame[mean], xerr=dm_STD, yerr=y_STD, label="Experiment", color='black',
                                 marker="D", markersize=self.__markerSize, ecolor='black', capsize=4, elinewidth=1.5, markeredgewidth=1.5)

                    c_lineWidth += 1
                if T1B:
                    T1 = title + T1A + str(round(i, T1I)) + T1B
                else:
                    T1 = title + T1A + str(round(i, T1I))
                if F1B:
                    F1 = title + F1A + str(round(i, F1I)) + F1B
                else:
                    F1 = title + F1A + str(round(i, F1I))

                ax1.grid(True, which='both', axis="both", alpha=0.5)
                ax1.set_xscale("log")
                ax1.set_yscale("log")
                ax1.set_xlabel(self.__xLabelMobDiameter, fontsize=self.__xLabelGeneralFontSize)
                ax1.set_ylabel(self.__yLabelPlotCrossSections, fontsize=self.__yLabelGeneralFontSize)
                ax1.set_xlim(self.__xAxisLimitsComp[0], self.__xAxisLimitsComp[1])

                ax1.set_title(T1, fontsize=self.__plotTitleGeneralFontSize)
                ax1.xaxis.set_major_formatter(FormatStrFormatter("%i"))
                ax1.xaxis.set_minor_formatter(FormatStrFormatter("%i"))
                ax1.legend(bbox_to_anchor=(1.00, 0.5), loc='center left', fontsize='large')
                Address = GF.getAddressTo(FolderName=self.__folderNameGraph + "\CrossSections_EXP", FileName=F1, Extension="jpg")
                plt.savefig(Address, format='jpg', dpi=self.__figureDPI, bbox_inches='tight')
                plt.clf()
                plt.close()

        except Exception as e:
            logging.exception(e)
            raise

    def PlotCrossSectionsIntermediate_EXP(self, clName, tlName, mean, STD):
        try:
            if self.__rhoEff100nmCTE_EXP == True:
                self.PlotCrossSectionsCore_EXP(A1=self.ser_rhoEff100nm, A1Name=self.dfMainInfo.AGG_EFF_RHO_100NM_CENTER, A2=self.ser_Dm, A2Name=self.dfMainInfo.AGG_EFF_DM_CENTER,
                                               column=clName, mean=mean, STD=STD, title=tlName, T1A=" for " + "$\\rho_{eff,100}$= ", T1I=1, T1B=" kg/m$^3$", F1A=" for " + "rho_eff=", F1I=0)

            # if self.__DmCTE == True:
            #     self.PlotSimpleLines(A1=self.DmSeries, A1Name=self.MainDF.AGG_EFF_DM_CENTER, A2=self.sigmaMobSeries, A2Name=self.MainDF.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER,
            #                                column=clName, title=tlName, T1A=" for " + "$D_m=$", T1I=2, F1A=" for " + "Dm=", F1I=2)
            #
            # if self.__sigmaMobCTE == True:
            #     self.PlotSimpleLines(A1=self.sigmaMobSeries, A1Name=self.MainDF.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER, A2=self.rhoEff100nmSeries,
            #                                A2Name=self.MainDF.AGG_EFF_RHO_100NM_CENTER,
            #                                column=clName, title=tlName, T1A=" for " + "$\sigma_p|d_m=$", T1I=2, F1A=" for " + "sigma_p=", F1I=2)

        except Exception as e:
            logging.exception(e)
            raise

    def PlotCrossSections_EXP(self):
        try:

            ######################## Absorption RDG
            clName = 'ABS_RDG'
            mean = 'ABS_CRS'
            STD = 'ABS_CRS_STD'
            tlName = "Experiment and RDG-FA Absorption Cross Section"
            self.PlotCrossSectionsIntermediate_EXP(clName=clName, tlName=tlName, mean=mean, STD=STD)
            ######################## Scattering RDG
            clName = 'SCA_RDG'
            mean = 'SCA_CRS'
            STD = 'SCA_CRS_STD'
            tlName = "Experiment and RDG-FA Scattering Cross Section"
            self.PlotCrossSectionsIntermediate_EXP(clName=clName, tlName=tlName, mean=mean, STD=STD)
            ######################## Absorption TMatrix
            clName = 'ABS_TMatrix'
            mean = 'ABS_CRS'
            STD = 'ABS_CRS_STD'
            tlName = "Experiment and TMatrix Absorption Cross Section"
            self.PlotCrossSectionsIntermediate_EXP(clName=clName, tlName=tlName, mean=mean, STD=STD)
            ######################## Scattering TMatrix
            clName = 'SCA_TMatrix'
            mean = 'SCA_CRS'
            STD = 'SCA_CRS_STD'
            tlName = "Experiment and TMatrix Scattering Cross Section"
            self.PlotCrossSectionsIntermediate_EXP(clName=clName, tlName=tlName, mean=mean, STD=STD)


        except Exception as e:
            logging.exception(e)
            raise

    ##########################################################################################
    ##########################################################################################
    ##########################################################################################
    ##########################################################################################
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
            x = 5
            c = 2.5
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

                    if sum >= 0.90:
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
            self.dfMainInfo = pd.read_csv(GF.getAddressTo(FolderName=self.__folderNameData, FileName='beacon', Extension='csv'))
            self.dictData = {}
            self.dictDataLabel = {}
            self.dictTotalOpticalProp = {}

            self.GetDataDefinition()

            for index, fileN in self.serMainFileName.iteritems():

                self.dictData[fileN] = pd.read_csv(GF.getAddressTo(FolderName=self.__folderNameData, FileName=fileN))
                self.dictData[fileN] = self.dictData[fileN].replace(0, np.nan)

                if self.__dpDeviation:
                    self.Check_dpRatio(fileName=fileN)

                self.CalcTotalMACMSC(self.dictData[fileN], fileN)

                self.CalcEfficiency(fileName=fileN)

                self.GetDataLabel(fileName=fileN, index=index)

            if self.__CRSGraphs:
                self.PlotCrossSections()

            if self.__ErrorAndRatioGraphs:
                self.PlotErrorAndRatio()

            if self.__SSAGraphs:
                self.PlotSSA()

            if self.__MACMSCGraphs:
                self.PlotMACMSC()

            if self.__EfficiencyGraphs:
                self.PlotEff()

            if self.__dpInfoGraphs:
                self.PlotdpInfo()

            if self.__dmDistributionGraphsBar:
                self.PlotTotalGraphsBar()

            if self.__dmDistributionGraphsLine:
                self.PlotTotalGraphsLine()

            if self.__dmDistributionRatioGraphsBar:
                self.PlotTotalRatioGraphsBar()

        except Exception as e:
            logging.exception(e)
            raise

    ###################################
    ###################################
    ###################################
    ###################################

    def PlotTotalRatioGraphsBar(self):
        try:
            folderName = 'RatioBarGraphs'
            ######################## ABS Cross Section
            columnName1 = 'ABS_RDG_Total'
            columnDetail1 = 'RDG-FA'
            columnName2 = 'ABS_TMatrix_Total'
            columnDetail2 = 'T-matrix'
            titleName = "T-matrix and RDG-FA Absorption Cross Section Ratio"
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
            titleName = "T-matrix and RDG-FA Scattering Cross Section Ratio"
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

    def PlotTotalGraphsLine(self):
        try:
            folderName = 'TotalLineGraphs'
            ######################## ABS Cross Section
            columnName1 = 'ABS_RDG'
            columnDetail1 = 'RDG-FA'
            columnName2 = 'ABS_TMatrix'
            columnDetail2 = 'T-matrix'
            titleName = "Absorption Cross Section"
            titleAppend = " in µm" + "$^{}$".format(2) + " Vs. Mobility-equivalent Diameter"
            self.PlotSettingTotalLines(colName1=columnName1, colDetail1=columnDetail1,
                                       colName2=columnName2, colDetail2=columnDetail2,
                                       titleName=titleName, titleA=titleAppend,
                                       shareY='all', folderName=folderName,
                                       yScale='log', yLabel=self.__yLabelPlotCrossSections,
                                       yAxisFormat='%1.1e')
            ######################## SCA Cross Section
            columnName1 = 'SCA_RDG'
            columnDetail1 = 'RDG-FA'
            columnName2 = 'SCA_TMatrix'
            columnDetail2 = 'T-matrix'
            titleName = "Scattering Cross Section"
            titleAppend = " in µm" + "$^{}$".format(2) + " Vs. Mobility-equivalent Diameter"
            self.PlotSettingTotalLines(colName1=columnName1, colDetail1=columnDetail1,
                                       colName2=columnName2, colDetail2=columnDetail2,
                                       titleName=titleName, titleA=titleAppend,
                                       shareY='all', folderName=folderName,
                                       yScale='log', yLabel=self.__yLabelPlotCrossSections,
                                       yAxisFormat='%1.1e')
            ######################## MAC Cross Section
            columnName1 = 'MAC_RDG'
            columnDetail1 = 'RDG-FA'
            columnName2 = 'MAC_TMatrix'
            columnDetail2 = 'T-matrix'
            titleName = "Mass-absorption Coefficient (MAC)"
            titleAppend = " in m$^2$/g" + " Vs. Mobility-equivalent Diameter"
            self.PlotSettingTotalLines(colName1=columnName1, colDetail1=columnDetail1,
                                       colName2=columnName2, colDetail2=columnDetail2,
                                       titleName=titleName, titleA=titleAppend,
                                       shareY='all', folderName=folderName,
                                       yScale='linear', yLabel=self.__yLabelMAC,
                                       yAxisFormat='%1.2f')
            ######################## MSC Cross Section
            columnName1 = 'MSC_RDG'
            columnDetail1 = 'RDG-FA'
            columnName2 = 'MSC_TMatrix'
            columnDetail2 = 'T-matrix'
            titleName = "Mass-scattering Coefficient (MSC)"
            titleAppend = " in m$^2$/g" + " Vs. Mobility-equivalent Diameter"
            self.PlotSettingTotalLines(colName1=columnName1, colDetail1=columnDetail1,
                                       colName2=columnName2, colDetail2=columnDetail2,
                                       titleName=titleName, titleA=titleAppend,
                                       shareY='all', folderName=folderName,
                                       yScale='linear', yLabel=self.__yLabelMSC,
                                       yAxisFormat='%1.2f')
            ######################## SSA
            columnName1 = 'SSA_RDG'
            columnDetail1 = 'RDG-FA'
            columnName2 = 'SSA_TMatrix'
            columnDetail2 = 'T-matrix'
            titleName = "Single-scattering Albedo"
            titleAppend = "" + " Vs. Mobility-equivalent Diameter"
            self.PlotSettingTotalLines(colName1=columnName1, colDetail1=columnDetail1,
                                       colName2=columnName2, colDetail2=columnDetail2,
                                       titleName=titleName, titleA=titleAppend,
                                       shareY='all', folderName=folderName,
                                       yScale='linear', yLabel=self.__yLabelSSA,
                                       yAxisFormat='%1.2f')
            ######################## Ratio ABS and SCA
            columnName1 = 'RatioABS'
            columnDetail1 = 'ABS Ratio'
            columnName2 = 'RatioSCA'
            columnDetail2 = 'SCA Ratio'
            titleName = "T-matrix and RDG-FA Cross Section Ratio"
            titleAppend = "" + " Vs. Mobility-equivalent Diameter"
            self.PlotSettingTotalLines(colName1=columnName1, colDetail1=columnDetail1,
                                       colName2=columnName2, colDetail2=columnDetail2,
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
            titleAppend = "" + " Vs. Mobility-equivalent Diameter"
            self.PlotSettingTotalLines(colName1=columnName1, colDetail1=columnDetail1,
                                       colName2=columnName2, colDetail2=columnDetail2,
                                       titleName=titleName, titleA=titleAppend,
                                       shareY='all', folderName=folderName,
                                       yScale='linear', yLabel=self.__yLabel_dp,
                                       yAxisFormat='%1.2f')
            ######################## Absorption Eff
            columnName1 = 'ABS_RDG_Eff'
            columnDetail1 = 'RDG'
            columnName2 = 'ABS_TMatrix_Eff'
            columnDetail2 = 'T-matrix'
            titleName = "T-matrix and RDG-FA Absorption Efficiency"
            titleAppend = "" + " Vs. Mobility-equivalent Diameter"
            self.PlotSettingTotalLines(colName1=columnName1, colDetail1=columnDetail1,
                                       colName2=columnName2, colDetail2=columnDetail2,
                                       titleName=titleName, titleA=titleAppend,
                                       shareY='all', folderName=folderName,
                                       yScale='linear', yLabel=self.__yLabelEFF,
                                       yAxisFormat='%1.2f')
            ######################## Scattering Eff
            columnName1 = 'SCA_RDG_Eff'
            columnDetail1 = 'RDG'
            columnName2 = 'SCA_TMatrix_Eff'
            columnDetail2 = 'T-matrix'
            titleName = "T-matrix and RDG-FA Scattering Efficiency"
            titleAppend = "" + " Vs. Mobility-equivalent Diameter"
            self.PlotSettingTotalLines(colName1=columnName1, colDetail1=columnDetail1,
                                       colName2=columnName2, colDetail2=columnDetail2,
                                       titleName=titleName, titleA=titleAppend,
                                       shareY='all', folderName=folderName,
                                       yScale='linear', yLabel=self.__yLabelEFF,
                                       yAxisFormat='%1.2f')
            A = 33



        except Exception as e:
            logging.exception(e)
            raise

    ###################################
    ###################################
    ###################################
    ###################################

    def PlotTotalGraphsBar(self):
        try:
            folderName = 'TotalBarGraphs'
            ######################## ABS Cross Section
            columnName1 = 'ABS_RDG_Total'
            columnDetail1 = 'RDG-FA'
            columnName2 = 'ABS_TMatrix_Total'
            columnDetail2 = 'T-matrix'
            titleName = "Total Absorption Cross Section"
            titleAppend = " in µm" + "$^{}$".format(2)
            self.PlotSettingTotalBars(colName1=columnName1, colDetail1=columnDetail1,
                                      colName2=columnName2, colDetail2=columnDetail2,
                                      titleName=titleName, titleA=titleAppend,
                                      shareY='all', folderName=folderName,
                                      yScale='linear', yLabel=self.__yLabelPlotCrossSections,
                                      showValue=True, valueFormat='{:.2E}', yAxisFormat='%1.1e')
            ######################## SCA Cross Section
            columnName1 = 'SCA_RDG_Total'
            columnDetail1 = 'RDG-FA'
            columnName2 = 'SCA_TMatrix_Total'
            columnDetail2 = 'T-matrix'
            titleName = "Total Scattering Cross Section"
            titleAppend = " in µm" + "$^{}$".format(2)
            self.PlotSettingTotalBars(colName1=columnName1, colDetail1=columnDetail1,
                                      colName2=columnName2, colDetail2=columnDetail2,
                                      titleName=titleName, titleA=titleAppend,
                                      shareY='all', folderName=folderName,
                                      yScale='linear', yLabel=self.__yLabelPlotCrossSections,
                                      showValue=True, valueFormat='{:.2E}', yAxisFormat='%1.1e')
            ######################## MAC
            columnName1 = 'MAC_RDG_Total'
            columnDetail1 = 'RDG-FA'
            columnName2 = 'MAC_TMatrix_Total'
            columnDetail2 = 'T-matrix'
            titleName = "Mass-absorption Coefficient (MAC)"
            titleAppend = " in m$^2$/g"
            self.PlotSettingTotalBars(colName1=columnName1, colDetail1=columnDetail1,
                                      colName2=columnName2, colDetail2=columnDetail2,
                                      titleName=titleName, titleA=titleAppend,
                                      shareY='all', folderName=folderName,
                                      yScale='linear', yLabel=self.__yLabelMAC,
                                      showValue=True, valueFormat='{:.2f}', yAxisFormat='%1.2f')
            ######################## MSC
            columnName1 = 'MSC_RDG_Total'
            columnDetail1 = 'RDG-FA'
            columnName2 = 'MSC_TMatrix_Total'
            columnDetail2 = 'T-matrix'
            titleName = "Mass-scattering Coefficient (MSC)"
            titleAppend = " in m$^2$/g"
            self.PlotSettingTotalBars(colName1=columnName1, colDetail1=columnDetail1,
                                      colName2=columnName2, colDetail2=columnDetail2,
                                      titleName=titleName, titleA=titleAppend,
                                      shareY='all', folderName=folderName,
                                      yScale='linear', yLabel=self.__yLabelMSC,
                                      showValue=True, valueFormat='{:.2f}', yAxisFormat='%1.2f')
            ######################## SSA
            columnName1 = 'SSA_RDG_Total'
            columnDetail1 = 'RDG-FA'
            columnName2 = 'SSA_TMatrix_Total'
            columnDetail2 = 'T-matrix'
            titleName = "Single-scattering Albedo"
            titleAppend = ""
            self.PlotSettingTotalBars(colName1=columnName1, colDetail1=columnDetail1,
                                      colName2=columnName2, colDetail2=columnDetail2,
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
            folderName = 'PrimaryParticle'
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
            folderName = 'OpticalEfficiency'
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
            folderName = 'MACs_MSCs'
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
            folderName = 'SSAs'
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
            folderName = 'ErrorsAndRatios'
            ######################## Absorption Percentage Error
            columnName = 'RealPercentErrorABS'
            titleName = "Percentage Error Between T-matrix and RDG-FA Absorption Cross Section"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='linear', yLabel=self.__yLabelErrorPercent)
            ######################## Absorption Ratio
            columnName = 'RatioABS'
            titleName = "T-matrix and RDG-FA Absorption Cross Section Ratio"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='linear', yLabel=self.__yLabelRatioCRS)
            ######################## Scattering Percentage Error
            columnName = 'RealPercentErrorSCA'
            titleName = "Percentage Error Between T-matrix and RDG-FA Scattering Cross Section"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='linear', yLabel=self.__yLabelErrorPercent)
            ######################## Scattering Ratio
            columnName = 'RatioSCA'
            titleName = "T-matrix and RDG-FA Scattering Cross Section Ratio"
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
            folderName = 'OpticalCrossSections'
            ######################## Absorption RDG
            columnName = 'ABS_RDG'
            titleName = "RDG-FA Absorption Cross Section"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='log', yLabel=self.__yLabelPlotCrossSections)
            ######################## Scattering RDG
            columnName = 'SCA_RDG'
            titleName = "RDG-FA Scattering Cross Section"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='log', yLabel=self.__yLabelPlotCrossSections)
            ######################## Absorption TMatrix
            columnName = 'ABS_TMatrix'
            titleName = "T-matrix Absorption Cross Section"
            self.PlotSetting_3CatLine(columnName=columnName, titleName=titleName, folderName=folderName, yScale='log', yLabel=self.__yLabelPlotCrossSections)
            ######################## Scattering TMatrix
            columnName = 'SCA_TMatrix'
            titleName = "T-matrix Scattering Cross Section"
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
                              column1, column1T,
                              column2, column2T,
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
                for j in A2.unique():
                    fileNames = self.dfMainInfo[(A1Name == i) & (A2Name == j)].AA_FileName

                    indexBarPlot = np.arange(3)

                    RDG_Height = []
                    TMatrix_Height = []

                    for index, fileN in fileNames.iteritems():
                        if fileN in self.dictData:
                            df = self.dictTotalOpticalProp[fileN]

                            selectedRow = df.loc[(df['sigma'] == sigma) & (df['median'] == median)]

                            if (len(selectedRow[column1]) != 1) or (len(selectedRow[column2]) != 1):
                                raise
                            else:
                                for c1 in selectedRow[column1]:
                                    for c2 in selectedRow[column2]:
                                        RDG_Height.append(c1)
                                        TMatrix_Height.append(c2)

                    AXES_A = ax1[rowCount, colCount].bar(indexBarPlot - self.__barWidth / 2, RDG_Height, color='red', width=self.__barWidth, edgecolor='white', label=column1T)
                    AXES_B = ax1[rowCount, colCount].bar(indexBarPlot + self.__barWidth / 2, TMatrix_Height, color='blue', width=self.__barWidth, edgecolor='white', label=column2T)

                    if showValue:
                        autoLabel(self, rects=AXES_A, row=rowCount, col=colCount, format=valueFormat)
                        autoLabel(self, rects=AXES_B, row=rowCount, col=colCount, format=valueFormat)

                    ax1[rowCount, colCount].set_xticks(indexBarPlot)
                    ax1[rowCount, colCount].set_xticklabels(("$D_m=$" + str(2.2), "$D_m=$" + str(2.5), "$D_m=$" + str(2.8)), fontsize=self.__xLabelGroupFontSizeTotal)
                    ax1[rowCount, colCount].grid(True, which='both', axis="y", alpha=0.5)
                    ax1[rowCount, colCount].set_yscale(yScale)
                    ax1[rowCount, colCount].yaxis.set_tick_params(labelsize=self.__yMajorTickLabelGeneralFontSize * 1.3)
                    ax1[rowCount, colCount].yaxis.set_major_formatter(FormatStrFormatter(yAxisFormat))

                    if colCount == 0:
                        ax1[rowCount, colCount].set_ylabel("$\\rho_{eff,100}$= " + str(i) + " kg/m$^3$", fontsize=self.__yLabelMultipleSubplotFontSize + 2)
                        ax1[rowCount, colCount].yaxis.labelpad = 39
                    if rowCount == (len(A1.unique()) - 1):
                        ax1[rowCount, colCount].set_xlabel("$\sigma_p|d_m=$" + str(j), fontsize=self.__xLabelMultipleSubplotFontSize)
                        ax1[rowCount, colCount].xaxis.labelpad = 34

                    colCount += 1
                rowCount += 1

            F1 = f"{title} for median={median}-sigma={sigma}"
            T1 = title + titleAppend + " for " + "$d_{m,g}=$" + str(round(median)) + " nm " + "and " + "$\sigma_g=$" + str(round(sigma, 2))

            fig.subplots_adjust(top=0.93, wspace=0.017, hspace=0.018)
            fig.suptitle(T1, fontsize=self.__plotTitleFontSizeTotal, x=0.65)

            if 'e' in yAxisFormat:
                fig.text(0.5, 0.06, self.__xLabelMassMob, ha='center', fontsize=18)
                fig.text(0.051, 0.5, yCommonLabel, va='center', rotation='vertical', fontsize=18)
            elif 'f' in yAxisFormat:
                fig.text(0.5, 0.06, self.__xLabelMassMob, ha='center', fontsize=18)
                fig.text(0.062, 0.5, yCommonLabel, va='center', rotation='vertical', fontsize=18)

            leg = plt.legend(bbox_to_anchor=(1.05, 1.7), markerscale=self.__legendMarkerScale * 1.33, fontsize=25, loc='lower left', borderaxespad=0.)
            ################# set the lineWidth of each legend object
            for legObj in leg.legendHandles:
                legObj.set_linewidth(3.0)
            #################
            self.SaveAndClosePlot(folderName=folderName, F1=F1)


        except Exception as e:
            logging.exception(e)
            raise

    def PlotSettingTotalBars(self, colName1, colDetail1, colName2, colDetail2,
                             titleName, titleA, shareY, folderName,
                             yScale, yLabel, showValue, valueFormat, yAxisFormat):
        try:
            for sigma in self.__dmSigmaLogN:
                for median in self.__dmMedianLogN:
                    self.PlotMultipleBarsTotal(A1=self.ser_rhoEff100nm, A1Name=self.dfMainInfo.AGG_EFF_RHO_100NM_CENTER,
                                               A2=self.ser_SigmaMob, A2Name=self.dfMainInfo.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER,
                                               column1=colName1, column1T=colDetail1,
                                               column2=colName2, column2T=colDetail2,
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
                               column1, column1T,
                               column2, column2T,
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

                        if fileN in self.dictData:
                            df = self.dictData[fileN]

                            ax1[rowCount, colCount].plot(df['dm'], df[column1], label=column1T + ', ' + Dmlist[DmCount],
                                                         color='red', linewidth=self.__lineWidth[c_lineWidth] * 1.33,
                                                         linestyle=self.__lineStyle[c_lineStyle], alpha=alphaMainLine,
                                                         marker=self.__markerStyle[c_markerStyle], markersize=self.__markerSize * 1.30)
                            ax1[rowCount, colCount].plot(df['dm'], df[column2], label=column2T + ', ' + Dmlist[DmCount],
                                                         color='blue', linewidth=self.__lineWidth[c_lineWidth] * 1.33,
                                                         linestyle=self.__lineStyle[c_lineStyle], alpha=alphaMainLine,
                                                         marker=self.__markerStyle[c_markerStyle], markersize=self.__markerSize * 1.30)
                            c_lineStyle += 1
                            c_markerStyle += 1
                            c_lineWidth += 1
                            DmCount += 1

                    ax1[rowCount, colCount].grid(True, which='both', axis="both", alpha=0.5)
                    ax1[rowCount, colCount].set_xscale('log')
                    ax1[rowCount, colCount].set_yscale(yScale)
                    ax1[rowCount, colCount].xaxis.set_major_formatter(FormatStrFormatter("%i"))
                    ax1[rowCount, colCount].xaxis.set_major_locator(plt.MaxNLocator(6))
                    ax1[rowCount, colCount].tick_params(axis='x', which='both', labelsize=self.__xMajorTickLabelGeneralFontSize * 1.30, labelrotation=25)
                    ax1[rowCount, colCount].tick_params(axis='y', which='both', labelsize=self.__yMajorTickLabelGeneralFontSize * 1.30)
                    ax1[rowCount, colCount].yaxis.set_major_formatter(FormatStrFormatter(yAxisFormat))

                    if colCount == 0:
                        ax1[rowCount, colCount].set_ylabel("$\\rho_{eff,100}$= " + str(i) + " kg/m$^3$", fontsize=self.__yLabelMultipleSubplotFontSize + 2)
                        ax1[rowCount, colCount].yaxis.labelpad = 39

                    if rowCount == (len(A1.unique()) - 1):
                        ax1[rowCount, colCount].set_xlabel("$\sigma_p|d_m=$" + str(j), fontsize=self.__xLabelMultipleSubplotFontSize)
                        ax1[rowCount, colCount].xaxis.labelpad = 34

                    colCount += 1
                rowCount += 1

            F1 = f"{title}"
            T1 = title + titleAppend
            fig.subplots_adjust(top=0.93, wspace=0.017, hspace=0.018)
            fig.suptitle(T1, fontsize=self.__plotTitleFontSizeTotal, x=0.65)

            if 'e' in yAxisFormat:
                fig.text(0.5, 0.06, self.__xLabelMobDiameter, ha='center', fontsize=18)
                fig.text(0.051, 0.5, yCommonLabel, va='center', rotation='vertical', fontsize=18)
            elif 'f' in yAxisFormat:
                fig.text(0.5, 0.06, self.__xLabelMobDiameter, ha='center', fontsize=18)
                fig.text(0.062, 0.5, yCommonLabel, va='center', rotation='vertical', fontsize=18)

            leg = plt.legend(bbox_to_anchor=(1.05, 1.7), markerscale=self.__legendMarkerScale * 1.33, fontsize=25, loc='lower left', borderaxespad=0.)
            ################# set the lineWidth of each legend object
            for legObj in leg.legendHandles:
                legObj.set_linewidth(3.0)
            #################
            self.SaveAndClosePlot(folderName=folderName, F1=F1)

        except Exception as e:
            logging.exception(e)
            raise

    def PlotSettingTotalLines(self, colName1, colDetail1, colName2, colDetail2,
                              titleName, titleA, shareY, folderName,
                              yScale, yLabel, yAxisFormat):
        try:

            self.PlotMultipleLinesTotal(A1=self.ser_rhoEff100nm, A1Name=self.dfMainInfo.AGG_EFF_RHO_100NM_CENTER,
                                        A2=self.ser_SigmaMob, A2Name=self.dfMainInfo.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER,
                                        column1=colName1, column1T=colDetail1,
                                        column2=colName2, column2T=colDetail2,
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

                    AXES_A = ax1[rowCount, colCount].bar(indexBarPlot - self.__barRatioWidth, B[0], color='red', width=self.__barRatioWidth, edgecolor='white', label=A[0])
                    AXES_B = ax1[rowCount, colCount].bar(indexBarPlot, B[1], color='blue', width=self.__barRatioWidth, edgecolor='white', label=A[1])
                    AXES_C = ax1[rowCount, colCount].bar(indexBarPlot + self.__barRatioWidth, B[2], color='green', width=self.__barRatioWidth, edgecolor='white', label=A[2])
                    if showValue:
                        autoLabel(self, rects=AXES_A, row=rowCount, col=colCount, format=valueFormat)
                        autoLabel(self, rects=AXES_B, row=rowCount, col=colCount, format=valueFormat)
                        autoLabel(self, rects=AXES_C, row=rowCount, col=colCount, format=valueFormat)

                    ax1[rowCount, colCount].set_xticks(indexBarPlot)
                    ax1[rowCount, colCount].set_xticklabels(("$D_m=$" + str(2.2), "$D_m=$" + str(2.5), "$D_m=$" + str(2.8)), fontsize=self.__xLabelGroupFontSizeTotal)
                    ax1[rowCount, colCount].grid(True, which='both', axis="y", alpha=0.5)
                    ax1[rowCount, colCount].set_yscale(yScale)
                    ax1[rowCount, colCount].yaxis.set_tick_params(labelsize=self.__yMajorTickLabelGeneralFontSize * 1.3)
                    ax1[rowCount, colCount].yaxis.set_major_formatter(FormatStrFormatter(yAxisFormat))
                    # ax1[rowCount, colCount].yaxis.set_minor_formatter(FormatStrFormatter(yAxisFormat))

                    if colCount == 0:
                        ax1[rowCount, colCount].set_ylabel("$\\rho_{eff,100}$= " + str(i) + " kg/m$^3$", fontsize=self.__yLabelMultipleSubplotFontSize + 2)
                        ax1[rowCount, colCount].yaxis.labelpad = 39
                    if rowCount == (len(A1.unique()) - 1):
                        ax1[rowCount, colCount].set_xlabel("$\sigma_p|d_m=$" + str(j), fontsize=self.__xLabelMultipleSubplotFontSize)
                        ax1[rowCount, colCount].xaxis.labelpad = 34

                    colCount += 1
                rowCount += 1

            if median:
                F1 = f"{title} for median={median}"
                T1 = title + titleAppend + " for " + "$d_{m,g}=$" + str(round(median)) + " nm "
            if sigma:
                F1 = f"{title} for sigma={sigma}"
                T1 = title + titleAppend + " for " + "$\sigma_g=$" + str(round(sigma, 2))

            fig.subplots_adjust(top=0.93, wspace=0.017, hspace=0.018)
            fig.suptitle(T1, fontsize=self.__plotTitleFontSizeTotal, x=0.65)

            if 'e' in yAxisFormat:
                fig.text(0.5, 0.06, self.__xLabelMassMob, ha='center', fontsize=18)
                fig.text(0.051, 0.5, yCommonLabel, va='center', rotation='vertical', fontsize=18)
            elif 'f' in yAxisFormat:
                fig.text(0.5, 0.06, self.__xLabelMassMob, ha='center', fontsize=18)
                fig.text(0.062, 0.5, yCommonLabel, va='center', rotation='vertical', fontsize=18)

            leg = plt.legend(bbox_to_anchor=(1.05, 1.7), markerscale=self.__legendMarkerScale * 1.33, fontsize=25, loc='lower left', borderaxespad=0.)
            ################# set the lineWidth of each legend object
            for legObj in leg.legendHandles:
                legObj.set_linewidth(3.0)
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
            Address = GF.getAddressTo(FolderName=self.__folderNameGraph + f"\{folderName}", FileName=F1, Extension="jpg")
            plt.savefig(Address, format='jpg', dpi=self.__figureDPI, bbox_inches='tight')

            if self.__WMF_SVGSaving:
                Address = GF.getAddressTo(FolderName=self.__folderNameGraph + f"\{folderName}", FileName=F1, Extension="svg")
                plt.savefig(Address, format='svg', dpi=self.__figureDPI, bbox_inches='tight')

                AddressWMF = GF.getAddressTo(FolderName=self.__folderNameGraph + f"\{folderName}", FileName=F1, Extension="wmf")
                subprocess.call([self.__inkScapePath, str(Address.resolve()), '--export-wmf', str(AddressWMF.resolve())])

            plt.clf()
            plt.close()



        except Exception as e:
            logging.exception(e)
            raise

    ######################################################################
    ######################################################################
    ######################################################################
    ######################################################################

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
            # plt.savefig(Address, format='jpg', dpi=Figure_DPI, bbox_inches='tight')
            plt.clf()
            plt.close()

        except Exception as e:
            logging.exception(e)
            raise

    def Fig_Plot_Save_Scatterplot_Matrix(self, Address, Dataframe, Figure_DPI=1000):
        try:
            sns.set(style='white')
            g = sns.pairplot(Dataframe)
            g.map_lower(sns.kdeplot, cmap="GnBu", shade=True)
            g.map_diag(plt.hist, edgecolor="b")
            g.map_upper(self.corrDots)

            xlabels, ylabels = [], []
            for ax in g.axes[-1, :]:
                xlabel = ax.xaxis.get_label_text()
                xlabels.append(xlabel)
            for ax in g.axes[:, 0]:
                ylabel = ax.yaxis.get_label_text()
                ylabels.append(ylabel)
            for i in range(len(xlabels)):
                for j in range(len(ylabels)):
                    g.axes[j, i].xaxis.set_label_text(xlabels[i])
                    g.axes[j, i].yaxis.set_label_text(ylabels[j])

            plt.savefig(Address, format='jpg', dpi=Figure_DPI)
            plt.clf()
            plt.close()
        except Exception as e:
            logging.exception(e)
            raise

    def Fig_Plot_Save_Scatter_X_Linear_Y_Linear(self, Address, X_Array, Y_array, tickLabelStyle='sci', X_Min=None, X_Max=None, Y_Min=None, Y_Max=None, X_Label=None, Y_label=None, Plot_Title=None,
                                                label_font_size=12, Plot_Title_Size=12, Figure_DPI=1000, alpha_Y=0.3, Marker_Size=3):
        try:

            fig, ax1 = plt.subplots()
            plt.ticklabel_format(style=tickLabelStyle, axis='x', scilimits=(0, 0))
            plt.ticklabel_format(style=tickLabelStyle, axis='y', scilimits=(0, 0))

            if X_Min == None:
                X_Min = float(min(X_Array))
                X_Min = X_Min - (abs(X_Min) * 0.2)
            if X_Max == None:
                X_Max = float(max(X_Array))
                X_Max = X_Max + (abs(X_Max) * 0.2)
            if Y_Min == None:
                Y_Min = float(min(Y_array))
                Y_Min = Y_Min - (abs(Y_Min) * 0.2)
            if Y_Max == None:
                Y_Max = float(max(Y_array))
                Y_Max = Y_Max + (abs(Y_Max) * 0.2)

            ax1.scatter(X_Array, Y_array, s=7, alpha=alpha_Y)
            if X_Label != None:
                ax1.set_xlabel(X_Label, fontsize=label_font_size)
            if Y_label != None:
                ax1.set_ylabel(Y_label, fontsize=label_font_size)
            ax1.set_xlim(X_Min, X_Max)
            ax1.set_ylim(Y_Min, Y_Max)
            ax1.grid(True, which='major', axis="both", alpha=0.5)
            if Plot_Title != None:
                plt.title(Plot_Title, fontsize=Plot_Title_Size, y=1.0)
            plt.savefig(Address, format='jpg', dpi=Figure_DPI, bbox_inches='tight')
            plt.clf()
            plt.close()

        except Exception as e:
            logging.exception(e)
            raise

    def corrDots(self, *args, **kwargs):
        try:
            corr_r = args[0].corr(args[1], 'pearson')
            corr_text = f"{corr_r:2.2f}".replace("0.", ".")
            ax = plt.gca()
            ax.set_axis_off()
            marker_size = abs(corr_r) * 10000
            ax.scatter([.5], [.5], marker_size, [corr_r], alpha=0.95, cmap="coolwarm", vmin=-1, vmax=1, transform=ax.transAxes)
            font_size = abs(corr_r) * 40 + 5
            ax.annotate(corr_text, [.5, .5, ], xycoords="axes fraction", ha='center', va='center', fontsize=font_size)

        except Exception as e:
            logging.exception(e)
            raise
