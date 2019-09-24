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
from scipy.optimize import fsolve
import subprocess

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
            self.__dmSigmaLogN = [1.2, 1.4, 1.6]
            self.__dmMedianLogN = [100, 200, 300]
            self.__saveHistogram = False
            #################################################
            #################################################
            ################################################# Comparison Graphs
            self.__EXPFull = True
            self.__CRSGraphsExpAdded = True
            self.__ErrorAndRatioGraphsExpAdded = True
            self.__MACMSCGraphsExpAdded = True
            self.__EfficiencyGraphsExpAdded = True
            self.__SSAGraphsExpAdded = True
            self.__dpInfoGraphsExpAdded = True
            self.__dmDistributionGraphsBarExpAdded = True
            self.__dmDistributionGraphsLineExpAdded = True
            self.__dmDistributionRatioGraphsBarExpAdded = True

            #################################################
            #################################################
            ################################################# Experiment Graphs
            self.__EXPAccurate = True
            self.__CRSGraphsEXPAccurate = True
            self.__ErrorAndRatioGraphsEXPAccurate = True
            self.__MACMSCGraphsEXPAccurate = True
            self.__SSAGraphsEXPAccurate = True
            self.__dpMedianGraphsEXPAccurate = True

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
            self.__A1LineWidth = [1.8, 1.2, 0.6]
            self.__A1AlphaMainLine = 0.35
            self.__A2Color = 'blue'
            self.__A2LineWidth = [4, 3, 2.5]
            self.__A2AlphaMainLine = 0.65
            self.__subplotLineStyle = ['-', '--', '-.']
            self.__subplotMarkerStyle = ["o", "X", "^"]
            self.__subplotMarkerSize = 8
            self.__xMajorTickLabelSubplotFontSize = 15
            self.__xTickLabelSubplotRotation = 40
            self.__yMajorTickLabelSubplotFontSize = 15
            self.__yLabelEachSubplotFontSize = 26
            self.__yLabelEachSubplotRotation = -90
            self.__yLabelEachSubplotPad = 40
            self.__xLabelEachSubplotFontSize = 26
            self.__xLabelEachSubplotPad = 17
            self.__subplotGridSetup = [0.88, 0.018, 0.018]
            self.__xLabelCommonSubplotAdjustment = 0.033  # lower to get down
            self.__xLabelCommonSubplotFontSize = 24
            self.__yLabelCommonSubplotFontSize = 24
            self.__yLabelCommonSubplotAdjustment = [0.026, 0.044]  # lower to get left
            self.__legendSubplotAdjustment = [1.24, 1.7]
            self.__legendSubplotMarkerScale = 2.25
            self.__legendSubplotFontSize = 20
            self.__legendSubplotLineWidth = 4
            #################################################
            self.__figureDPI = 500
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
            self.__barWidth = 0.34
            self.__barRatioWidth = 0.25
            self.__plotTitleFontSizeTotal = 30
            self.__xLabelGroupFontSizeTotal = 18
            self.__xLabelCommonSubplotFontSize = 27

        except Exception as e:
            logging.exception(e)
            raise

    def ExecuteMobility(self, Dm, DmSTD,
                        rhoEff, rhoEffSTD,
                        temperature, temperatureSTD,
                        pressure, pressureSTD,
                        arraySize, randomSize,
                        fileName):
        try:
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

            def _calcMonteCarlo_Viscosity(length, temperature, pressure):
                try:
                    results = np.empty(shape=int(length))
                    for i in range(int(length)):
                        results[i] = calcViscosity(Temperature=temperature[i], Pressure=pressure[i])
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

            def __calcMobilityDiameter(aerodynamicDiameter_nm, meanFreePath_m, eff_Dm, eff_Rho_100nm):
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

            def _calcMonteCarlo_MobilityDiameter(length, aerodynamicDiameter_nm, meanFreePath_m, eff_Dm, eff_Rho_100nm):
                try:
                    results = np.empty(shape=int(length))
                    for i in range(int(length)):
                        results[i] = __calcMobilityDiameter(aerodynamicDiameter_nm=aerodynamicDiameter_nm[i],
                                                            meanFreePath_m=meanFreePath_m[i],
                                                            eff_Dm=eff_Dm[i],
                                                            eff_Rho_100nm=eff_Rho_100nm[i])
                    return results

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

            def __calcAbsorptionEff(absorptionCS, mobilityDiameter):
                try:

                    A = absorptionCS / (((mobilityDiameter / 1000) ** 2) * pi / 4)
                    return A

                except Exception as e:
                    logging.exception(e)
                    raise

            def _calcMonteCarlo_AbsorptionEff(length, absorptionCS, mobilityDiameter):
                try:
                    results = np.empty(shape=int(length))
                    for i in range(int(length)):
                        results[i] = __calcAbsorptionEff(absorptionCS=absorptionCS[i], mobilityDiameter=mobilityDiameter[i])
                    return results

                except Exception as e:
                    logging.exception(e)
                    raise

            ###################
            ###################

            def __calcScatteringEff(scatteringCS, mobilityDiameter):
                try:

                    A = scatteringCS / (((mobilityDiameter / 1000) ** 2) * pi / 4)
                    return A

                except Exception as e:
                    logging.exception(e)
                    raise

            def _calcMonteCarlo_ScatteringEff(length, scatteringCS, mobilityDiameter):
                try:
                    results = np.empty(shape=int(length))
                    for i in range(int(length)):
                        results[i] = __calcScatteringEff(scatteringCS=scatteringCS[i], mobilityDiameter=mobilityDiameter[i])
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
            dm_AVE = []
            dm_STD = []
            ABS_CRS_AVE = []
            ABS_CRS_STD = []
            SCA_CRS_AVE = []
            SCA_CRS_STD = []
            SSA_AVE = []
            SSA_STD = []
            ABS_Eff_AVE = []
            ABS_Eff_STD = []
            SCA_Eff_AVE = []
            SCA_Eff_STD = []
            ##################
            ##################
            for index, aerodynamicDiam in serAerodynamicDiameter.iteritems():

                arrAgg_Aerodynamic_Diameter = _createRandomNormalArr(center=serAerodynamicDiameter[index],
                                                                     width=ser_AAC_Diam_STD[index],
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

                arrPAX_Absorption_Coefficient = _createRandomNormalArr(center=ser_PAX_ABS_Coef_Ave[index],
                                                                       width=ser_PAX_ABS_Coef_STD[index],
                                                                       number=arraySize)

                arrPAX_Scattering_Coefficient = _createRandomNormalArr(center=ser_PAX_SCA_Coef_Ave[index],
                                                                       width=ser_PAX_SCA_Coef_STD[index],
                                                                       number=arraySize)

                arrCPC_Number_Density = _createRandomNormalArr(center=ser_CPC_Ave[index],
                                                               width=ser_CPC_STD[index],
                                                               number=arraySize)

                ######################################################

                arrAgg_Aerodynamic_Diameter_Random = _getRandomFromArr(array=arrAgg_Aerodynamic_Diameter, number=randomSize)
                arrAgg_Eff_Dm_Random = _getRandomFromArr(array=arrAgg_Eff_Dm, number=randomSize)
                arrAgg_Eff_Rho_100nm_Random = _getRandomFromArr(array=arrAgg_Eff_Rho_100nm, number=randomSize)
                arrAir_Temperature_Random = _getRandomFromArr(array=arrAir_Temperature, number=randomSize)
                arrAir_Pressure_Random = _getRandomFromArr(array=arrAir_Pressure, number=randomSize)
                arrPAX_Absorption_Coefficient_Random = _getRandomFromArr(array=arrPAX_Absorption_Coefficient, number=randomSize)
                arrPAX_Scattering_Coefficient_Random = _getRandomFromArr(array=arrPAX_Scattering_Coefficient, number=randomSize)
                arrCPC_Number_Density_Random = _getRandomFromArr(array=arrCPC_Number_Density, number=randomSize)

                ####################### Calculation

                arrMonteCarloViscosity = _calcMonteCarlo_Viscosity(length=randomSize,
                                                                   temperature=arrAir_Temperature_Random,
                                                                   pressure=arrAir_Pressure_Random)

                arrMonteCarloMeanFreePath = _calcMonteCarlo_MeanFreePath(length=randomSize,
                                                                         viscosity=arrMonteCarloViscosity,
                                                                         temperature=arrAir_Temperature_Random,
                                                                         pressure=arrAir_Pressure_Random)

                # arrViscosityAverageSTD = _getAverageAndSTD(arrMonteCarloViscosity)
                # arrMeanFreePathAverageSTD = _getAverageAndSTD(arrMonteCarloMeanFreePath)

                # Mobility diameter in nm
                arrMonteCarloMobilityDiameter = _calcMonteCarlo_MobilityDiameter(length=randomSize,
                                                                                 aerodynamicDiameter_nm=arrAgg_Aerodynamic_Diameter_Random,
                                                                                 meanFreePath_m=arrMonteCarloMeanFreePath,
                                                                                 eff_Dm=arrAgg_Eff_Dm_Random,
                                                                                 eff_Rho_100nm=arrAgg_Eff_Rho_100nm_Random)

                arrMobilityDiameterAverageSTD = _getAverageAndSTD(arrMonteCarloMobilityDiameter)

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
                arrMonteCarloAbsorptionEFF = _calcMonteCarlo_AbsorptionEff(length=randomSize,
                                                                           absorptionCS=arrMonteCarloAbsorptionCS,
                                                                           mobilityDiameter=arrMonteCarloMobilityDiameter)

                arrMonteCarloScatteringEFF = _calcMonteCarlo_ScatteringEff(length=randomSize,
                                                                           scatteringCS=arrMonteCarloScatteringCS,
                                                                           mobilityDiameter=arrMonteCarloMobilityDiameter)

                arrAbsorptionEffAverageSTD = _getAverageAndSTD(arrMonteCarloAbsorptionEFF)
                arrScatteringEffAverageSTD = _getAverageAndSTD(arrMonteCarloScatteringEFF)
                #######################

                dm_AVE.append(arrMobilityDiameterAverageSTD[0])
                dm_STD.append(arrMobilityDiameterAverageSTD[1])

                ABS_CRS_AVE.append(arrAbsorptionCSAverageSTD[0])
                ABS_CRS_STD.append(arrAbsorptionCSAverageSTD[1])

                SCA_CRS_AVE.append(arrScatteringCSAverageSTD[0])
                SCA_CRS_STD.append(arrScatteringCSAverageSTD[1])

                SSA_AVE.append(arrSSAAverageSTD[0])
                SSA_STD.append(arrSSAAverageSTD[1])

                ABS_Eff_AVE.append(arrAbsorptionEffAverageSTD[0])
                ABS_Eff_STD.append(arrAbsorptionEffAverageSTD[1])

                SCA_Eff_AVE.append(arrScatteringEffAverageSTD[0])
                SCA_Eff_STD.append(arrScatteringEffAverageSTD[1])
                ##########################
                if self.__saveHistogram:
                    ToSaveHistogram(self, folderName='ExpHistogram' + f"\{fileName}", fileName='mobilityDiameter' + f"_{str(round(aerodynamicDiam))}", array=arrMonteCarloMobilityDiameter)
                    ToSaveHistogram(self, folderName='ExpHistogram' + f"\{fileName}", fileName='absCRS' + f"_{str(round(aerodynamicDiam))}", array=arrMonteCarloAbsorptionCS)
                    ToSaveHistogram(self, folderName='ExpHistogram' + f"\{fileName}", fileName='scaCRS' + f"_{str(round(aerodynamicDiam))}", array=arrMonteCarloScatteringCS)
                    ToSaveHistogram(self, folderName='ExpHistogram' + f"\{fileName}", fileName='ssa' + f"_{str(round(aerodynamicDiam))}", array=arrMonteCarloSSA)
                    ToSaveHistogram(self, folderName='ExpHistogram' + f"\{fileName}", fileName='absEFF' + f"_{str(round(aerodynamicDiam))}", array=arrMonteCarloAbsorptionEFF)
                    ToSaveHistogram(self, folderName='ExpHistogram' + f"\{fileName}", fileName='scaEFF' + f"_{str(round(aerodynamicDiam))}", array=arrMonteCarloScatteringEFF)

            ExpDict = {'dm_AVE': dm_AVE,
                       'dm_STD': dm_STD,
                       'ABS_CRS_AVE': ABS_CRS_AVE,
                       'ABS_CRS_STD': ABS_CRS_STD,
                       'SCA_CRS_AVE': SCA_CRS_AVE,
                       'SCA_CRS_STD': SCA_CRS_STD,
                       'SSA_AVE': SSA_AVE,
                       'SSA_STD': SSA_STD,
                       'ABS_Eff_AVE': ABS_Eff_AVE,
                       'ABS_Eff_STD': ABS_Eff_STD,
                       'SCA_Eff_AVE': SCA_Eff_AVE,
                       'SCA_Eff_STD': SCA_Eff_STD,
                       }

            dfExperiment = pd.DataFrame(ExpDict)
            return dfExperiment

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
            ##################################
            self.__dpDeviationPercent = 25
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

    def GetDataLabel(self, fileName, index):
        try:
            self.dictDataLabel[fileName] = "$D_m=$" + str(round(self.ser_DmAve.loc[index], 2)) + ", " + \
                                           "$\\rho_{eff,100}$=" + str(round(self.ser_rhoEff100nmAve.loc[index], 1)) + ", " + \
                                           "$\sigma_p|d_m=$" + str(round(self.ser_SigmaMob.loc[index], 2)) + ", " + \
                                           "$D_{TEM}=$" + str(abs(round(self.ser_D_TEM.loc[index], 2))) + ", " + \
                                           "$d_{p,100}=$" + str(round(self.ser_dp100_nano.loc[index], 1))

        except Exception as e:
            logging.exception(e)
            raise

    def Experiment_RDG_TMatrixComparisonGraphs(self):
        try:
            if self.__EXPAccurate:

                self.dfMainInfoAcu = pd.read_csv(GF.GetAddressTo(folderName=self.__folderNameData, fileName='beaconACU', extension='csv'))
                self.dictData = {}
                self.dictDataLabel = {}
                self.dictTotalOpticalProp = {}
                self.dfEXP = {}
                self.GetDataDefinitionEXP(mainDf=self.dfMainInfoAcu)

                for index, fileN in self.serMainFileName.iteritems():
                    self.dfEXP[fileN] = self.ExecuteMobility(Dm=self.ser_DmAve.loc[index], DmSTD=self.ser_DmSTD.loc[index],
                                                             rhoEff=self.ser_rhoEff100nmAve.loc[index], rhoEffSTD=self.ser_rhoEff100nmSTD.loc[index],
                                                             temperature=self.serTemperatureAve.loc[index], temperatureSTD=self.serTemperatureSTD.loc[index],
                                                             pressure=self.serPressureAve.loc[index], pressureSTD=self.serPressureSTD.loc[index],
                                                             arraySize=self.serArraySize.loc[index], randomSize=self.serRandomSize.loc[index],
                                                             fileName=fileN)

                    self.dictData[fileN] = pd.read_csv(GF.GetAddressTo(folderName=self.__folderNameData, fileName=fileN))
                    self.dictData[fileN] = self.dictData[fileN].replace(0, np.nan)

                    if self.__dpDeviation:
                        self.Check_dpRatio(fileName=fileN)

                    self.CalcTotalMACMSC(self.dictData[fileN], fileN)

                    self.CalcEfficiency(fileName=fileN)

                    self.GetDataLabel(fileName=fileN, index=index)

                if self.__CRSGraphsEXPAccurate:
                    pass
                if self.__ErrorAndRatioGraphsEXPAccurate:
                    pass
                if self.__SSAGraphsEXPAccurate:
                    pass
                if self.__MACMSCGraphsEXPAccurate:
                    pass
                if self.__dpMedianGraphsEXPAccurate:
                    pass

            ###########################################
            ###########################################
            ###########################################
            ###########################################
            if self.__EXPFull:
                self.dfMainInfo = pd.read_csv(GF.GetAddressTo(folderName=self.__folderNameData, fileName='beacon', extension='csv'))
                self.dictData = {}
                self.dictDataLabel = {}
                self.dictTotalOpticalProp = {}
                self.dfEXP = {}

                self.GetDataDefinitionEXP(mainDf=self.dfMainInfo)

                for index, fileN in self.serMainFileName.iteritems():
                    self.dfEXP[fileN] = self.ExecuteMobility(Dm=self.ser_DmAve.loc[index], DmSTD=self.ser_DmSTD.loc[index],
                                                             rhoEff=self.ser_rhoEff100nmAve.loc[index], rhoEffSTD=self.ser_rhoEff100nmSTD.loc[index],
                                                             temperature=self.serTemperatureAve.loc[index], temperatureSTD=self.serTemperatureSTD.loc[index],
                                                             pressure=self.serPressureAve.loc[index], pressureSTD=self.serPressureSTD.loc[index],
                                                             arraySize=self.serArraySize.loc[index], randomSize=self.serRandomSize.loc[index],
                                                             fileName=fileN)

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
                Address = GF.GetAddressTo(folderName=self.__folderNameGraph + "\SSAs_EXP", fileName=F1, extension="jpg")
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
                Address = GF.GetAddressTo(folderName=self.__folderNameGraph + "\CrossSections_EXP", fileName=F1, extension="jpg")
                plt.savefig(Address, format='jpg', dpi=self.__figureDPI, bbox_inches='tight')
                plt.clf()
                plt.close()

        except Exception as e:
            logging.exception(e)
            raise

    def PlotCrossSectionsIntermediate_EXP(self, clName, tlName, mean, STD):
        try:
            if self.__rhoEff100nmCTE_EXPAccurate == True:
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

    def PlotTotalRatioGraphsBar(self):
        try:
            folderName = 'RatioBarGraphs'
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

    def PlotTotalGraphsLine(self):
        try:
            folderName = 'TotalLineGraphs'
            ######################## ABS Cross Section
            columnName1 = 'ABS_RDG'
            columnDetail1 = 'RDG-FA'
            columnName2 = 'ABS_TMatrix'
            columnDetail2 = 'T-matrix'
            titleName = "Absorption Cross-section"
            titleAppend = " Vs. Mobility-equivalent Diameter"
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
            titleName = "Scattering Cross-section"
            titleAppend = " Vs. Mobility-equivalent Diameter"
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
            titleAppend = " Vs. Mobility-equivalent Diameter"
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
            titleAppend = " Vs. Mobility-equivalent Diameter"
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
            titleName = "Single-scattering Albedo (SSA)"
            titleAppend = " Vs. Mobility-equivalent Diameter"
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
            titleAppend = " Vs. Mobility-equivalent Diameter"
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
            titleAppend = " Vs. Mobility-equivalent Diameter"
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
            titleAppend = " Vs. Mobility-equivalent Diameter"
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
            titleAppend = " Vs. Mobility-equivalent Diameter"
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
            titleName = "Total Absorption Cross-section"
            titleAppend = ""
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
            titleName = "Total Scattering Cross-section"
            titleAppend = ""
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
            titleAppend = ""
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
            titleAppend = ""
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
            folderName = 'OpticalCrossSections'
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

                    AXES_A = ax1[rowCount, colCount].bar(indexBarPlot - self.__barWidth / 2, RDG_Height, color=self.__A1Color, width=self.__barWidth, edgecolor='white', hatch=self.__patterns[0],
                                                         label=column1T)
                    AXES_B = ax1[rowCount, colCount].bar(indexBarPlot + self.__barWidth / 2, TMatrix_Height, color=self.__A2Color, width=self.__barWidth, edgecolor='white', hatch=self.__patterns[4],
                                                         label=column2T)

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
                                                         color=self.__A1Color, linewidth=self.__A1LineWidth[c_lineWidth],
                                                         linestyle=self.__subplotLineStyle[c_lineStyle], alpha=self.__A1AlphaMainLine,
                                                         marker=self.__subplotMarkerStyle[c_markerStyle], markersize=self.__subplotMarkerSize)
                            ax1[rowCount, colCount].plot(df['dm'], df[column2], label=column2T + ', ' + Dmlist[DmCount],
                                                         color=self.__A2Color, linewidth=self.__A2LineWidth[c_lineWidth],
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

    ######################################################################
    ######################################################################
    ######################################################################
    ######################################################################
