from ConfigParserM import logging
import GeneralFunctions as GF
import numpy as np
from matplotlib import rcParams
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.ticker import FormatStrFormatter
from math import pi
from math import exp
from scipy.optimize import fsolve

####### Plotting Parameters
rcParams['mathtext.fontset'] = 'stix'
rcParams['font.family'] = 'STIXGeneral'


##############

class GraphTools:
    def __init__(self, FolderInfo):
        try:
            self.__CRSGraphs = True
            self.__ErrorAndRatioGraphs = True
            self.__MACMSCGraphs = True
            self.__SSAGraphs = True
            self.__dpMedianGraphs = True
            #################################################
            self.__DmCTE = True
            self.__rhoEff100nmCTE = True
            self.__sigmaMobCTE = True
            #################################################
            self.__CRSGraphs_EXP = True
            self.__ErrorAndRatioGraphs_EXP = True
            self.__MACMSCGraphs_EXP = True
            self.__SSAGraphs_EXP = True
            self.__dpMedianGraphs_EXP = True
            #################################################
            self.__DmCTE_EXP = True
            self.__rhoEff100nmCTE_EXP = True
            self.__sigmaMobCTE_EXP = True
            #################################################
            self.__folderNameGraph = FolderInfo['FOLDER_NAME_GRAPH']
            self.__folderNameData = FolderInfo['FOLDER_NAME_DATA']
            self.__xLabelPlotCrossSections = "Mobility-equivalent Diameter (nm)"
            self.__yLabelPlotCrossSections = "Cross Section (µm" + "$^{}$".format(2) + ")"
            self.__yLabelErrorAndRatio1 = "Error Percent (%)"
            self.__yLabelErrorAndRatio2 = "Cross Section Ratio (TMatrix/RDG)"
            self.__yLabelSSA = "SSA"
            self.__yLabelMAC = "MAC (m$^2$/g)"
            self.__yLabelMSC = "MSC (m$^2$/g)"
            self.__yLabeldpMedian = "Primary Median Diameter (nm)"
            self.__yLabelCalcNumber = "#"
            #################################################
            self.__plotTitle_Size = 12
            self.__labelFont_size = 12
            self.__figureDPI = 750
            self.__markerSize = 2
            self.__alphaMainLine = 0.5
            self.__lineColor = ['red', 'blue', 'green']
            self.__lineStyle = ['-', '--', ':']
            self.__markerStyle = ["o", "X", "^"]
            self.__lineWidth = [1.5, 1, 0.5]
            self.__mobilityDiamLimits = [82, 850]

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

            B = 3

            self.MainDF = pd.read_csv(GF.getAddressTo(FolderName=self.__folderNameData, FileName='beacon', Extension='csv'))
            self.full_data_src_dict = {}
            self.full_data_label_dict = {}
            self.mainFileNames = self.MainDF.loc[self.MainDF['AA_Plot'] == 1, 'AA_FileName']
            self.DmSeries = self.MainDF.loc[self.MainDF['AA_Plot'] == 1, 'AGG_EFF_DM_CENTER']
            self.rhoEff100nmSeries = self.MainDF.loc[self.MainDF['AA_Plot'] == 1, 'AGG_EFF_RHO_100NM_CENTER']
            self.sigmaMobSeries = self.MainDF.loc[self.MainDF['AA_Plot'] == 1, 'AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER']
            self.D_TEMSeries = self.MainDF.loc[self.MainDF['AA_Plot'] == 1, 'D_TEM']
            self.dp100_nanoSeries = self.MainDF.loc[self.MainDF['AA_Plot'] == 1, 'dp100_nano']

            for index, item in self.mainFileNames.iteritems():
                self.full_data_src_dict[item] = pd.read_csv(GF.getAddressTo(FolderName=self.__folderNameData, FileName=item))
                self.full_data_src_dict[item] = self.full_data_src_dict[item].replace(0, np.nan)
                self.full_data_label_dict[item] = "$D_m=$" + str(round(self.DmSeries.loc[index], 2)) + ", " + "$\\rho_{eff,100}$=" + str(round(self.rhoEff100nmSeries.loc[index], 1)) \
                                                  + ", " + "$\sigma_p|d_m=$" + str(round(self.sigmaMobSeries.loc[index], 2)) + ", " + "$D_{TEM}=$" + str(round(self.D_TEMSeries.loc[index], 2)) \
                                                  + ", " + "$d_{p,100}=$" + str(round(self.dp100_nanoSeries.loc[index], 1))

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

                    fileNames = self.MainDF[(A1Name == i) & (A2Name == j)].AA_FileName

                    c_lineStyle = 0
                    c_markerStyle = 0
                    c_lineColor = 0

                    for index, item in fileNames.iteritems():
                        if item in self.full_data_src_dict:
                            df = self.full_data_src_dict[item]
                            ax1.plot(df['dm'], df[column], label=self.full_data_label_dict[item], color=self.__lineColor[c_lineColor], linewidth=self.__lineWidth[c_lineWidth],
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
                ax1.set_xlabel(self.__xLabelPlotCrossSections, fontsize=self.__labelFont_size)
                ax1.set_ylabel(yLabel, fontsize=self.__labelFont_size)
                ax1.set_xlim(self.__mobilityDiamLimits[0], self.__mobilityDiamLimits[1])

                ax1.set_title(T1, fontsize=self.__plotTitle_Size)
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
                self.PlotSSACore_EXP(A1=self.rhoEff100nmSeries, A1Name=self.MainDF.AGG_EFF_RHO_100NM_CENTER, A2=self.DmSeries, A2Name=self.MainDF.AGG_EFF_DM_CENTER,
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

                    fileNames = self.MainDF[(A1Name == i) & (A2Name == j)].AA_FileName

                    c_lineStyle = 0
                    c_markerStyle = 0
                    c_lineColor = 0

                    for index, item in fileNames.iteritems():
                        if item in self.full_data_src_dict:
                            df = self.full_data_src_dict[item]
                            ax1.plot(df['dm'], df[column], label=self.full_data_label_dict[item], color=self.__lineColor[c_lineColor], linewidth=self.__lineWidth[c_lineWidth],
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
                ax1.set_xlabel(self.__xLabelPlotCrossSections, fontsize=self.__labelFont_size)
                ax1.set_ylabel(self.__yLabelPlotCrossSections, fontsize=self.__labelFont_size)
                ax1.set_xlim(self.__mobilityDiamLimits[0], self.__mobilityDiamLimits[1])

                ax1.set_title(T1, fontsize=self.__plotTitle_Size)
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
                self.PlotCrossSectionsCore_EXP(A1=self.rhoEff100nmSeries, A1Name=self.MainDF.AGG_EFF_RHO_100NM_CENTER, A2=self.DmSeries, A2Name=self.MainDF.AGG_EFF_DM_CENTER,
                                               column=clName, mean=mean, STD=STD, title=tlName, T1A=" for " + "$\\rho_{eff,100}$= ", T1I=1, T1B=" kg/m$^3$", F1A=" for " + "rho_eff=", F1I=0)

            # if self.__DmCTE == True:
            #     self.PlotCrossSectionsCore(A1=self.DmSeries, A1Name=self.MainDF.AGG_EFF_DM_CENTER, A2=self.sigmaMobSeries, A2Name=self.MainDF.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER,
            #                                column=clName, title=tlName, T1A=" for " + "$D_m=$", T1I=2, F1A=" for " + "Dm=", F1I=2)
            #
            # if self.__sigmaMobCTE == True:
            #     self.PlotCrossSectionsCore(A1=self.sigmaMobSeries, A1Name=self.MainDF.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER, A2=self.rhoEff100nmSeries,
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
            a = 3

        except Exception as e:
            logging.exception(e)
            raise

    ###################################
    ###################################
    ###################################
    ###################################
    def RDG_TMatrixComparisonGraphs(self):
        try:
            self.MainDF = pd.read_csv(GF.getAddressTo(FolderName=self.__folderNameData, FileName='beacon', Extension='csv'))
            self.full_data_src_dict = {}
            self.full_data_label_dict = {}
            self.mainFileNames = self.MainDF.loc[self.MainDF['AA_Plot'] == 1, 'AA_FileName']
            self.DmSeries = self.MainDF.loc[self.MainDF['AA_Plot'] == 1, 'AGG_EFF_DM_CENTER']
            self.rhoEff100nmSeries = self.MainDF.loc[self.MainDF['AA_Plot'] == 1, 'AGG_EFF_RHO_100NM_CENTER']
            self.sigmaMobSeries = self.MainDF.loc[self.MainDF['AA_Plot'] == 1, 'AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER']
            self.D_TEMSeries = self.MainDF.loc[self.MainDF['AA_Plot'] == 1, 'D_TEM']
            self.dp100_nanoSeries = self.MainDF.loc[self.MainDF['AA_Plot'] == 1, 'dp100_nano']

            for index, item in self.mainFileNames.iteritems():
                self.full_data_src_dict[item] = pd.read_csv(GF.getAddressTo(FolderName=self.__folderNameData, FileName=item))
                self.full_data_src_dict[item] = self.full_data_src_dict[item].replace(0, np.nan)
                self.full_data_label_dict[item] = "$D_m=$" + str(round(self.DmSeries.loc[index], 2)) + ", " + "$\\rho_{eff,100}$=" + str(round(self.rhoEff100nmSeries.loc[index], 1)) \
                                                  + ", " + "$\sigma_p|d_m=$" + str(round(self.sigmaMobSeries.loc[index], 2)) + ", " + "$D_{TEM}=$" + str(round(self.D_TEMSeries.loc[index], 2)) \
                                                  + ", " + "$d_{p,100}=$" + str(round(self.dp100_nanoSeries.loc[index], 1))

            if self.__CRSGraphs:
                self.PlotCrossSections()
            if self.__ErrorAndRatioGraphs:
                self.PlotErrorAndRatio()
            if self.__SSAGraphs:
                self.PlotSSA()
            if self.__MACMSCGraphs:
                self.PlotMACMSC()
            if self.__dpMedianGraphs:
                self.PlotdpMedian()

        except Exception as e:
            logging.exception(e)
            raise

    ###################################
    ###################################
    ###################################
    ###################################
    def PlotdpMedianCore(self, A1, A1Name, A2, A2Name, T1A, T1I, F1A, F1I, title, column, yLabel, T1B=None, F1B=None):
        try:
            alphaMainLine = self.__alphaMainLine

            for i in A1.unique():

                fig, ax1 = plt.subplots()

                c_lineColor = 0
                c_lineWidth = 0

                for j in A2.unique():

                    fileNames = self.MainDF[(A1Name == i) & (A2Name == j)].AA_FileName

                    c_lineStyle = 0
                    c_markerStyle = 0

                    for index, item in fileNames.iteritems():
                        if item in self.full_data_src_dict:
                            df = self.full_data_src_dict[item]
                            ax1.plot(df['dm'], df[column], label=self.full_data_label_dict[item], color=self.__lineColor[c_lineColor], linewidth=self.__lineWidth[c_lineWidth],
                                     linestyle=self.__lineStyle[c_lineStyle], alpha=alphaMainLine, marker=self.__markerStyle[c_markerStyle], markersize=self.__markerSize)

                        c_lineStyle += 1
                        c_markerStyle += 1

                    c_lineColor += 1
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
                ax1.set_xlabel(self.__xLabelPlotCrossSections, fontsize=self.__labelFont_size)
                ax1.set_ylabel(yLabel, fontsize=self.__labelFont_size)
                ax1.set_xlim(self.__mobilityDiamLimits[0], self.__mobilityDiamLimits[1])

                ax1.set_title(T1, fontsize=self.__plotTitle_Size)
                ax1.xaxis.set_major_formatter(FormatStrFormatter("%i"))
                ax1.xaxis.set_minor_formatter(FormatStrFormatter("%i"))
                ax1.legend(bbox_to_anchor=(1.00, 0.5), loc='center left', fontsize='large')
                Address = GF.getAddressTo(FolderName=self.__folderNameGraph + "\PrimaryParticle", FileName=F1, Extension="jpg")
                plt.savefig(Address, format='jpg', dpi=self.__figureDPI, bbox_inches='tight')
                plt.clf()
                plt.close()

        except Exception as e:
            logging.exception(e)
            raise

    def PlotdpMedianIntermediate(self, clName, tlName, yLabel):
        try:
            if self.__rhoEff100nmCTE == True:
                self.PlotdpMedianCore(A1=self.rhoEff100nmSeries, A1Name=self.MainDF.AGG_EFF_RHO_100NM_CENTER, A2=self.DmSeries, A2Name=self.MainDF.AGG_EFF_DM_CENTER,
                                      column=clName, yLabel=yLabel, title=tlName, T1A=" for " + "$\\rho_{eff,100}$= ", T1I=1, T1B=" kg/m$^3$", F1A=" for " + "rho_eff=", F1I=0)

            if self.__DmCTE == True:
                self.PlotdpMedianCore(A1=self.DmSeries, A1Name=self.MainDF.AGG_EFF_DM_CENTER, A2=self.sigmaMobSeries, A2Name=self.MainDF.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER,
                                      column=clName, yLabel=yLabel, title=tlName, T1A=" for " + "$D_m=$", T1I=2, F1A=" for " + "Dm=", F1I=2)

            if self.__sigmaMobCTE == True:
                self.PlotdpMedianCore(A1=self.sigmaMobSeries, A1Name=self.MainDF.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER, A2=self.rhoEff100nmSeries,
                                      A2Name=self.MainDF.AGG_EFF_RHO_100NM_CENTER,
                                      column=clName, yLabel=yLabel, title=tlName, T1A=" for " + "$\sigma_p|d_m=$", T1I=2, F1A=" for " + "sigma_p=", F1I=2)

        except Exception as e:
            logging.exception(e)
            raise

    def PlotdpMedian(self):
        try:
            ######################## MAC RDG
            clName = 'dp_median'
            tlName = "Primary Median Diameter for each Mobility Diameter"
            self.PlotdpMedianIntermediate(clName=clName, tlName=tlName, yLabel=self.__yLabeldpMedian)
            ######################## MSC RDG
            clName = 'NumberOfCalcs'
            tlName = "Number of Calculations for each Mobility Diameter"
            self.PlotdpMedianIntermediate(clName=clName, tlName=tlName, yLabel=self.__yLabelCalcNumber)

        except Exception as e:
            logging.exception(e)
            raise

    ###################################
    ###################################
    ###################################
    ###################################
    def PlotMACMSCCore(self, A1, A1Name, A2, A2Name, T1A, T1I, F1A, F1I, title, column, yLabel, T1B=None, F1B=None):
        try:
            alphaMainLine = self.__alphaMainLine

            for i in A1.unique():

                fig, ax1 = plt.subplots()

                c_lineColor = 0
                c_lineWidth = 0

                for j in A2.unique():

                    fileNames = self.MainDF[(A1Name == i) & (A2Name == j)].AA_FileName

                    c_lineStyle = 0
                    c_markerStyle = 0

                    for index, item in fileNames.iteritems():
                        if item in self.full_data_src_dict:
                            df = self.full_data_src_dict[item]
                            ax1.plot(df['dm'], df[column], label=self.full_data_label_dict[item], color=self.__lineColor[c_lineColor], linewidth=self.__lineWidth[c_lineWidth],
                                     linestyle=self.__lineStyle[c_lineStyle], alpha=alphaMainLine, marker=self.__markerStyle[c_markerStyle], markersize=self.__markerSize)

                        c_lineStyle += 1
                        c_markerStyle += 1

                    c_lineColor += 1
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
                ax1.set_xlabel(self.__xLabelPlotCrossSections, fontsize=self.__labelFont_size)
                ax1.set_ylabel(yLabel, fontsize=self.__labelFont_size)
                ax1.set_xlim(self.__mobilityDiamLimits[0], self.__mobilityDiamLimits[1])

                ax1.set_title(T1, fontsize=self.__plotTitle_Size)
                ax1.xaxis.set_major_formatter(FormatStrFormatter("%i"))
                ax1.xaxis.set_minor_formatter(FormatStrFormatter("%i"))
                ax1.legend(bbox_to_anchor=(1.00, 0.5), loc='center left', fontsize='large')
                Address = GF.getAddressTo(FolderName=self.__folderNameGraph + "\MAC_MSC", FileName=F1, Extension="jpg")
                plt.savefig(Address, format='jpg', dpi=self.__figureDPI, bbox_inches='tight')
                plt.clf()
                plt.close()

        except Exception as e:
            logging.exception(e)
            raise

    def PlotMACMSCIntermediate(self, clName, tlName, yLabel):
        try:
            if self.__rhoEff100nmCTE == True:
                self.PlotMACMSCCore(A1=self.rhoEff100nmSeries, A1Name=self.MainDF.AGG_EFF_RHO_100NM_CENTER, A2=self.DmSeries, A2Name=self.MainDF.AGG_EFF_DM_CENTER,
                                    column=clName, yLabel=yLabel, title=tlName, T1A=" for " + "$\\rho_{eff,100}$= ", T1I=1, T1B=" kg/m$^3$", F1A=" for " + "rho_eff=", F1I=0)

            if self.__DmCTE == True:
                self.PlotMACMSCCore(A1=self.DmSeries, A1Name=self.MainDF.AGG_EFF_DM_CENTER, A2=self.sigmaMobSeries, A2Name=self.MainDF.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER,
                                    column=clName, yLabel=yLabel, title=tlName, T1A=" for " + "$D_m=$", T1I=2, F1A=" for " + "Dm=", F1I=2)

            if self.__sigmaMobCTE == True:
                self.PlotMACMSCCore(A1=self.sigmaMobSeries, A1Name=self.MainDF.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER, A2=self.rhoEff100nmSeries,
                                    A2Name=self.MainDF.AGG_EFF_RHO_100NM_CENTER,
                                    column=clName, yLabel=yLabel, title=tlName, T1A=" for " + "$\sigma_p|d_m=$", T1I=2, F1A=" for " + "sigma_p=", F1I=2)

        except Exception as e:
            logging.exception(e)
            raise

    def PlotMACMSC(self):
        try:
            ######################## MAC RDG
            clName = 'MAC_RDG'
            tlName = "MAC for RDG-FA"
            self.PlotMACMSCIntermediate(clName=clName, tlName=tlName, yLabel=self.__yLabelMAC)
            ######################## MSC RDG
            clName = 'MSC_RDG'
            tlName = "MSC for RDG-FA"
            self.PlotMACMSCIntermediate(clName=clName, tlName=tlName, yLabel=self.__yLabelMSC)
            ######################## MAC TMatrix
            clName = 'MAC_TMatrix'
            tlName = "MAC for TMatrix"
            self.PlotMACMSCIntermediate(clName=clName, tlName=tlName, yLabel=self.__yLabelMAC)
            ######################## MSC TMatrix
            clName = 'MSC_TMatrix'
            tlName = "MSC for TMatrix"
            self.PlotMACMSCIntermediate(clName=clName, tlName=tlName, yLabel=self.__yLabelMSC)



        except Exception as e:
            logging.exception(e)
            raise

    ###################################
    ###################################
    ###################################
    ###################################
    def PlotSSACore(self, A1, A1Name, A2, A2Name, T1A, T1I, F1A, F1I, title, column, yLabel, T1B=None, F1B=None):
        try:
            alphaMainLine = self.__alphaMainLine

            for i in A1.unique():

                fig, ax1 = plt.subplots()

                c_lineColor = 0
                c_lineWidth = 0

                for j in A2.unique():

                    fileNames = self.MainDF[(A1Name == i) & (A2Name == j)].AA_FileName

                    c_lineStyle = 0
                    c_markerStyle = 0

                    for index, item in fileNames.iteritems():
                        if item in self.full_data_src_dict:
                            df = self.full_data_src_dict[item]
                            ax1.plot(df['dm'], df[column], label=self.full_data_label_dict[item], color=self.__lineColor[c_lineColor], linewidth=self.__lineWidth[c_lineWidth],
                                     linestyle=self.__lineStyle[c_lineStyle], alpha=alphaMainLine, marker=self.__markerStyle[c_markerStyle], markersize=self.__markerSize)

                        c_lineStyle += 1
                        c_markerStyle += 1

                    c_lineColor += 1
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
                ax1.set_xlabel(self.__xLabelPlotCrossSections, fontsize=self.__labelFont_size)
                ax1.set_ylabel(yLabel, fontsize=self.__labelFont_size)
                ax1.set_xlim(self.__mobilityDiamLimits[0], self.__mobilityDiamLimits[1])

                ax1.set_title(T1, fontsize=self.__plotTitle_Size)
                ax1.xaxis.set_major_formatter(FormatStrFormatter("%i"))
                ax1.xaxis.set_minor_formatter(FormatStrFormatter("%i"))
                ax1.legend(bbox_to_anchor=(1.00, 0.5), loc='center left', fontsize='large')
                Address = GF.getAddressTo(FolderName=self.__folderNameGraph + "\SSAs", FileName=F1, Extension="jpg")
                plt.savefig(Address, format='jpg', dpi=self.__figureDPI, bbox_inches='tight')
                plt.clf()
                plt.close()

        except Exception as e:
            logging.exception(e)
            raise

    def PlotSSAIntermediate(self, clName, tlName, yLabel):
        try:
            if self.__rhoEff100nmCTE == True:
                self.PlotSSACore(A1=self.rhoEff100nmSeries, A1Name=self.MainDF.AGG_EFF_RHO_100NM_CENTER, A2=self.DmSeries, A2Name=self.MainDF.AGG_EFF_DM_CENTER,
                                 column=clName, yLabel=yLabel, title=tlName, T1A=" for " + "$\\rho_{eff,100}$= ", T1I=1, T1B=" kg/m$^3$", F1A=" for " + "rho_eff=", F1I=0)

            if self.__DmCTE == True:
                self.PlotSSACore(A1=self.DmSeries, A1Name=self.MainDF.AGG_EFF_DM_CENTER, A2=self.sigmaMobSeries, A2Name=self.MainDF.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER,
                                 column=clName, yLabel=yLabel, title=tlName, T1A=" for " + "$D_m=$", T1I=2, F1A=" for " + "Dm=", F1I=2)

            if self.__sigmaMobCTE == True:
                self.PlotSSACore(A1=self.sigmaMobSeries, A1Name=self.MainDF.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER, A2=self.rhoEff100nmSeries,
                                 A2Name=self.MainDF.AGG_EFF_RHO_100NM_CENTER,
                                 column=clName, yLabel=yLabel, title=tlName, T1A=" for " + "$\sigma_p|d_m=$", T1I=2, F1A=" for " + "sigma_p=", F1I=2)

        except Exception as e:
            logging.exception(e)
            raise

    def PlotSSA(self):
        try:
            ######################## SSA TMatrix
            clName = 'SSA_TMatrix'
            tlName = "SSA for TMatrix"
            self.PlotSSAIntermediate(clName=clName, tlName=tlName, yLabel=self.__yLabelSSA)

            ######################## SSA RDG
            clName = 'SSA_RDG'
            tlName = "SSA for RDG-FA"
            self.PlotSSAIntermediate(clName=clName, tlName=tlName, yLabel=self.__yLabelSSA)



        except Exception as e:
            logging.exception(e)
            raise

    ###################################
    ###################################
    ###################################
    ###################################
    def PlotErrorAndRatioCore(self, A1, A1Name, A2, A2Name, T1A, T1I, F1A, F1I, title, column, yLabel, T1B=None, F1B=None):
        try:
            alphaMainLine = self.__alphaMainLine

            for i in A1.unique():

                fig, ax1 = plt.subplots()

                c_lineColor = 0
                c_lineWidth = 0

                for j in A2.unique():

                    fileNames = self.MainDF[(A1Name == i) & (A2Name == j)].AA_FileName

                    c_lineStyle = 0
                    c_markerStyle = 0

                    for index, item in fileNames.iteritems():
                        if item in self.full_data_src_dict:
                            df = self.full_data_src_dict[item]
                            ax1.plot(df['dm'], df[column], label=self.full_data_label_dict[item], color=self.__lineColor[c_lineColor], linewidth=self.__lineWidth[c_lineWidth],
                                     linestyle=self.__lineStyle[c_lineStyle], alpha=alphaMainLine, marker=self.__markerStyle[c_markerStyle], markersize=self.__markerSize)

                        c_lineStyle += 1
                        c_markerStyle += 1

                    c_lineColor += 1
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
                ax1.set_xlabel(self.__xLabelPlotCrossSections, fontsize=self.__labelFont_size)
                ax1.set_ylabel(yLabel, fontsize=self.__labelFont_size)
                ax1.set_xlim(self.__mobilityDiamLimits[0], self.__mobilityDiamLimits[1])

                ax1.set_title(T1, fontsize=self.__plotTitle_Size)
                ax1.xaxis.set_major_formatter(FormatStrFormatter("%i"))
                ax1.xaxis.set_minor_formatter(FormatStrFormatter("%i"))
                ax1.legend(bbox_to_anchor=(1.00, 0.5), loc='center left', fontsize='large')
                Address = GF.getAddressTo(FolderName=self.__folderNameGraph + "\ErrorAndRatio", FileName=F1, Extension="jpg")
                plt.savefig(Address, format='jpg', dpi=self.__figureDPI, bbox_inches='tight')
                plt.clf()
                plt.close()

        except Exception as e:
            logging.exception(e)
            raise

    def PlotErrorAndRatioIntermediate(self, clName, tlName, yLabel):
        try:
            if self.__rhoEff100nmCTE == True:
                self.PlotErrorAndRatioCore(A1=self.rhoEff100nmSeries, A1Name=self.MainDF.AGG_EFF_RHO_100NM_CENTER, A2=self.DmSeries, A2Name=self.MainDF.AGG_EFF_DM_CENTER,
                                           column=clName, yLabel=yLabel, title=tlName, T1A=" for " + "$\\rho_{eff,100}$= ", T1I=1, T1B=" kg/m$^3$", F1A=" for " + "rho_eff=", F1I=0)

            if self.__DmCTE == True:
                self.PlotErrorAndRatioCore(A1=self.DmSeries, A1Name=self.MainDF.AGG_EFF_DM_CENTER, A2=self.sigmaMobSeries, A2Name=self.MainDF.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER,
                                           column=clName, yLabel=yLabel, title=tlName, T1A=" for " + "$D_m=$", T1I=2, F1A=" for " + "Dm=", F1I=2)

            if self.__sigmaMobCTE == True:
                self.PlotErrorAndRatioCore(A1=self.sigmaMobSeries, A1Name=self.MainDF.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER, A2=self.rhoEff100nmSeries,
                                           A2Name=self.MainDF.AGG_EFF_RHO_100NM_CENTER,
                                           column=clName, yLabel=yLabel, title=tlName, T1A=" for " + "$\sigma_p|d_m=$", T1I=2, F1A=" for " + "sigma_p=", F1I=2)

        except Exception as e:
            logging.exception(e)
            raise

    def PlotErrorAndRatio(self):
        try:
            ######################## Absorption Percentage Error
            clName = 'RealPercentErrorABS'
            tlName = "Percentage Error for TMatrix and RDG-FA Absorption Cross Section"
            self.PlotErrorAndRatioIntermediate(clName=clName, tlName=tlName, yLabel=self.__yLabelErrorAndRatio1)
            ######################## Absorption Absolute Percentage Error
            # clName = 'AbsolutePercentErrorABS'
            # tlName = "Absolute Percentage Error in TMatrix and RDG-FA Absorption Cross Section"
            # self.PlotErrorAndRatioIntermediate(clName=clName, tlName=tlName, yLabel=self.__yLabelErrorAndRatio1)
            ######################## Absorption Ratio
            clName = 'RatioABS'
            tlName = "TMatrix and RDG-FA Absorption Cross Section Ratio"
            self.PlotErrorAndRatioIntermediate(clName=clName, tlName=tlName, yLabel=self.__yLabelErrorAndRatio2)
            ######################## Scattering Percentage Error
            clName = 'RealPercentErrorSCA'
            tlName = "Percentage Error for TMatrix and RDG-FA Scattering Cross Section"
            self.PlotErrorAndRatioIntermediate(clName=clName, tlName=tlName, yLabel=self.__yLabelErrorAndRatio1)
            ######################## Scattering Absolute Percentage Error
            # clName = 'AbsolutePercentErrorSCA'
            # tlName = "Absolute Percentage Error in TMatrix and RDG-FA Scattering Cross Section"
            # self.PlotErrorAndRatioIntermediate(clName=clName, tlName=tlName, yLabel=self.__yLabelErrorAndRatio1)
            ######################## Scattering Ratio
            clName = 'RatioSCA'
            tlName = "TMatrix and RDG-FA Scattering Cross Section Ratio"
            self.PlotErrorAndRatioIntermediate(clName=clName, tlName=tlName, yLabel=self.__yLabelErrorAndRatio2)


        except Exception as e:
            logging.exception(e)
            raise

    ###################################
    ###################################
    ###################################
    ###################################
    def PlotCrossSectionsCore(self, A1, A1Name, A2, A2Name, T1A, T1I, F1A, F1I, title, column, T1B=None, F1B=None):
        try:
            alphaMainLine = self.__alphaMainLine

            for i in A1.unique():

                fig, ax1 = plt.subplots()

                c_lineColor = 0
                c_lineWidth = 0

                for j in A2.unique():

                    fileNames = self.MainDF[(A1Name == i) & (A2Name == j)].AA_FileName

                    c_lineStyle = 0
                    c_markerStyle = 0

                    for index, item in fileNames.iteritems():
                        if item in self.full_data_src_dict:
                            df = self.full_data_src_dict[item]
                            ax1.plot(df['dm'], df[column], label=self.full_data_label_dict[item], color=self.__lineColor[c_lineColor], linewidth=self.__lineWidth[c_lineWidth],
                                     linestyle=self.__lineStyle[c_lineStyle], alpha=alphaMainLine, marker=self.__markerStyle[c_markerStyle], markersize=self.__markerSize)

                        c_lineStyle += 1
                        c_markerStyle += 1

                    c_lineColor += 1
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
                ax1.set_xlabel(self.__xLabelPlotCrossSections, fontsize=self.__labelFont_size)
                ax1.set_ylabel(self.__yLabelPlotCrossSections, fontsize=self.__labelFont_size)
                ax1.set_xlim(self.__mobilityDiamLimits[0], self.__mobilityDiamLimits[1])

                ax1.set_title(T1, fontsize=self.__plotTitle_Size)
                ax1.xaxis.set_major_formatter(FormatStrFormatter("%i"))
                ax1.xaxis.set_minor_formatter(FormatStrFormatter("%i"))
                ax1.legend(bbox_to_anchor=(1.00, 0.5), loc='center left', fontsize='large')
                Address = GF.getAddressTo(FolderName=self.__folderNameGraph + "\CrossSections", FileName=F1, Extension="jpg")
                plt.savefig(Address, format='jpg', dpi=self.__figureDPI, bbox_inches='tight')
                plt.clf()
                plt.close()

        except Exception as e:
            logging.exception(e)
            raise

    def PlotCrossSectionsIntermediate(self, clName, tlName):
        try:
            if self.__rhoEff100nmCTE == True:
                self.PlotCrossSectionsCore(A1=self.rhoEff100nmSeries, A1Name=self.MainDF.AGG_EFF_RHO_100NM_CENTER, A2=self.DmSeries, A2Name=self.MainDF.AGG_EFF_DM_CENTER,
                                           column=clName, title=tlName, T1A=" for " + "$\\rho_{eff,100}$= ", T1I=1, T1B=" kg/m$^3$", F1A=" for " + "rho_eff=", F1I=0)

            if self.__DmCTE == True:
                self.PlotCrossSectionsCore(A1=self.DmSeries, A1Name=self.MainDF.AGG_EFF_DM_CENTER, A2=self.sigmaMobSeries, A2Name=self.MainDF.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER,
                                           column=clName, title=tlName, T1A=" for " + "$D_m=$", T1I=2, F1A=" for " + "Dm=", F1I=2)

            if self.__sigmaMobCTE == True:
                self.PlotCrossSectionsCore(A1=self.sigmaMobSeries, A1Name=self.MainDF.AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER, A2=self.rhoEff100nmSeries,
                                           A2Name=self.MainDF.AGG_EFF_RHO_100NM_CENTER,
                                           column=clName, title=tlName, T1A=" for " + "$\sigma_p|d_m=$", T1I=2, F1A=" for " + "sigma_p=", F1I=2)

        except Exception as e:
            logging.exception(e)
            raise

    def PlotCrossSections(self):
        try:
            ######################## Absorption RDG
            clName = 'ABS_RDG'
            tlName = "RDG-FA Absorption Cross Section"
            self.PlotCrossSectionsIntermediate(clName=clName, tlName=tlName)
            ######################## Scattering RDG
            clName = 'SCA_RDG'
            tlName = "RDG-FA Scattering Cross Section"
            self.PlotCrossSectionsIntermediate(clName=clName, tlName=tlName)
            ######################## Absorption TMatrix
            clName = 'ABS_TMatrix'
            tlName = "TMatrix Absorption Cross Section"
            self.PlotCrossSectionsIntermediate(clName=clName, tlName=tlName)
            ######################## Scattering TMatrix
            clName = 'SCA_TMatrix'
            tlName = "TMatrix Scattering Cross Section"
            self.PlotCrossSectionsIntermediate(clName=clName, tlName=tlName)

        except Exception as e:
            logging.exception(e)
            raise

    ###################################
    ###################################
    ###################################
    ###################################

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
