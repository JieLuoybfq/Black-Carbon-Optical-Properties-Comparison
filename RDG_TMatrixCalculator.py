# Keyhan Babaee, https://github.com/KeyhanB
# V1.2
# Oct 2019
import os

from matplotlib.ticker import FormatStrFormatter
from sklearn.metrics import explained_variance_score, r2_score

from ConfigReaderModule import logging
from math import *
import numpy as np
import AppFunctions as FN
from BoundaryFinder import BCDBBoundaryCheck as BCheck
from decimal import Decimal
import pandas as pd
import GeneralFunctions as GF
import gmpy2
from scipy.ndimage import gaussian_filter1d
import CalculationCore
import TotalCrossSectionCalculator
from matplotlib import pyplot as plt
from matplotlib import rcParams

####### Plotting Parameters
rcParams['mathtext.fontset'] = 'stix'
rcParams['font.family'] = 'STIXGeneral'


class KeyhanV2:
    def __init__(self, inputDict):
        try:
            ######################################################################################
            time_now = GF.GetDateAndTimeUTCNow()
            temp = {}
            temp['AA_FileName'] = f"{time_now}.csv"
            temp['AA_Plot'] = 1
            self.infoDict = {**temp, **inputDict}
            ######################################################################################
            self.__AGG_EFF_DM_CENTER = inputDict['AGG_EFF_DM_CENTER']
            self.__AGG_EFF_RHO_100NM_CENTER = inputDict['AGG_EFF_RHO_100NM_CENTER']
            self.__AGG_MATERIAL_DENSITY_CENTER = inputDict['AGG_MATERIAL_DENSITY_CENTER']
            self.__AGG_MOBILITY_DIAMETER_CENTER_MAX = inputDict['AGG_MOBILITY_DIAMETER_CENTER_MAX']
            self.__AGG_MOBILITY_DIAMETER_CENTER_MIN = inputDict['AGG_MOBILITY_DIAMETER_CENTER_MIN']
            self.__AGG_PREFACTOR_PROJECTED_AREA_COEFFICIENT_CENTER = inputDict['AGG_PREFACTOR_PROJECTED_AREA_COEFFICIENT_CENTER']
            self.__AGG_EXPONENT_PROJECTED_AREA_COEFFICIENT_CENTER = inputDict['AGG_EXPONENT_PROJECTED_AREA_COEFFICIENT_CENTER']
            self.__AGG_FRACTAL_DIMENSION_CENTER = inputDict['AGG_FRACTAL_DIMENSION_CENTER']
            self.__AGG_FRACTAL_PREFACTOR_CENTER = inputDict['AGG_FRACTAL_PREFACTOR_CENTER']
            self.__AGG_RI_REAL_CENTER = inputDict['AGG_RI_REAL_CENTER']
            self.__AGG_RI_IMAG_CENTER = inputDict['AGG_RI_IMAG_CENTER']
            self.__AGG_WLENGTH_CENTER = inputDict['AGG_WLENGTH_CENTER']
            self.__AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER = inputDict['AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER']
            ######################################################################################
            self.__Sample_Sigma_Bins = 23  # Number of bins
            self.__Primary_Sigma_dm_CTE_Bound = 3  # Number of Sigma G to cover
            self.__Primary_Sigma_dm_CTE_Nt = 20  # Number of bins
            ######################################################################################
            self.infoDict['MobilityBins'] = self.__Sample_Sigma_Bins
            self.infoDict['NumberOfSigma'] = self.__Primary_Sigma_dm_CTE_Bound
            self.infoDict['PrimaryBins'] = self.__Primary_Sigma_dm_CTE_Nt
            self.infoDict['FittedSigma'] = 0
            self.infoDict['FittedMedian'] = 0
            ######################################################################################
            self.__PlotDetails = True
            self.__dmSelectForPlot = 4  # Number of bins
            self.__folderNameGraph = 'Detail Graphs'
            self.__OnlySVGSaving = False
            self.__figureDPI = 400
            ######################################################################################
            self.__dmCalc = True
            self.__dmSigma = 1.4
            self.__dmMedian = 130
            self.__dpMin = 1
            self.__dpMax = 105
            self.__dpBin = 50
            ######################################################################################
            self.__GoBeyondBound = True
            self.__TMatrixActive = False
            self.__RDGActive = True
            self.__TraditionalDistribution = True
            ######################################################################################
        except Exception as e:
            logging.exception(e)
            raise

    def KeyhanV2Calc(self, DB_Info):
        try:
            self.arrMobilityDiamNano = self.CalcMobilityDiamBins()
            self.dict_dpMedianNano = self.CalcPrimaryParticleSizeMedianNano()
            self.dict_dpMobDistribNano, self.dict_dpMobDistribChance = self.CalcPrimDiamBinAtEachMob()
            self.dict_NpMobDistrib = self.CalcPrimaryParticleNumber()

            dictSuggested = self.CreateDictForEvaluation()

            dictChecked = self.CheckDictWithBoundary(dict=dictSuggested, Beyond=self.__GoBeyondBound)
            dictConvertedTR = self.ConvertDictToArray(dict=dictChecked)

            dict_dpMedianReal, dictCheckedRealNpdp = self.CalcRealNpdpDistribution(dictChecked=dictChecked, dictSuggested=dictSuggested)
            '''
            dictConvertedRE = self.ConvertDictToArray(dict=dictCheckedRealNpdp)

            self.RDGTMCore(DB_Info=DB_Info, dictConverted=dictConvertedTR, dictChecked=dictChecked, mediandpDict=self.dict_dpMedianNano, fileAppend="TR", RDGActive=self.__RDGActive,
                           TmatrixActive=self.__TMatrixActive)

            if self.__TraditionalDistribution:
                self.RDGTMCore(DB_Info=DB_Info, dictConverted=dictConvertedRE, dictChecked=dictCheckedRealNpdp, mediandpDict=dict_dpMedianReal, fileAppend="RE", RDGActive=self.__RDGActive,
                               TmatrixActive=self.__TMatrixActive)
            '''
            if os.path.exists('TMatrix_RDG_Result\Beacon.csv'):
                dfInfoDB = pd.read_csv(f"TMatrix_RDG_Result\Beacon.csv")
                dfInfoDB.loc[len(dfInfoDB)] = self.infoDict
                dfInfoDB.to_csv(f"TMatrix_RDG_Result\Beacon.csv", index=False)
            else:
                newInfoDf = pd.DataFrame([self.infoDict])
                newInfoDf.to_csv(f"TMatrix_RDG_Result\Beacon.csv", index=False)
                

        except Exception as e:
            logging.exception(e)
            raise

    ######################################################################################
    ######################################################################################
    ######################################################################################
    ######################################################################################
    ######################################################################################
    ######################################################################################
    ######################################################################################

    def RDGTMCore(self, DB_Info, dictConverted, dictChecked, mediandpDict, fileAppend, RDGActive, TmatrixActive):
        try:
            CC = CalculationCore.CalculationCore(infoDict=DB_Info, calcDict=dictConverted, RDGActive=RDGActive, TmatrixActive=TmatrixActive)
            dictOutputT1, dictOutputR1 = CC.Calc()

            TCSC = TotalCrossSectionCalculator.TotalCrossSectionCalculator(infoDict=self.infoDict, outTMatrix=dictOutputT1, outRDG=dictOutputR1, checkedDict=dictChecked,
                                                                           dict_dpMedianNano=mediandpDict, RDGActive=RDGActive, TmatrixActive=TmatrixActive)
            dfResult = TCSC.Calc()
            dfResult.to_csv(f"TMatrix_RDG_Result\{fileAppend}_{self.infoDict['AA_FileName']}", index=False)

        except Exception as e:
            logging.exception(e)
            raise

    def PlotdpDistribution(self, dict, resR, resC, dmFalse_resR, dmFalse_resC, dictMedian, title):
        try:

            yLimits = [1, 100]
            xLimits = [49, 1500]
            lineSt = ['-', '--', '-.']
            ######## creating dp array
            bound_D_Max = yLimits[1]
            bound_D_Min = yLimits[0]
            total_Number_Bins = 49
            D_Ratio = (bound_D_Max / bound_D_Min) ** (1 / (total_Number_Bins - 1))
            dp_Nano_tot = []
            for i in range(0, total_Number_Bins):
                d1 = bound_D_Min * (D_Ratio ** i)
                dp_Nano_tot.append(round(d1, 3))

            ############ getting median diameter for the dict
            dpMedian = []
            for dm in self.arrMobilityDiamNano:
                dpMedian.append(dictMedian[dm])

            fig, ax1 = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(12, 6))
            ##########################################
            ##########################################
            ########################################## chance and Np

            ax1[0].plot(self.arrMobilityDiamNano, dpMedian, color='black', label="Median diameter")
            ax1[0].set_xscale('log')
            ##########################################
            factor = 1.8
            dm_Base = self.dmRandom[0]
            counter = 0
            for dm in self.dmRandom:
                l1 = len(dict[dm]['Df'])
                if l1 == 1:
                    ##### one point of random
                    ax1[0].plot([dm], [dict[dm]['dp'][0]], label=f'd$_p$ for dm:{round(dm)} nm', color='red', marker='*', markersize=9)
                else:
                    dp = dict[dm]['dp']
                    dpChance = dict[dm]['chance']
                    Np = dict[dm]['Np']
                    x_adj = []
                    x_base = []
                    for i in range(len(dpChance)):
                        x_adj.append(dm - (float(dpChance[i]) * float(Np[i]) * factor ** (dm_Base / dm)))
                        x_base.append(dm)
                    ax1[0].plot(x_adj, dp, label=f'd$_p$ distribution for dm:{round(dm)} nm', color='green', linestyle=lineSt[counter])
                    counter += 1
                    ax1[0].plot(x_base, dp, color='gray', alpha=0.75)

            #####################  total distribution: dm is calculated
            logNormal_with_dm = []
            sum = 0
            if resR['Sigma_G'] != 1:
                for i in range(0, total_Number_Bins):
                    if i != 0:
                        l1 = self._calcLogNDistribPDF(resR['D_Median'], resR['Sigma_G'], dp_Nano_tot[i], dp_Nano_tot[i - 1])
                    else:
                        l1 = 0
                    sum += l1
                    logNormal_with_dm.append(l1)
            else:
                dp_Nano_tot_with_dm = []
                dp_Nano_tot_with_dm.append(resR['D_Median'])
                logNormal_with_dm.append(0.5)
            #####################  total distribution: dm is not calculated
            logNormal_without_dm = []
            sumNO = 0
            if dmFalse_resR['Sigma_G'] != 1:
                for i in range(0, total_Number_Bins):
                    if i != 0:
                        l1 = self._calcLogNDistribPDF(dmFalse_resR['D_Median'], dmFalse_resR['Sigma_G'], dp_Nano_tot[i], dp_Nano_tot[i - 1])
                    else:
                        l1 = 0
                    sumNO += l1
                    logNormal_without_dm.append(l1)
            else:
                dp_Nano_tot_without_dm = []
                dp_Nano_tot_without_dm.append(dmFalse_resR['D_Median'])
                logNormal_without_dm.append(1)
            ##################################### Main dm distribution
            logNormal_Main_dm = []
            sumMain = 0
            if self.__dmSigma != 1:
                for i in range(0, self.__Sample_Sigma_Bins):
                    if i != 0:
                        l1 = 400 * self._calcLogNDistribPDF(self.__dmMedian, self.__dmSigma, self.arrMobilityDiamNano[i], self.arrMobilityDiamNano[i - 1])
                    else:
                        l1 = 0
                    sumMain += l1
                    logNormal_Main_dm.append(l1)
            else:
                dp_Nano_Main_dm = []
                dp_Nano_Main_dm.append(self.__dmMedian)
                logNormal_Main_dm.append(55)
            #####################################

            last_dm = self.arrMobilityDiamNano[-1] * 1.5
            x_adj = []
            x_base = []
            if resR['Sigma_G'] != 1:
                for i in range(len(logNormal_with_dm)):
                    x_adj.append(last_dm - (float(logNormal_with_dm[i]) * (last_dm / dm_Base) * 250))
                    x_base.append(last_dm)
                ax1[0].plot(x_adj, dp_Nano_tot, label=f'Total distribution averaged over d$_m$', color='black')
                ax1[0].plot(x_base, dp_Nano_tot, color='gray')
            else:
                x_adj = last_dm
                ax1[0].plot(x_adj, dp_Nano_tot_with_dm, label=f'Total distribution averaged over d$_m$', color='black', marker='*', markersize=9)

            #####################################
            x_adj = []
            x_base = []
            if dmFalse_resR['Sigma_G'] != 1:
                for i in range(len(logNormal_without_dm)):
                    x_adj.append(last_dm - (float(logNormal_without_dm[i]) * (last_dm / dm_Base) * 250))
                    x_base.append(last_dm)
                ax1[0].plot(x_adj, dp_Nano_tot, label=f'Total distribution with constant d$_m$', color='gray', linestyle=':', alpha=0.85)
                ax1[0].plot(x_base, dp_Nano_tot, color='gray')
            else:
                x_adj = last_dm
                ax1[0].plot(x_adj, dp_Nano_tot_without_dm, label=f'Total distribution with constant d$_m$', color='gray', marker='s', markersize=11, alpha=0.85)
            #####################################

            if self.__dmSigma != 1:
                ax1[0].plot(self.arrMobilityDiamNano, logNormal_Main_dm, label=f'd$_m$ distribution, median:{self.__dmMedian}, sigma:{self.__dmSigma}', color='rosybrown', linestyle='--', alpha=0.8)
            else:
                ax1[0].plot(dp_Nano_Main_dm, logNormal_Main_dm, label=f'd$_m$ distribution, median:{self.__dmMedian}, sigma:{self.__dmSigma}', color='rosybrown', marker='*', markersize=9)

            #####################################
            ax1[0].xaxis.set_major_formatter(FormatStrFormatter("%i"))
            # ax1[0].xaxis.set_major_locator(plt.MaxNLocator(6))
            ax1[0].xaxis.set_tick_params(labelsize=12)
            ax1[0].yaxis.set_tick_params(labelsize=12)
            ax1[0].set_ylim(yLimits)
            ax1[0].set_xlim(xLimits)
            ax1[0].set_title(f"Considering [Np*Likelihood], Median:{round(resR['D_Median'], 1)}, Sigma:{round(resR['Sigma_G'], 2)}")
            # ax1[0].legend()

            ####################################################################################
            ####################################################################################
            ####################################################################################
            ####################################################################################
            ####################################################################################
            ####################################################################################
            ##########################################
            ##########################################
            ##########################################
            ##########################################
            ##########################################
            ##########################################
            ########## only chance
            ax1[1].plot(self.arrMobilityDiamNano, dpMedian, color='black', label="Median Diameter")
            ax1[1].set_xscale('log')
            ##########################################
            factor = 2.1
            dm_Base = self.dmRandom[0]
            counter = 0
            for dm in self.dmRandom:
                l1 = len(dict[dm]['Df'])
                if l1 == 1:
                    ##### one point of random
                    ax1[1].plot([dm], [dict[dm]['dp'][0]], label=f'd$_p$ for dm:{round(dm)} nm', color='red', marker='*', markersize=9)
                else:
                    dp = dict[dm]['dp']
                    dpChance = dict[dm]['chance']
                    # Np = dict[dm]['Np']
                    x_adj = []
                    x_base = []
                    for i in range(len(dpChance)):
                        x_adj.append(dm - (float(dpChance[i]) * 80 * factor ** (dm / dm_Base)))
                        x_base.append(dm)
                    ax1[1].plot(x_adj, dp, label=f'd$_p$ distribution for dm:{round(dm)} nm', color='green', linestyle=lineSt[counter])
                    counter += 1
                    ax1[1].plot(x_base, dp, color='gray', alpha=0.75)

            #####################  total distribution: dm is calculated
            logNormal_with_dm = []
            sum = 0
            if resC['Sigma_G'] != 1:
                for i in range(0, total_Number_Bins):
                    if i != 0:
                        l1 = self._calcLogNDistribPDF(resC['D_Median'], resC['Sigma_G'], dp_Nano_tot[i], dp_Nano_tot[i - 1])
                    else:
                        l1 = 0
                    sum += l1
                    logNormal_with_dm.append(l1)
            else:
                dp_Nano_tot_with_dm = []
                dp_Nano_tot_with_dm.append(resC['D_Median'])
                logNormal_with_dm.append(0.5)
            #####################  total distribution: dm is not calculated
            logNormal_without_dm = []
            sumNO = 0
            if dmFalse_resC['Sigma_G'] != 1:
                for i in range(0, total_Number_Bins):
                    if i != 0:
                        l1 = self._calcLogNDistribPDF(dmFalse_resC['D_Median'], dmFalse_resC['Sigma_G'], dp_Nano_tot[i], dp_Nano_tot[i - 1])
                    else:
                        l1 = 0
                    sumNO += l1
                    logNormal_without_dm.append(l1)
            else:
                dp_Nano_tot_without_dm = []
                dp_Nano_tot_without_dm.append(dmFalse_resC['D_Median'])
                logNormal_without_dm.append(1)
            ##################################### Main dm distribution
            logNormal_Main_dm = []
            sumMain = 0
            if self.__dmSigma != 1:
                for i in range(0, self.__Sample_Sigma_Bins):
                    if i != 0:
                        l1 = 400 * self._calcLogNDistribPDF(self.__dmMedian, self.__dmSigma, self.arrMobilityDiamNano[i], self.arrMobilityDiamNano[i - 1])
                    else:
                        l1 = 0
                    sumMain += l1
                    logNormal_Main_dm.append(l1)
            else:
                dp_Nano_Main_dm = []
                dp_Nano_Main_dm.append(self.__dmMedian)
                logNormal_Main_dm.append(55)
            #####################################

            last_dm = self.arrMobilityDiamNano[-1] * 1.5
            x_adj = []
            x_base = []
            if resC['Sigma_G'] != 1:
                for i in range(len(logNormal_with_dm)):
                    x_adj.append(last_dm - (float(logNormal_with_dm[i]) * (last_dm / dm_Base) * 250))
                    x_base.append(last_dm)
                ax1[1].plot(x_adj, dp_Nano_tot, label=f'Total distribution averaged over d$_m$', color='black')
                ax1[1].plot(x_base, dp_Nano_tot, color='gray')
            else:
                x_adj = last_dm
                ax1[1].plot(x_adj, dp_Nano_tot_with_dm, label=f'Total distribution averaged over d$_m$', color='black', marker='*', markersize=9)

            #####################################
            x_adj = []
            x_base = []
            if dmFalse_resC['Sigma_G'] != 1:
                for i in range(len(logNormal_without_dm)):
                    x_adj.append(last_dm - (float(logNormal_without_dm[i]) * (last_dm / dm_Base) * 250))
                    x_base.append(last_dm)
                ax1[1].plot(x_adj, dp_Nano_tot, label=f'Total distribution with constant d$_m$', color='gray', linestyle=':', alpha=0.85)
                ax1[1].plot(x_base, dp_Nano_tot, color='gray')
            else:
                x_adj = last_dm
                ax1[1].plot(x_adj, dp_Nano_tot_without_dm, label=f'Total distribution with constant d$_m$', color='gray', marker='s', markersize=11, alpha=0.85)
            #####################################

            if self.__dmSigma != 1:
                ax1[1].plot(self.arrMobilityDiamNano, logNormal_Main_dm, label=f'd$_m$ distribution, median:{self.__dmMedian}, sigma:{self.__dmSigma}', color='rosybrown', linestyle='--', alpha=0.8)
            else:
                ax1[1].plot(dp_Nano_Main_dm, logNormal_Main_dm, label=f'd$_m$ distribution, median:{self.__dmMedian}, sigma:{self.__dmSigma}', color='rosybrown', marker='*', markersize=9)

            ####################################
            ax1[1].xaxis.set_major_formatter(FormatStrFormatter("%i"))
            # ax1[1].xaxis.set_major_locator(plt.MaxNLocator(6))
            ax1[1].xaxis.set_tick_params(labelsize=12)
            ax1[1].yaxis.set_tick_params(labelsize=12)
            # ax1[1].xaxis.set_major_formatter(FormatStrFormatter('%1.2f'))
            # ax1[1].yaxis.set_major_formatter(FormatStrFormatter('%1.2f'))
            ax1[1].set_ylim(yLimits)
            ax1[1].set_xlim(xLimits)
            ax1[1].set_title(f"Considering Likelihood Only, Median:{round(resC['D_Median'], 1)}, Sigma:{round(resC['Sigma_G'], 2)}")
            # ax1[1].legend()
            ########################
            ########################
            ########################
            fig.suptitle(f"{title}_Dm:{self.__AGG_EFF_DM_CENTER}_rhoeff:{self.__AGG_EFF_RHO_100NM_CENTER}_Sigma_p:{self.__AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER}")
            fig.subplots_adjust(top=0.88, wspace=0.03, hspace=0.03)

            fig.text(0.5, 0.023, 'Mobility Diameter (nm)', ha='center', fontsize=14)
            fig.text(0.08, 0.5, 'Primary Particle Diameter (nm)', va='center', rotation='vertical', fontsize=14)
            leg = plt.legend(bbox_to_anchor=(1.05, 0.5), markerscale=1, fontsize=12, loc='lower left', borderaxespad=0.)
            ################# set the lineWidth of each legend object
            for legObj in leg.legendHandles:
                legObj.set_linewidth(2)

            self.SaveAndClosePlot(
                folderName="Np_dp_Distribution",
                F1=f"{title}_Dm_{self.__AGG_EFF_DM_CENTER}_rhoeff_{self.__AGG_EFF_RHO_100NM_CENTER}_Sigma_p_{self.__AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER}")

        except Exception as e:
            logging.exception(e)
            raise

    def PlotdpDistributionTot(self, dict, resR, dictMedian, title, resFK):
        try:

            yLimits = [1, 60]
            xLimits = [49, 1500]
            lineSt = ['-', '--', '-.']
            ######## creating dp array
            bound_D_Max = yLimits[1]
            bound_D_Min = yLimits[0]
            total_Number_Bins = 49
            D_Ratio = (bound_D_Max / bound_D_Min) ** (1 / (total_Number_Bins - 1))
            dp_Nano_tot = []
            for i in range(0, total_Number_Bins):
                d1 = bound_D_Min * (D_Ratio ** i)
                dp_Nano_tot.append(round(d1, 3))

            ############ getting median diameter for the dict
            dpMedian = []
            for dm in self.arrMobilityDiamNano:
                dpMedian.append(dictMedian[dm])

            fig, ax1 = plt.subplots()
            ##########################################
            ##########################################
            ########################################## chance and Np

            ax1.plot(self.arrMobilityDiamNano, dpMedian, color='black', label="Median diameter")
            ax1.set_xscale('log')
            ##########################################
            factor = 1.6
            dm_Base = self.dmRandom[0]
            counter = 0
            for dm in self.dmRandom:
                l1 = len(dict[dm]['Df'])
                if l1 == 1:
                    ##### one point of random
                    ax1.plot([dm], [dict[dm]['dp'][0]], label=f'd$_p$ for dm:{round(dm)} nm', color='red', marker='*', markersize=9)
                else:
                    dp = dict[dm]['dp']
                    dpChance = dict[dm]['chance']
                    Np = dict[dm]['Np']
                    x_adj = []
                    x_base = []
                    for i in range(len(dpChance)):
                        x_adj.append(dm - (float(dpChance[i]) * float(Np[i]) * factor ** (dm_Base / dm)))
                        x_base.append(dm)
                    ax1.plot(x_adj, dp, label=f'd$_p$ distribution for dm:{round(dm)} nm', color='green', linestyle=lineSt[counter])
                    counter += 1
                    ax1.plot(x_base, dp, color='gray', alpha=0.75)

            #####################
            logNormal_with_dm = []
            sum = 0
            if resR['Sigma_G'] != 1:
                for i in range(0, total_Number_Bins):
                    if i != 0:
                        l1 = self._calcLogNDistribPDF(resR['D_Median'], resR['Sigma_G'], dp_Nano_tot[i], dp_Nano_tot[i - 1])
                    else:
                        l1 = 0
                    sum += l1
                    logNormal_with_dm.append(l1)
            else:
                dp_Nano_tot_with_dm = []
                dp_Nano_tot_with_dm.append(resR['D_Median'])
                logNormal_with_dm.append(0.5)

            #####################################

            last_dm = self.arrMobilityDiamNano[-1] * 1.5
            x_adj = []
            x_base = []
            if resR['Sigma_G'] != 1:
                for i in range(len(logNormal_with_dm)):
                    x_adj.append(last_dm - (float(logNormal_with_dm[i]) * (last_dm / dm_Base) * 250))
                    x_base.append(last_dm)
                ax1.plot(x_adj, dp_Nano_tot, label=f'Total distribution averaged over d$_m$', color='black')
                ax1.plot(x_base, dp_Nano_tot, color='gray')
            else:
                x_adj = last_dm
                ax1.plot(x_adj, dp_Nano_tot_with_dm, label=f'Total distribution averaged over d$_m$', color='black', marker='*', markersize=9)

            #####################################

            #####################################
            ax1.xaxis.set_major_formatter(FormatStrFormatter("%i"))
            # ax1[0].xaxis.set_major_locator(plt.MaxNLocator(6))
            ax1.xaxis.set_tick_params(labelsize=12)
            ax1.yaxis.set_tick_params(labelsize=12)
            ax1.set_ylim(yLimits)
            ax1.set_xlim(xLimits)
            ax1.set_title(f"Considering [Np*Likelihood], Median:{round(resFK['D_Median'], 1)}, Sigma:{round(resFK['Sigma_G'], 2)}")
            # ax1[0].legend()

            ########################
            fig.suptitle(f"{title}_Dm:{self.__AGG_EFF_DM_CENTER}_rhoeff:{self.__AGG_EFF_RHO_100NM_CENTER}_Sigma_p:{self.__AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER}")
            fig.subplots_adjust(top=0.88, wspace=0.03, hspace=0.03)

            fig.text(0.5, 0.023, 'Mobility Diameter (nm)', ha='center', fontsize=14)
            fig.text(0.04, 0.5, 'Primary Particle Diameter (nm)', va='center', rotation='vertical', fontsize=14)
            leg = plt.legend(bbox_to_anchor=(1.05, 0.5), markerscale=1, fontsize=12, loc='lower left', borderaxespad=0.)
            ################# set the lineWidth of each legend object
            for legObj in leg.legendHandles:
                legObj.set_linewidth(2)

            self.SaveAndClosePlot(
                folderName="Np_dp_Distribution",
                F1=f"{title}_Dm_{self.__AGG_EFF_DM_CENTER}_rhoeff_{self.__AGG_EFF_RHO_100NM_CENTER}_Sigma_p_{self.__AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER}")

        except Exception as e:
            logging.exception(e)
            raise

    def LogNormalFit(self, x, y):
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
                Sigma_G = round(Sigma_G, 3)
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
                ######### check goodness of fitting

                logNormalFitted = []
                sum = 0
                if Sigma_G != 1:
                    for i in range(len(x)):
                        if i != 0:
                            l1 = self._calcLogNDistribPDF(D_Median, Sigma_G, x[i], x[i - 1])
                        else:
                            l1 = 0
                        sum += l1
                        logNormalFitted.append(l1)

                else:
                    logNormalFitted.append(1)

                ratio = max(y) / max(logNormalFitted)

                fittedCorrected = []
                for i in logNormalFitted:
                    fittedCorrected.append(i * ratio)

                # plt.plot(x,y)
                # plt.plot(x, fittedCorrected)
                # plt.show()
                if Sigma_G != 1:
                    explainedVarianceScore = explained_variance_score(y, fittedCorrected)
                    R2score = r2_score(y, fittedCorrected)
                else:
                    explainedVarianceScore = 1
                    R2score = 1

                result = {'D_G': round(D_G, 1), 'Sigma_G': round(Sigma_G, 2), 'Total_Conc': round(Total_Conc, 1), 'D_Median': round(D_Median, 1), 'R2score': round(R2score, 2),
                          'explainedVarianceScore': round(explainedVarianceScore, 2)}
                return result

            else:
                logging.exception(f"Lognormal fit error: array lengths are not equal {len(x)},{len(y)}")
            ##################

        except Exception as e:
            logging.exception(e)
            raise

    def _CalcRealNpDp(self, dict, dmCalc, Np_multiple):
        try:

            if dmCalc:
                logNormalPDFdm = []
                sum = 0
                if self.__dmSigma != 1:
                    for i in range(0, self.__Sample_Sigma_Bins):
                        if i != 0:
                            l1 = self._calcLogNDistribPDF(self.__dmMedian, self.__dmSigma, self.arrMobilityDiamNano[i], self.arrMobilityDiamNano[i - 1])
                        else:
                            l1 = 0
                        sum += l1
                        logNormalPDFdm.append(l1)

                else:
                    logNormalPDFdm.append(1)

            #################################################
            #################################################
            bound_D_MaxCC = self.__dpMax
            bound_D_MinCC = self.__dpMin
            total_Number_BinsCC = self.__dpBin
            #################################################
            #################################################
            dpArr = []
            NpChanceArr = []
            counter = 0
            for dm in dict:

                if dmCalc == True:
                    dmChance = logNormalPDFdm[counter] * 100
                    counter += 1
                else:
                    dmChance = 1

                # number of aggregates in the specific dm
                l1 = len(dict[dm]['Df'])

                if l1 == 1:
                    Np = dict[dm]['Np'][:1]
                    dp = dict[dm]['dp'][:1]
                    chance = dict[dm]['chance'][:1]
                    for i in range(len(dp)):
                        dpArr.append(dp[i])
                        if Np_multiple:
                            NpChanceArr.append(chance[i] * Np[i] * Decimal(dmChance))
                        else:
                            NpChanceArr.append(chance[i] * Decimal(1) * Decimal(dmChance))

                else:
                    Np = dict[dm]['Np']
                    dp = dict[dm]['dp']
                    chance = dict[dm]['chance']
                    for i in range(len(dp)):
                        dpArr.append(dp[i])
                        if Np_multiple:
                            NpChanceArr.append(chance[i] * Np[i] * Decimal(dmChance))
                        else:
                            NpChanceArr.append(chance[i] * Decimal(1) * Decimal(dmChance))

            diameter_Nano = []
            NpChanceBinArr = []
            D_Ratio = (bound_D_MaxCC / bound_D_MinCC) ** (1 / (total_Number_BinsCC - 1))

            for i in range(0, total_Number_BinsCC):
                d1 = bound_D_MinCC * (D_Ratio ** i)
                diameter_Nano.append(round(d1, 3))
                NpChanceBinArr.append(0)

            for i in range(1, len(diameter_Nano)):
                for p in range(len(dpArr)):
                    if dpArr[p] < diameter_Nano[i] and dpArr[p] > diameter_Nano[i - 1]:
                        NpChanceBinArr[i] = NpChanceBinArr[i] + float(NpChanceArr[p]) * 10000

            res = self.LogNormalFit(diameter_Nano, NpChanceBinArr)
            logging.info(
                f"LogNormal fit res:___D_G:{res['D_G']}___Sigma_G:{res['Sigma_G']}___Total_Conc:{res['Total_Conc']}___D_Median:{res['D_Median']}___R2:{res['R2score']}___explainedVarianceScore:{res['explainedVarianceScore']}")
            return res

        except Exception as e:
            logging.exception(e)
            raise

    def RealNpDp(self, dict):
        try:

            resR = self._CalcRealNpDp(dict=dict, dmCalc=self.__dmCalc, Np_multiple=True)
            resC = self._CalcRealNpDp(dict=dict, dmCalc=self.__dmCalc, Np_multiple=False)

            dmFalse_resR = self._CalcRealNpDp(dict=dict, dmCalc=False, Np_multiple=True)
            dmFalse_resC = self._CalcRealNpDp(dict=dict, dmCalc=False, Np_multiple=False)

            return resR, resC, dmFalse_resR, dmFalse_resC
        except Exception as e:
            logging.exception(e)
            raise

    def CalcRealNpdpDistribution(self, dictChecked, dictSuggested):
        try:
            # resSuggestedR, resSuggestedC = self.RealNpDp(dictSuggested)
            resCheckedR, resCheckedC, dmFalse_resR, dmFalse_resC = self.RealNpDp(dictChecked)

            if self.__PlotDetails:
                self.PlotdpDistribution(dictChecked, resCheckedR, resCheckedC, dmFalse_resR, dmFalse_resC, self.dict_dpMedianNano, "Np_dp")
                # self.PlotdpDistribution(dictSuggested, resSuggestedR, resSuggestedC, self.dict_dpMedianNano, "Obs-Suggested")

            ############################
            res = resCheckedR
            ############################
            self.infoDict['FittedSigma'] = res['Sigma_G']
            self.infoDict['FittedMedian'] = res['D_Median']
            self.infoDict['R2score'] = res['R2score']
            self.infoDict['explainedVarianceScore'] = res['explainedVarianceScore']
            ############################
            ############################
            ############################
            bound_D_Max = res['D_Median'] * (res['Sigma_G'] ** self.__Primary_Sigma_dm_CTE_Bound)
            bound_D_Min = res['D_Median'] * (res['Sigma_G'] ** (-1 * self.__Primary_Sigma_dm_CTE_Bound))
            total_Number_Bins = self.__Primary_Sigma_dm_CTE_Nt

            diameter_Nano = []
            logNormalPDF = []
            D_Ratio = (bound_D_Max / bound_D_Min) ** (1 / (total_Number_Bins - 1))

            for i in range(0, total_Number_Bins):
                d1 = bound_D_Min * (D_Ratio ** i)
                diameter_Nano.append(round(d1, 3))

            sum = 0
            if res['Sigma_G'] != 1:
                for i in range(0, total_Number_Bins):
                    if i != 0:
                        l1 = self._calcLogNDistribPDF(res['D_Median'], res['Sigma_G'], diameter_Nano[i], diameter_Nano[i - 1])
                    else:
                        l1 = 0
                    sum += l1
                    logNormalPDF.append(Decimal(l1))

            else:
                diameter_Nano = []
                diameter_Nano.append(res['D_Median'])
                logNormalPDF.append(Decimal(1))

            primaryDiam = diameter_Nano
            primaryChance = logNormalPDF
            ################################
            ################################
            dictNp = {}
            dictChp = {}
            for dm in self.arrMobilityDiamNano:
                arrNp = []
                arrNpChance = []
                dpArr = primaryDiam
                chanceArr = primaryChance
                for i in range(len(dpArr)):
                    np1 = self._calcPPN(dm=dm, dp=dpArr[i])
                    arrNp.append(np1)
                    arrNpChance.append(100000 * float(chanceArr[i] / np1))
                res1 = self.LogNormalFit(x=dpArr, y=arrNpChance)
                #################
                realChance = []
                sum = 0
                if res1['Sigma_G'] != 1:
                    for i in range(len(dpArr)):
                        if i != 0:
                            l1 = self._calcLogNDistribPDF(res1['D_Median'], res1['Sigma_G'], dpArr[i], dpArr[i - 1])
                        else:
                            l1 = 0
                        sum += l1
                        realChance.append(Decimal(l1))
                else:
                    # diameter_Nano.append(res['D_Median'])
                    realChance.append(Decimal(1))
                ##################
                dictNp[dm] = arrNp
                dictChp[dm] = realChance
                #######################
                if len(dictNp[dm]) != len(dictChp[dm]):
                    raise
                # npdp=[]
                # for i in range(len(dictNp[dm])):
                #    npdp.append(float(10000*dictNp[dm][i]*dictChp[dm][i]))
                # rescheck = self.LogNormalFit(x=dpArr, y=npdp)
                # k=3

            ################################
            ################################
            dpMedianDict = {}
            suggestedDict = {}
            for dm in self.arrMobilityDiamNano:
                dpMedianDict[dm] = res['D_Median']
                if res['Sigma_G'] != 1:
                    pointNumber = self.__Primary_Sigma_dm_CTE_Nt
                else:
                    pointNumber = 1
                Df = FN.CreateConstantArray(Number=self.__AGG_FRACTAL_DIMENSION_CENTER, HowMany=pointNumber)
                kf = FN.CreateConstantArray(Number=self.__AGG_FRACTAL_PREFACTOR_CENTER, HowMany=pointNumber)
                RI_Real = FN.CreateConstantArray(Number=self.__AGG_RI_REAL_CENTER, HowMany=pointNumber)
                RI_Imag = FN.CreateConstantArray(Number=self.__AGG_RI_IMAG_CENTER, HowMany=pointNumber)
                wavelength = FN.CreateConstantArray(Number=self.__AGG_WLENGTH_CENTER, HowMany=pointNumber)
                S = {}
                S['Df'] = Df
                S['Kf'] = kf
                S['RI_Real'] = RI_Real
                S['RI_Imag'] = RI_Imag
                S['Np'] = dictNp[dm]
                S['chance'] = dictChp[dm]
                S['wL'] = wavelength

                monometerParameter = []
                dpDecimal = []
                for dp in primaryDiam:
                    dpDecimal.append(round(Decimal(dp), 6))
                    monometerParameter.append(pi * dp / self.__AGG_WLENGTH_CENTER)
                S['dp'] = dpDecimal
                S['MP'] = monometerParameter
                suggestedDict[dm] = S

            dictCheckedNpdp = self.CheckDictWithBoundary(dict=suggestedDict, Beyond=self.__GoBeyondBound)

            if self.__PlotDetails:
                # resSuggestedR, resSuggestedC = self.RealNpDp(suggestedDict)
                resCheckedR = self._CalcRealNpDp(dict=dictCheckedNpdp, dmCalc=False, Np_multiple=True)
                # for dm in self.arrMobilityDiamNano:
                #   dpMedianDict[dm] = resCheckedR['D_Median']-1
                # resCheckedR, resCheckedC = self.RealNpDp(dictCheckedNpdp)
                # self.PlotdpDistribution(suggestedDict, resSuggestedR, resSuggestedC, dpMedianDict, "Trad-Suggested")
                self.PlotdpDistributionTot(dictCheckedNpdp, resCheckedR, dpMedianDict, "Total", res)

            return dpMedianDict, dictCheckedNpdp

        except Exception as e:
            logging.exception(e)
            raise

    def ConvertDictToArray(self, dict):
        try:
            conDict = {}
            for dm in dict:
                L = dict[dm]
                X = []
                for i in range(len(L['Df'])):
                    S = []
                    S.append(L['Df'][i])
                    S.append(L['Kf'][i])
                    S.append(L['RI_Real'][i])
                    S.append(L['RI_Imag'][i])
                    S.append(L['wL'][i])
                    S.append(L['dp'][i])
                    S.append(L['Np'][i])
                    X.append(S)
                conDict[dm] = X
            return conDict

        except Exception as e:
            logging.exception(e)
            raise

    def CheckDictWithBoundary(self, dict, Beyond):
        try:
            checked = {}
            for dm in dict:
                checking = BCheck(Beyond)
                checked[dm] = checking.CheckInputDict(inputDict=dict[dm])
            return checked

        except Exception as e:
            logging.exception(e)
            raise

    def CalcPrimaryParticleNumber(self):
        try:
            dictNp = {}
            for dm in self.arrMobilityDiamNano:
                arrNp = []
                dpArr = self.dict_dpMobDistribNano[dm]
                for dp in dpArr:
                    arrNp.append(self._calcPPN(dm=dm, dp=dp))
                dictNp[dm] = arrNp
            return dictNp

        except Exception as e:
            logging.exception(e)
            raise

    def CreateDictForEvaluation(self):
        try:
            suggestedDict = {}
            for dm in self.arrMobilityDiamNano:
                if self.__AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER != 1:
                    pointNumber = self.__Primary_Sigma_dm_CTE_Nt
                else:
                    pointNumber = 1
                Df = FN.CreateConstantArray(Number=self.__AGG_FRACTAL_DIMENSION_CENTER, HowMany=pointNumber)
                kf = FN.CreateConstantArray(Number=self.__AGG_FRACTAL_PREFACTOR_CENTER, HowMany=pointNumber)
                RI_Real = FN.CreateConstantArray(Number=self.__AGG_RI_REAL_CENTER, HowMany=pointNumber)
                RI_Imag = FN.CreateConstantArray(Number=self.__AGG_RI_IMAG_CENTER, HowMany=pointNumber)
                wavelength = FN.CreateConstantArray(Number=self.__AGG_WLENGTH_CENTER, HowMany=pointNumber)
                S = {}
                S['Df'] = Df
                S['Kf'] = kf
                S['RI_Real'] = RI_Real
                S['RI_Imag'] = RI_Imag
                S['Np'] = self.dict_NpMobDistrib[dm]
                S['chance'] = self.dict_dpMobDistribChance[dm]
                S['wL'] = wavelength

                monometerParameter = []
                dpDecimal = []
                for dp in self.dict_dpMobDistribNano[dm]:
                    dpDecimal.append(round(Decimal(dp), 6))
                    monometerParameter.append(pi * dp / self.__AGG_WLENGTH_CENTER)
                S['dp'] = dpDecimal
                S['MP'] = monometerParameter
                suggestedDict[dm] = S

            return suggestedDict
        except Exception as e:
            logging.exception(e)
            raise

    def CalcMobilityDiamBins(self):
        try:
            bound_D_Max = self.__AGG_MOBILITY_DIAMETER_CENTER_MAX
            bound_D_Min = self.__AGG_MOBILITY_DIAMETER_CENTER_MIN
            total_Number_Bins = self.__Sample_Sigma_Bins
            diameter_Nano = []

            D_Ratio = (bound_D_Max / bound_D_Min) ** (1 / (total_Number_Bins - 1))

            for i in range(0, total_Number_Bins):
                d1 = bound_D_Min * (D_Ratio ** i)
                diameter_Nano.append(round(d1, 3))

            return sorted(diameter_Nano, key=float)

        except Exception as e:
            logging.exception(e)
            raise

    def CalcPrimDiamBinAtEachMob(self):
        try:
            primaryDiam = {}
            primaryChance = {}

            lengthq = len(self.arrMobilityDiamNano)
            qu = int(lengthq / self.__dmSelectForPlot)
            self.dmRandom = []
            for i in range(1, self.__dmSelectForPlot):
                self.dmRandom.append(self.arrMobilityDiamNano[3 + qu * i])

            for dm in self.arrMobilityDiamNano:

                dpMedian = self.dict_dpMedianNano[dm]

                bound_D_Max = dpMedian * (self.__AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER ** self.__Primary_Sigma_dm_CTE_Bound)
                bound_D_Min = dpMedian * (self.__AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER ** (-1 * self.__Primary_Sigma_dm_CTE_Bound))
                total_Number_Bins = self.__Primary_Sigma_dm_CTE_Nt

                diameter_Nano = []
                logNormalPDF = []

                A = np.linspace(-1.2, 1.2, total_Number_Bins)
                B = []
                for i in A:
                    B.append(float(gmpy2.root(pow(i, 7), 5)))

                x_Down = dpMedian - bound_D_Min
                x_Up = bound_D_Max - dpMedian

                for i in range(len(B)):
                    if B[i] < 0:
                        diameter_Nano.append(dpMedian - x_Down * abs(B[i]))
                    if B[i] > 0:
                        diameter_Nano.append(dpMedian + x_Up * abs(B[i]))

                diameter_Nano = sorted(diameter_Nano, key=float)
                diameter_Nano = gaussian_filter1d(diameter_Nano, 8)
                diameter_Nano = diameter_Nano.tolist()

                sum = 0

                if self.__AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER != 1:
                    for i in range(0, total_Number_Bins):
                        if i != 0:
                            l1 = self._calcLogNDistribPDF(dpMedian, self.__AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER, diameter_Nano[i], diameter_Nano[i - 1])
                        else:
                            l1 = 0
                        sum += l1
                        logNormalPDF.append(Decimal(l1))
                else:
                    diameter_Nano.append(dpMedian)
                    logNormalPDF.append(Decimal(1))

                primaryDiam[dm] = diameter_Nano
                primaryChance[dm] = logNormalPDF

            return primaryDiam, primaryChance

        except Exception as e:
            logging.exception(e)
            raise

    def CalcPrimaryParticleSizeMedianNano(self):
        try:
            dp100_nano = self._calc_dp100nm_nano()
            D_TEM = self._calc_D_TEM_From_EffDens_D_Alpha()
            self.infoDict['dp100_nano'] = dp100_nano
            self.infoDict['D_TEM'] = D_TEM
            dm_dp = self._calcPPSM(dp100_nano, D_TEM)
            return dm_dp

        except Exception as e:
            logging.exception(e)
            raise

    def _calc_D_TEM_From_EffDens_D_Alpha(self):
        try:
            Primary_D_Alpha = self.__AGG_EXPONENT_PROJECTED_AREA_COEFFICIENT_CENTER
            Eff_Dm = self.__AGG_EFF_DM_CENTER

            D_TEM = (2 * Primary_D_Alpha - Eff_Dm) / (2 * Primary_D_Alpha - 3)
            return D_TEM

        except Exception as e:
            logging.exception(e)
            raise

    def _calc_dp100nm_nano(self):
        try:
            Eff_rho_100nm = self.__AGG_EFF_RHO_100NM_CENTER
            Primary_K_Alpha = self.__AGG_PREFACTOR_PROJECTED_AREA_COEFFICIENT_CENTER
            Soot_Material_Density = self.__AGG_MATERIAL_DENSITY_CENTER
            Primary_D_Alpha = self.__AGG_EXPONENT_PROJECTED_AREA_COEFFICIENT_CENTER

            K = ((Eff_rho_100nm / (Primary_K_Alpha * Soot_Material_Density)) ** (1 / (3 - 2 * Primary_D_Alpha))) * (100)
            return K

        except Exception as e:
            logging.exception(e)
            raise

    def _calcPPN(self, dm, dp):
        try:
            N = self.__AGG_PREFACTOR_PROJECTED_AREA_COEFFICIENT_CENTER * ((dm / dp) ** (2 * self.__AGG_EXPONENT_PROJECTED_AREA_COEFFICIENT_CENTER))
            return round(Decimal(N), 3)

        except Exception as e:
            logging.exception(e)
            raise

    def _calcPPSM(self, dp100_nano, D_TEM):
        # Primary particle size for each mobility diameter
        try:
            dm_dp = {}
            for dm in self.arrMobilityDiamNano:
                dp = (dp100_nano * ((dm / 100) ** D_TEM))
                dm_dp[dm] = round(dp, 3)
            return dm_dp

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

    def SaveAndClosePlot(self, folderName, F1):
        try:
            if self.__OnlySVGSaving:

                Address = GF.GetAddressTo(folderName=self.__folderNameGraph + f"\{folderName}", fileName=F1, extension="svg")
                plt.savefig(Address, format='svg', dpi=self.__figureDPI, bbox_inches='tight')

            else:
                Address = GF.GetAddressTo(folderName=self.__folderNameGraph + f"\{folderName}", fileName=F1, extension="jpg")
                plt.savefig(Address, format='jpg', dpi=self.__figureDPI, bbox_inches='tight')

            plt.clf()
            plt.close()

        except Exception as e:
            logging.exception(e)
            raise
