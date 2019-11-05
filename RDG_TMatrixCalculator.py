# Keyhan Babaee, https://github.com/KeyhanB
# V1.2
# Oct 2019
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
            self.__Sample_Sigma_Bins = 39  # Number of bins
            self.__Primary_Sigma_dm_CTE_Bound = 3  # Number of Sigma G to cover
            self.__Primary_Sigma_dm_CTE_Nt = 36  # Number of bins
            ######################################################################################
            self.infoDict['MobilityBins'] = self.__Sample_Sigma_Bins
            self.infoDict['NumberOfSigma'] = self.__Primary_Sigma_dm_CTE_Bound
            self.infoDict['PrimaryBins'] = self.__Primary_Sigma_dm_CTE_Nt
            self.infoDict['FittedSigma'] = 0
            self.infoDict['FittedMedian'] = 0
            ######################################################################################
            self.__PlotDetails = False
            self.__dmSelectForPlot = 4  # Number of bins
            self.__folderNameGraph = 'Detail Graphs'
            self.__OnlySVGSaving = False
            self.__figureDPI = 400
            ######################################################################################
            self.__dmCalc = False
            self.__dmSigma = 1.4
            self.__dmMedian = 130
            self.__GoBeyondBound = True

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
            dictConvertedRE = self.ConvertDictToArray(dict=dictCheckedRealNpdp)

            self.RDGTMCore(DB_Info=DB_Info, dictConverted=dictConvertedTR, dictChecked=dictChecked, mediandpDict=self.dict_dpMedianNano, fileAppend="TR", RDGActive=True, TmatrixActive=False)
            self.RDGTMCore(DB_Info=DB_Info, dictConverted=dictConvertedRE, dictChecked=dictCheckedRealNpdp, mediandpDict=dict_dpMedianReal, fileAppend="RE", RDGActive=True, TmatrixActive=False)

            # newInfoDf = pd.DataFrame([self.infoDict])
            # newInfoDf.to_csv(f"TMatrix_RDG_Result\Beacon.csv", index=False)
            dfInfoDB = pd.read_csv(f"TMatrix_RDG_Result\Beacon.csv")
            dfInfoDB.loc[len(dfInfoDB)] = self.infoDict
            dfInfoDB.to_csv(f"TMatrix_RDG_Result\Beacon.csv", index=False)
            j = 1

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

    def PlotdpDistribution(self, dict, resR, resC, dictMedian, title):
        try:

            resR['D_Median'] = round(resR['D_Median'], 1)
            resR['Sigma_G'] = round(resR['Sigma_G'], 2)
            resC['D_Median'] = round(resC['D_Median'], 1)
            resC['Sigma_G'] = round(resC['Sigma_G'], 2)
            yLimits = [1, 100]
            xLimits = [49, 1500]

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
            ax1[0].plot(self.arrMobilityDiamNano, dpMedian, label="MedianDiameter")
            ax1[0].set_xscale('log')
            factor = 1
            dm_Base = self.dmRandom[0]

            for dm in self.dmRandom:
                l1 = len(dict[dm]['Df'])
                if l1 == 1:
                    ##### one point of random
                    ax1[0].plot([dm], [dict[dm]['dp'][0]], label=round(dm), color='red', marker='*', markersize=9)
                else:
                    dp = dict[dm]['dp']
                    dpChance = dict[dm]['chance']
                    Np = dict[dm]['Np']
                    x_adj = []
                    x_base = []
                    for i in range(len(dpChance)):
                        x_adj.append(dm - (float(dpChance[i]) * float(Np[i]) * (dm_Base / dm) * factor))
                        x_base.append(dm)
                    ax1[0].plot(x_adj, dp, label=round(dm), color='green')
                    ax1[0].plot(x_base, dp, color='yellow')

            #### total
            logNormalPDF = []
            last_dm = self.arrMobilityDiamNano[-1] * 1.5
            sum = 0

            if resR['Sigma_G'] != 1:
                for i in range(0, total_Number_Bins):
                    if i != 0:
                        l1 = self._calcLogNDistribPDF(resR['D_Median'], resR['Sigma_G'], dp_Nano_tot[i], dp_Nano_tot[i - 1])
                    else:
                        l1 = 0
                    sum += l1
                    logNormalPDF.append(l1)
            else:
                dp_Nano_tot1 = []
                dp_Nano_tot1.append(resR['D_Median'])
                logNormalPDF.append(1)

            x_adj = []
            x_base = []

            if resR['Sigma_G'] != 1:
                for i in range(len(logNormalPDF)):
                    x_adj.append(last_dm - (float(logNormalPDF[i]) * (last_dm / dm_Base) * 300))
                    x_base.append(last_dm)
                ax1[0].plot(x_adj, dp_Nano_tot, label='Total', color='black')
                ax1[0].plot(x_base, dp_Nano_tot, color='yellow')
            else:
                x_adj = last_dm
                ax1[0].plot(x_adj, dp_Nano_tot1, label='Total', color='black', marker='*', markersize=9)
            ax1[0].set_ylim(yLimits)
            ax1[0].set_xlim(xLimits)
            ax1[0].set_title(f"Considering [Np*Chance], Median:{round(resR['D_Median'], 1)}, Sigma:{round(resR['Sigma_G'], 2)}")
            ax1[0].legend()
            ##########################################
            ##########################################
            ########## only chance
            ax1[1].plot(self.arrMobilityDiamNano, dpMedian, label="MedianDiameter")
            ax1[1].set_xscale('log')
            factor = 300
            dm_Base = self.dmRandom[0]

            for dm in self.dmRandom:
                l1 = len(dict[dm]['Df'])
                if l1 == 1:
                    ax1[1].plot([dm], [dict[dm]['dp'][0]], label=round(dm), color='red', marker='*', markersize=9)
                else:
                    dp = dict[dm]['dp']
                    dpChance = dict[dm]['chance']
                    # Np = dict[dm]['Np']
                    x_adj = []
                    x_base = []
                    for i in range(len(dpChance)):
                        x_adj.append(dm - (float(dpChance[i]) * (dm / dm_Base) * factor))
                        x_base.append(dm)
                    ax1[1].plot(x_adj, dp, label=round(dm), color='green')
                    ax1[1].plot(x_base, dp, color='yellow')

            logNormalPDF = []
            last_dm = self.arrMobilityDiamNano[-1] * 1.5
            sum = 0
            if resC['Sigma_G'] != 1:
                for i in range(0, total_Number_Bins):
                    if i != 0:
                        l1 = self._calcLogNDistribPDF(resC['D_Median'], resC['Sigma_G'], dp_Nano_tot[i], dp_Nano_tot[i - 1])
                    else:
                        l1 = 0
                    sum += l1
                    logNormalPDF.append(l1)
            else:
                dp_Nano_tot1 = []
                dp_Nano_tot1.append(resC['D_Median'])
                logNormalPDF.append(1)

            x_adj = []
            x_base = []
            if resC['Sigma_G'] != 1:
                for i in range(len(logNormalPDF)):
                    x_adj.append(last_dm - (float(logNormalPDF[i]) * (last_dm / dm_Base) * 300))
                    x_base.append(last_dm)
                ax1[1].plot(x_adj, dp_Nano_tot, label='Total', color='black')
                ax1[1].plot(x_base, dp_Nano_tot, color='yellow')
            else:
                x_adj = last_dm
                ax1[1].plot(x_adj, dp_Nano_tot1, label='Total', color='black', marker='*', markersize=9)
                # ax1[0].plot(x_base, diameter_Nano_tot1, color='yellow')

            ax1[1].set_ylim(yLimits)
            ax1[1].set_xlim(xLimits)
            ax1[1].set_title(f"Considering Chance Only, Median:{round(resC['D_Median'], 1)}, Sigma:{round(resC['Sigma_G'], 2)}")
            ax1[1].legend()
            ########################
            ########################
            ########################
            fig.suptitle(f"{title}_Dm:{self.__AGG_EFF_DM_CENTER}_rhoeff:{self.__AGG_EFF_RHO_100NM_CENTER}_Sigma_p:{self.__AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER}")
            self.SaveAndClosePlot(
                folderName=f"Dp_Dm_Distribution_Dm_{self.__AGG_EFF_DM_CENTER}_rhoeff_{self.__AGG_EFF_RHO_100NM_CENTER}_Sigma_p_{self.__AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER}",
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

    def RealNpDp(self, dict):
        try:
            #################################################

            if self.__dmCalc:
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
            bound_D_MaxCC = 100
            bound_D_MinCC = 1
            total_Number_BinsCC = self.__Primary_Sigma_dm_CTE_Nt
            #################################################
            #################################################
            ############ chance multiplied
            dpArr = []
            NpChanceArr = []
            counter = 0
            for dm in dict:
                if self.__dmCalc == True:
                    dmChance = logNormalPDFdm[counter] * 100
                    counter += 1
                else:
                    dmChance = 1

                l1 = len(dict[dm]['Df'])
                if l1 == 1:
                    Np = dict[dm]['Np'][:1]
                    dp = dict[dm]['dp'][:1]
                    chance = dict[dm]['chance'][:1]
                    for i in range(len(dp)):
                        dpArr.append(dp[i])
                        NpChanceArr.append(chance[i] * Np[i] * Decimal(dmChance))
                else:
                    Np = dict[dm]['Np']
                    dp = dict[dm]['dp']
                    chance = dict[dm]['chance']
                    for i in range(len(dp)):
                        dpArr.append(dp[i])
                        NpChanceArr.append(chance[i] * Np[i] * Decimal(dmChance))

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
                        NpChanceBinArr[i] = NpChanceBinArr[i] + float(NpChanceArr[p])

            resR = self.LogNormalFit(diameter_Nano, NpChanceBinArr)
            logging.info(f"LogNormal fit [Chance*Np] res:___D_G:{resR['D_G']}___Sigma_G:{resR['Sigma_G']}___Total_Conc:{resR['Total_Conc']}___D_Median:{resR['D_Median']}")
            #########################
            #########################
            ######################### Only Chance
            dpArrC = []
            NpChanceArrC = []
            counter = 0
            for dm in dict:
                if self.__dmCalc == True:
                    dmChance = logNormalPDFdm[counter] * 100
                    counter += 1
                else:
                    dmChance = 1

                l1 = len(dict[dm]['Df'])
                if l1 == 1:
                    # Np = dict[dm]['Np'][:1]
                    dp = dict[dm]['dp'][:1]
                    chance = dict[dm]['chance'][:1]
                    for i in range(len(dp)):
                        dpArrC.append(dp[i])
                        NpChanceArrC.append(chance[i] * Decimal(dmChance))
                else:
                    # Np = dict[dm]['Np']
                    dp = dict[dm]['dp']
                    chance = dict[dm]['chance']
                    for i in range(len(dp)):
                        dpArrC.append(dp[i])
                        NpChanceArrC.append(chance[i] * Decimal(dmChance))

            diameter_NanoC = []
            NpChanceBinArrC = []

            D_Ratio = (bound_D_MaxCC / bound_D_MinCC) ** (1 / (total_Number_BinsCC - 1))

            for i in range(0, total_Number_BinsCC):
                d1 = bound_D_MinCC * (D_Ratio ** i)
                diameter_NanoC.append(round(d1, 3))
                NpChanceBinArrC.append(0)

            for i in range(1, len(diameter_NanoC)):
                for p in range(len(dpArrC)):
                    if dpArrC[p] < diameter_NanoC[i] and dpArrC[p] > diameter_NanoC[i - 1]:
                        NpChanceBinArrC[i] = NpChanceBinArrC[i] + float(NpChanceArrC[p])

            resC = self.LogNormalFit(diameter_NanoC, NpChanceBinArrC)
            logging.info(f"LogNormal fit [Only Chance] res:___D_G:{resC['D_G']}___Sigma_G:{resC['Sigma_G']}___Total_Conc:{resC['Total_Conc']}___D_Median:{resC['D_Median']}")

            return resR, resC
        except Exception as e:
            logging.exception(e)
            raise

    def CalcRealNpdpDistribution(self, dictChecked, dictSuggested):
        try:
            resSuggestedR, resSuggestedC = self.RealNpDp(dictSuggested)
            resCheckedR, resCheckedC = self.RealNpDp(dictChecked)
            if self.__PlotDetails:
                self.PlotdpDistribution(dictSuggested, resSuggestedR, resSuggestedC, self.dict_dpMedianNano, "Obs-Suggested")
                self.PlotdpDistribution(dictChecked, resCheckedR, resCheckedC, self.dict_dpMedianNano, "Obs-Checked")

            ############################
            res = resCheckedR
            ############################
            self.infoDict['FittedSigma'] = res['Sigma_G']
            self.infoDict['FittedMedian'] = res['D_Median']
            ############################
            ############################
            ############################
            bound_D_Max = res['D_Median'] * (res['Sigma_G'] ** self.__Primary_Sigma_dm_CTE_Bound)
            bound_D_Min = res['D_Median'] * (res['Sigma_G'] ** (-1 * self.__Primary_Sigma_dm_CTE_Bound))
            total_Number_Bins = self.__Primary_Sigma_dm_CTE_Nt

            diameter_Nano = []
            logNormalPDF = []

            A = np.linspace(-1.2, 1.2, total_Number_Bins)
            B = []
            for i in A:
                B.append(float(gmpy2.root(pow(i, 7), 5)))

            x_Down = res['D_Median'] - bound_D_Min
            x_Up = bound_D_Max - res['D_Median']

            for i in range(len(B)):
                if B[i] < 0:
                    diameter_Nano.append(res['D_Median'] - x_Down * abs(B[i]))
                if B[i] > 0:
                    diameter_Nano.append(res['D_Median'] + x_Up * abs(B[i]))

            diameter_Nano = sorted(diameter_Nano, key=float)
            diameter_Nano = gaussian_filter1d(diameter_Nano, 8)
            diameter_Nano = diameter_Nano.tolist()

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
                diameter_Nano.append(res['D_Median'])
                logNormalPDF.append(Decimal(1))

            primaryDiam = diameter_Nano
            primaryChance = logNormalPDF
            ################################
            ################################
            dictNp = {}
            for dm in self.arrMobilityDiamNano:
                arrNp = []
                dpArr = primaryDiam
                for dp in dpArr:
                    arrNp.append(self._calcPPN(dm=dm, dp=dp))
                dictNp[dm] = arrNp
            ### add something

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
                S['chance'] = primaryChance
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
                resSuggestedR, resSuggestedC = self.RealNpDp(suggestedDict)
                resCheckedR, resCheckedC = self.RealNpDp(dictCheckedNpdp)
                self.PlotdpDistribution(suggestedDict, resSuggestedR, resSuggestedC, dpMedianDict, "Trad-Suggested")
                self.PlotdpDistribution(dictCheckedNpdp, resCheckedR, resCheckedC, dpMedianDict, "Trad-Checked")
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
