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
            self.__Sample_Sigma_Bins = 4  # Number of bins
            self.__Primary_Sigma_dm_CTE_Bound = 3  # Number of Sigma G to cover
            self.__Primary_Sigma_dm_CTE_Nt = 4  # Number of bins
            ######################################################################################
            self.infoDict['MobilityBins'] = self.__Sample_Sigma_Bins
            self.infoDict['NumberOfSigma'] = self.__Primary_Sigma_dm_CTE_Bound
            self.infoDict['PrimaryBins'] = self.__Primary_Sigma_dm_CTE_Nt
            self.infoDict['FittedSigma'] = 0
            self.infoDict['FittedMedian'] = 0
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

            dictChecked = self.CheckDictWithBoundary(dict=dictSuggested)
            dictConvertedTR = self.ConvertDictToArray(dict=dictChecked)

            dictdpMedianRealDict, dictCheckedRealNpdp = self.CalcRealNpdpDistribution(checkedDict=dictChecked)
            dictConvertedRE = self.ConvertDictToArray(dict=dictCheckedRealNpdp)

            self.RDGTMCore(DB_Info=DB_Info, dictConverted=dictConvertedTR, dictChecked=dictChecked, mediandpDict=self.dict_dpMedianNano, fileAppend="TR")
            self.RDGTMCore(DB_Info=DB_Info, dictConverted=dictConvertedRE, dictChecked=dictCheckedRealNpdp, mediandpDict=dictdpMedianRealDict, fileAppend="RE")

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

    def RDGTMCore(self, DB_Info, dictConverted, dictChecked, mediandpDict, fileAppend):
        try:
            CC = CalculationCore.CalculationCore(infoDict=DB_Info, calcDict=dictConverted)
            dictOutputT1, dictOutputR1 = CC.Calc()

            TCSC = TotalCrossSectionCalculator.TotalCrossSectionCalculator(infoDict=self.infoDict, outTMatrix=dictOutputT1, outRDG=dictOutputR1, checkedDict=dictChecked,
                                                                           dict_dpMedianNano=mediandpDict)
            dfResult = TCSC.Calc()
            dfResult.to_csv(f"TMatrix_RDG_Result\{fileAppend}_{self.infoDict['AA_FileName']}", index=False)

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

    def CalcRealNpdpDistribution(self, checkedDict):
        try:
            dpArr = []
            NpChanceArr = []
            for dm in checkedDict:
                Np = checkedDict[dm]['Np']
                dp = checkedDict[dm]['dp']
                chance = checkedDict[dm]['chance']
                for i in range(len(dp)):
                    dpArr.append(dp[i])
                    NpChanceArr.append(chance[i] * Np[i])

            diameter_Nano = []
            NpChanceBinArr = []
            bound_D_Max = 100
            bound_D_Min = 5
            total_Number_Bins = 43
            D_Ratio = (bound_D_Max / bound_D_Min) ** (1 / (total_Number_Bins - 1))

            for i in range(0, total_Number_Bins):
                d1 = bound_D_Min * (D_Ratio ** i)
                diameter_Nano.append(round(d1, 3))
                NpChanceBinArr.append(0)

            for i in range(1, len(diameter_Nano)):
                for p in range(len(dpArr)):
                    if dpArr[p] < diameter_Nano[i] and dpArr[p] > diameter_Nano[i - 1]:
                        NpChanceBinArr[i] = NpChanceBinArr[i] + float(NpChanceArr[p])

            res = self.LogNormalFit(diameter_Nano, NpChanceBinArr)
            logging.info(f"Lognormal fit res:___D_G:{res['D_G']}___Sigma_G:{res['Sigma_G']}___Total_Conc:{res['Total_Conc']}___D_Median:{res['D_Median']}")
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
            # fitted = [0]
            # fittedV2 = [0]
            # for i in range(1, total_Number_Bins):
            #    fitted.append(res['Total_Conc'] * self._calcLogNDistribPDF(median=res['D_Median'], sigmaG=res['Sigma_G'], D2=diameter_Nano[i], D1=diameter_Nano[i - 1]))
            #    fittedV2.append(res['Total_Conc'] * self._calcLogNDistribPDF(median=res['D_G'], sigmaG=res['Sigma_G'], D2=diameter_Nano[i], D1=diameter_Nano[i - 1]))
            # plt.plot(diameter_Nano, NpChanceBinArr, label="bins")
            # plt.plot(diameter_Nano, fitted, label="D_Median")
            # plt.plot(diameter_Nano, fittedV2, label="D_G")
            # plt.legend()
            # plt.show()
            dictCheckedNpdp = self.CheckDictWithBoundary(dict=suggestedDict)
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

    def CheckDictWithBoundary(self, dict):
        try:
            checked = {}
            for dm in dict:
                checking = BCheck()
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
            dp100_nano = self._calcPrimaryDiameter100nm_nano()
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

    def _calcPrimaryDiameter100nm_nano(self):
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
