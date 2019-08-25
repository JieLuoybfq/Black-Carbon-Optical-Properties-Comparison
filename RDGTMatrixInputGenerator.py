# Keyhan Babaee, https://github.com/KeyhanB
# V1.2
# Aug 2019
from ConfigParserM import logging
from math import log
from math import exp
from math import pi
import AppFunctions as FN
from BoundaryFinder import BCDBBoundaryCheck as BCheck
from decimal import Decimal
from TMatrix import TMatrixCalculation
from FSAC_RDG import RDGCalculation
from GraphFunctions import GraphMethods
import csv


class KeyhanV1:
    def __init__(self, inputDict):
        try:
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

            self.__Sample_Sigma_Bins = 50  # Number of bins
            self.__Primary_Sigma_dm_CTE_Bound = 3  # Number of Sigma G to cover
            self.__Primary_Sigma_dm_CTE_Nt = 30  # Number of bins

        except Exception as e:
            logging.exception(e)
            raise

    def D_TEM_FromEffDens_D_Alpha(self):
        try:
            Primary_D_Alpha = self.__AGG_EXPONENT_PROJECTED_AREA_COEFFICIENT_CENTER
            Eff_Dm = self.__AGG_EFF_DM_CENTER

            D_TEM = (2 * Primary_D_Alpha - Eff_Dm) / (2 * Primary_D_Alpha - 3)
            return D_TEM

        except Exception as e:
            logging.exception(e)
            raise

    def PrimaryDiameter100nm_nano(self):
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

    def KeyhanV1Calc(self, DB_Info):
        try:
            self.mobilityDiam_nano = self.MobilityDiamBinGenerator()
            self.dp_Median_nano = self.PrimaryParticleSizeMedian_nano()
            self.dp_Mobility_Distributed_nano, self.dp_Mobility_Distributed_chance = self.PrimaryDiamBinGenerator()
            self.Np_Mobility_Distributed = self.PrimaryParticleNumber()
            suggestedCalcDict = self.CreateDictForEvaluation()
            checkedDict = self.CheckDictWithBoundary(dict=suggestedCalcDict)
            convertedDict = self.ConvertDictToArray(Dict=checkedDict)
            T1 = TMatrixCalculation(DBInfo=DB_Info)
            R1 = RDGCalculation()
            inputT1 = {}
            outputT1 = {}
            inputR1 = {}
            outputR1 = {}
            errorReal = {}
            for dm in convertedDict:
                inputT1[dm], outputT1[dm] = T1.TMatrixCalc(TMatrix_Planned_Input=convertedDict[dm])
                inputR1[dm], outputR1[dm] = R1.RDGCalc(RDG_Planned_Input=convertedDict[dm])
            res_ABS_CRS_T1, res_SCA_CRS_T1, res_ABS_CRS_R1, res_SCA_CRS_R1 = self.CalcTotalCrossSection(outTMatrix=outputT1, outRDG=outputR1, checkedDict=checkedDict)
            A = self.ErrorEvaluation(ABS_CRS_T=res_ABS_CRS_T1, SCA_CRS_T=res_SCA_CRS_T1, ABS_CRS_R=res_ABS_CRS_R1, SCA_CRS_R=res_SCA_CRS_R1)
            csv_file = "Names.csv"
            with open(csv_file, 'w') as csvfile:
                writer = csv.writer(csvfile)
                # writer.writeheader()

                writer.writerows(A)

            # G=GraphMethods()
            # G.toSaveHistogram()
            b = 3

        except Exception as e:
            logging.exception(e)
            raise

    def ErrorEvaluation(self, ABS_CRS_T, SCA_CRS_T, ABS_CRS_R, SCA_CRS_R):
        try:
            out = []
            dm2 = []
            errorRealABS = {}
            errorRealSCA = {}
            errorAbsoluteABS = {}
            errorAbsoluteSCA = {}
            errorPercentABS = {}
            errorPercentSCA = {}
            errorNormalABS = {}
            errorNormalSCA = {}
            for dm in ABS_CRS_T:
                dm2.append(dm)
                errorRealABS[dm] = ABS_CRS_T[dm] - ABS_CRS_R[dm]
                errorRealSCA[dm] = SCA_CRS_T[dm] - SCA_CRS_R[dm]
                errorAbsoluteABS[dm] = abs(ABS_CRS_T[dm] - ABS_CRS_R[dm])
                errorAbsoluteSCA[dm] = abs(SCA_CRS_T[dm] - SCA_CRS_R[dm])
                if ABS_CRS_R[dm] != 0:
                    errorPercentABS[dm] = 100 * abs(ABS_CRS_T[dm] - ABS_CRS_R[dm]) / ABS_CRS_R[dm]
                    errorPercentSCA[dm] = 100 * abs(SCA_CRS_T[dm] - SCA_CRS_R[dm]) / SCA_CRS_R[dm]
                    errorNormalABS[dm] = ABS_CRS_T[dm] / ABS_CRS_R[dm]
                    errorNormalSCA[dm] = SCA_CRS_T[dm] / SCA_CRS_R[dm]
                else:
                    errorPercentABS[dm] = 0
                    errorPercentSCA[dm] = 0
                    errorNormalABS[dm] = 0
                    errorNormalSCA[dm] = 0

                out.append([errorRealABS[dm], errorRealSCA[dm], errorAbsoluteABS[dm], errorAbsoluteSCA[dm]
                               , errorPercentABS[dm], errorPercentSCA[dm], errorNormalABS[dm], errorNormalSCA[dm]])
            # out={'errorRealABS':errorRealABS,
            #      'errorRealSCA': errorRealSCA,
            #      'errorAbsoluteABS': errorAbsoluteABS,
            #      'errorAbsoluteSCA': errorAbsoluteSCA,
            #      'errorPercentABS': errorPercentABS,
            #      'errorPercentSCA':errorPercentSCA,
            #      'errorNormalABS': errorNormalABS,
            #      'errorNormalSCA': errorNormalSCA,
            #      }
            csv_file = "dm.csv"
            with open(csv_file, 'w') as csvfile:
                writer = csv.writer(csvfile)
                # writer.writeheader()
                writer.writerow(dm2)
            return out
        except Exception as e:
            logging.exception(e)
            raise

    def CalcTotalCrossSection(self, outTMatrix, outRDG, checkedDict):
        try:

            res_ABS_T = {}
            res_SCA_T = {}
            res_ABS_R = {}
            res_SCA_R = {}
            for dm in outTMatrix:
                ABS_T = 0
                SCA_T = 0
                ABS_R = 0
                SCA_R = 0
                for i in range(len(outTMatrix[dm])):
                    ABS_T += outTMatrix[dm][i][0] * checkedDict[dm]['chance'][i]
                    SCA_T += outTMatrix[dm][i][1] * checkedDict[dm]['chance'][i]
                    ABS_R += outRDG[dm][i][0] * checkedDict[dm]['chance'][i]
                    SCA_R += outRDG[dm][i][1] * checkedDict[dm]['chance'][i]
                res_ABS_T[dm] = ABS_T
                res_SCA_T[dm] = SCA_T
                res_ABS_R[dm] = ABS_R
                res_SCA_R[dm] = SCA_R

            return res_ABS_T, res_SCA_T, res_ABS_R, res_SCA_R
        except Exception as e:
            logging.exception(e)
            raise

    def ConvertDictToArray(self, Dict):
        try:
            conDict = {}
            for dm in Dict:
                L = Dict[dm]
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

    def PrimaryParticleNumber(self):
        try:
            Np = {}
            for dm in self.mobilityDiam_nano:
                NArr = []
                dpArr = self.dp_Mobility_Distributed_nano[dm]
                for dp in dpArr:
                    NArr.append(self.PPN_Calc(dm=dm, dp=dp))
                Np[dm] = NArr
            return Np

        except Exception as e:
            logging.exception(e)
            raise

    def PPN_Calc(self, dm, dp):
        try:
            N = self.__AGG_PREFACTOR_PROJECTED_AREA_COEFFICIENT_CENTER * ((dm / dp) ** (2 * self.__AGG_EXPONENT_PROJECTED_AREA_COEFFICIENT_CENTER))
            return round(Decimal(N), 3)

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

    def CreateDictForEvaluation(self):
        try:
            suggestedDict = {}
            for dm in self.mobilityDiam_nano:
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
                S['Np'] = self.Np_Mobility_Distributed[dm]
                S['chance'] = self.dp_Mobility_Distributed_chance[dm]
                S['wL'] = wavelength

                monometerParameter = []
                dpDecimal = []
                for dp in self.dp_Mobility_Distributed_nano[dm]:
                    dpDecimal.append(round(Decimal(dp), 3))
                    monometerParameter.append(pi * dp / self.__AGG_WLENGTH_CENTER)
                S['dp'] = dpDecimal
                S['MP'] = monometerParameter
                suggestedDict[dm] = S

            return suggestedDict
        except Exception as e:
            logging.exception(e)
            raise

    def MobilityDiamBinGenerator(self):
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

    def PrimaryDiamBinGenerator(self):
        try:
            primaryDiam = {}
            primaryChance = {}
            for dm in self.mobilityDiam_nano:

                dpMedian = self.dp_Median_nano[dm]

                bound_D_Max = dpMedian * (self.__AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER ** self.__Primary_Sigma_dm_CTE_Bound)
                bound_D_Min = dpMedian * (self.__AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER ** (-1 * self.__Primary_Sigma_dm_CTE_Bound))
                total_Number_Bins = self.__Primary_Sigma_dm_CTE_Nt
                diameter_Nano = []
                logNormalPDF = []
                # Keyhan = 10 ** 9
                sum = 0
                if self.__AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER != 1:
                    D_Ratio = (bound_D_Max / bound_D_Min) ** (1 / (total_Number_Bins - 1))
                    for i in range(0, total_Number_Bins):
                        d1 = bound_D_Min * (D_Ratio ** i)
                        diameter_Nano.append(round(d1, 3))
                        if i != 0:
                            l1 = self.LogN_Distribution(dpMedian, self.__AGG_POLYDISPERSITY_SIGMA_EACH_MOBILITY_CENTER, diameter_Nano[i], diameter_Nano[i - 1])
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

    def PrimaryParticleSizeMedian_nano(self):
        try:
            self.dp100_nano = self.PrimaryDiameter100nm_nano()
            self.D_TEM = self.D_TEM_FromEffDens_D_Alpha()
            dm_dp = self.PPSM_Calc()
            return dm_dp

        except Exception as e:
            logging.exception(e)
            raise

    def PPSM_Calc(self):
        try:
            dm_dp = {}
            for dm in self.mobilityDiam_nano:
                dp = (self.dp100_nano * ((dm / 100) ** self.D_TEM))
                dm_dp[dm] = round(dp, 3)
            return dm_dp

        except Exception as e:
            logging.exception(e)
            raise

    def LogN_Distribution(self, Median, SigmaG, Dp2, Dp1):  # return number between 0 to 1
        try:
            if SigmaG != 1:
                A = (1 / (log(SigmaG) * (2 * pi) ** 0.5)) * exp(-1 * (((log(Dp1) - log(Median)) ** 2) / (2 * (log(SigmaG)) ** 2))) * (log(Dp2) - log(Dp1))
            elif SigmaG == 1:
                if Dp2 >= Median and Median > Dp1:
                    return 1
                else:
                    return 0

            return A

        except Exception as e:
            logging.exception(e)
            raise
