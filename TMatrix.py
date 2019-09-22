# Keyhan Babaee, https://github.com/KeyhanB
# V1.2
# Aug 2019
import AppFunctions as FN
from ConfigParserModule import logging
from scipy.interpolate import griddata
import GeneralFunctions as GF
from math import pi
from decimal import Decimal
from DBManagement import MySQLManagement
from concurrent.futures import ThreadPoolExecutor
from functools import partial


class TMatrixCalculation:

    def __init__(self, DBInfo):
        try:
            logging.info("T-Matrix method started.")
            TMatrix_DB_Main_TableName = 'Raw_V1'
            DB = MySQLManagement(DBInfo)
            TMatrix_DB_Main_Column_Name, TMatrix_DB_Main_Data_Full = DB.ReadAllRowsfromTable(TableName=TMatrix_DB_Main_TableName)
            self.__TMatrix_DB_Main_Column_Name = GF.selectColumnsList(ColumnIndex=[1, 2, 3, 4, 5], List=TMatrix_DB_Main_Column_Name, Dimension=1)
            self.__TMatrix_DB_Main_SCT_Coeff_Full = GF.selectColumnsList(ColumnIndex=[7], List=TMatrix_DB_Main_Data_Full)
            self.__TMatrix_DB_Main_ABS_Coeff_Full = GF.selectColumnsList(ColumnIndex=[8], List=TMatrix_DB_Main_Data_Full)
            self.__TMatrix_DB_Main_Data_Full = GF.selectColumnsList(ColumnIndex=[1, 2, 3, 4, 5], List=TMatrix_DB_Main_Data_Full)
            self.__TMatrix_DB_Main_Unique_Values = FN.uniqueEntry(self.__TMatrix_DB_Main_Data_Full)

        except Exception as e:
            logging.exception(e)
            raise

    def TMatrixCalc(self, TMatrix_Planned_Input, thread=3):
        try:
            logging.info("T-Matrix calculation started.")
            division = thread
            TMatrix_Planned_Input_Chopped = GF.divideArray(NumberofDivisions=division, List=TMatrix_Planned_Input)
            func = partial(self.TMatrixInterpolator, self.__TMatrix_DB_Main_Data_Full, self.__TMatrix_DB_Main_Unique_Values, self.__TMatrix_DB_Main_ABS_Coeff_Full,
                           self.__TMatrix_DB_Main_SCT_Coeff_Full)
            arrInTMatrix = []
            arrOutTMatrix = []
            with ThreadPoolExecutor() as executor:
                for inputs, outputs in zip(TMatrix_Planned_Input_Chopped, executor.map(func, TMatrix_Planned_Input_Chopped)):
                    arrInTMatrix.append(inputs)
                    arrOutTMatrix.append(outputs)
            ####################################
            self.__TMatrix_Interpolation_Input = []
            self.__TMatrix_Interpolation_Output = []
            for i in range(division):
                for j in range(len(arrInTMatrix[i])):
                    self.__TMatrix_Interpolation_Input.append(arrInTMatrix[i][j])
                for j in range(len(arrOutTMatrix[i][0])):
                    arrT = []
                    arrT.append(arrOutTMatrix[i][0][j])
                    arrT.append(arrOutTMatrix[i][1][j])
                    self.__TMatrix_Interpolation_Output.append(arrT)
            logging.info("T-Matrix calculation finished.")
            return self.__TMatrix_Interpolation_Input, self.__TMatrix_Interpolation_Output

        except Exception as e:
            logging.exception(e)
            raise

    def TMatrixInterpolator(self, FullMainDB, MainDBUniques, ABS_MainDB, SCA_MainDB, TargetArray):
        try:

            Tolerance = [1, 3, 2, 2, 2]
            ABS_CS = []
            SCA_CS = []
            refinedInput, Area = self.TmatrixRefiner(TargetArray)
            for i in range(len(refinedInput)):
                input = refinedInput[i][:]
                XXX, index = FN.getToleratedArray(Array=FullMainDB, Input=input, Tolerance=Tolerance, uniques=MainDBUniques)
                ABS_Coeff_Full = GF.selectColumnsList(index, ABS_MainDB, Dimension=1)
                SCA_Coeff_Full = GF.selectColumnsList(index, SCA_MainDB, Dimension=1)
                ABS_Interpolated = griddata(XXX, ABS_Coeff_Full, input, rescale=True)
                SCA_Interpolated = griddata(XXX, SCA_Coeff_Full, input, rescale=True)
                ABS_CS.append(Decimal(ABS_Interpolated[0]) * Area[i])
                SCA_CS.append(Decimal(SCA_Interpolated[0]) * Area[i])

            return ABS_CS, SCA_CS

        except Exception as e:
            logging.exception(e)
            raise

    def TmatrixRefiner(self, RawPoints):
        try:
            rows = len(RawPoints)
            Area = []
            out = []
            for i in range(rows):
                A = []
                A.append(RawPoints[i][0])
                A.append(RawPoints[i][2])
                A.append(RawPoints[i][3])
                A.append(Decimal(pi) * RawPoints[i][5] / RawPoints[i][4])
                A.append(RawPoints[i][6])
                Area.append(Decimal(pi) * (((RawPoints[i][6] ** (Decimal(1 / 3))) * (RawPoints[i][5] / Decimal(1000))) ** Decimal(2)) / Decimal(4))
                out.append(A)

            return out, Area

        except Exception as e:
            logging.exception(e)
            raise
