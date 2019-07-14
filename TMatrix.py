# Keyhan Babaee, https://github.com/KeyhanB
# V1.1
# July 2019
import AppFunctions as FN
from ConfigParserM import logging
from scipy.interpolate import griddata
import GeneralFunctions as GF
from math import pi
from decimal import Decimal


def TmatrixInterpolator(FullMainDB, MainDBUniques, ABS_MainDB, SCA_MainDB, TargetArray):
    try:

        Tolerance = [1, 2, 2, 2, 2]
        ABS_CS = []
        SCA_CS = []
        refinedInput, Area = TmatrixRefiner(TargetArray)
        for i in range(len(refinedInput)):
            input = refinedInput[i][:]
            XXX, index = FN.getToleratedArray(Array=FullMainDB, Input=input, Tolerance=Tolerance, uniques=MainDBUniques)
            ABS_Coeff_Full = GF.selectColumnsList(index, ABS_MainDB, Dimension=1)
            SCA_Coeff_Full = GF.selectColumnsList(index, SCA_MainDB, Dimension=1)
            ABS_Interpolated = griddata(XXX, ABS_Coeff_Full, input, rescale=True)
            SCA_Interpolated = griddata(XXX, SCA_Coeff_Full, input, rescale=True)
            ABS_CS.append(Decimal(ABS_Interpolated[0]) * Area[i])
            SCA_CS.append(Decimal(SCA_Interpolated[0]) * Area[i])
            # print(input, "---", TargetArray[i], "---", "Abs:", ABS_CS[i], "Sca:", SCA_CS[i], "Area:", Area[i], "\n")

        return ABS_CS, SCA_CS

    except Exception as e:
        logging.exception(e)
        raise


def TmatrixRefiner(RawPoints):
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
