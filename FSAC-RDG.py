# Source: https://github.com/keyhanB/FlareNet-Scattering-Absorption-Calculator

import numpy as np
from math import pi
from math import sin
from ConfigParserM import logging
from decimal import Decimal


def FSAC_RDG(Df, kf, R_RI, I_RL, WaveL_nm, dp, Np, Eff_Dm, Eff_rho_100nm, sigma_Each_Mobility_Diameter, Primary_K_Alpha, Primary_D_Alpha, Soot_Material_Density):
    try:

        # Primary_Sigma_da_CTE = sigma_Each_Mobility_Diameter
        # Eff_k = Effective_Density_K_FromRho100nm(Eff_Dm, Eff_rho_100nm)
        # Primary_D_TEM = DTEM_FromEffectiveDensity_D_Alpha(Primary_D_Alpha, Eff_Dm)
        # Primary_Diameter_100nm = Primary_dp100nm(Eff_rho_100nm, Primary_K_Alpha, Soot_Material_Density, Primary_D_Alpha)
        # alpha=2*pi*radius/lambda

        Theta_Number = 180  # Number of bins in Theta
        Phi_Number = 180  # Number of bins in Phi
        Theta_Start = 0  # Angle Start for Theta
        Theta_Finish = pi  # Angle Finish for Theta
        Phi_Start = 0  # Angle Start for Phi
        Phi_Finish = 2 * pi  # Angle Finish for Phi
        Theta_Radian = np.linspace(Theta_Start, Theta_Finish, num=Theta_Number)  # Theta List
        Phi_Radian = np.linspace(Phi_Start, Phi_Finish, Phi_Number)  # Phi List
        Theta_Diff = abs((Theta_Finish - Theta_Start) / Theta_Number)  # Delta Theta
        Phi_Diff = abs((Phi_Finish - Phi_Start) / Phi_Number)  # Delta Phi
        dp_meter = dp * Decimal(10 ** (-9))
        WaveL_meter = WaveL_nm * Decimal(10 ** (-9))
        Wave_Number = Decimal(2) * Decimal(pi) / WaveL_meter
        Soot_Refractive_Index = complex(R_RI, -1 * I_RL)
        Soot_Complex = (((Soot_Refractive_Index ** 2) - 1) / ((Soot_Refractive_Index ** 2) + 2))
        Soot_FM = (abs(Soot_Complex)) ** 2
        Soot_EM = Soot_Complex.imag
        ########### Total Absorption
        Absorption_Cross_Section_Agg_meter2 = RDG_Absorption(K=Wave_Number, N=Np, Dp=dp_meter, E=Soot_EM)  # within  aggregate
        Absorption_Cross_Section_Agg_um2 = Absorption_Cross_Section_Agg_meter2 * Decimal(10 ** (12))
        ########### Differential Scattering
        # qDp_Full = []
        # qDp_Temp = []
        Differential_Scattering_Cross_Section_Agg_meter = 0
        for t in range(Theta_Number):
            q = Scattering_Wave_Vector(WaveLength_meter=WaveL_meter, Theta_radian=Theta_Radian[t])
            # qDp_Temp.append(q * dp_meter)
            Differential_Scattering_Cross_Section_Agg = RDG_Def_Scattering(K=Wave_Number, N=Np, Dp=dp_meter, q=q, F=Soot_FM, D_RDG=Df, K_RDG=kf, Formula=2)
            Differential_Scattering_Cross_Section_Agg_meter += Diff_Integral_Phi(Differential_Scattering_Cross_Section_Agg, Phi_Radian, Theta_Radian[t], Theta_Diff, Phi_Diff)
        # qDp_Full.append(qDp_Temp)
        Differential_Scattering_Cross_Section_Agg_um2 = Differential_Scattering_Cross_Section_Agg_meter * Decimal(10 ** (12))

        ########### Total Scattering
        Scattering_Cross_Section_Total_Agg = RDG_Total_Scattering(K=Wave_Number, N=Np, Dp=dp_meter, F=Soot_FM, D_RDG=Df, K_RDG=kf, Formula=2)

        return Absorption_Cross_Section_Agg_um2, Differential_Scattering_Cross_Section_Agg_um2


    except Exception as e:
        logging.exception(e)
        raise


def RDG_Total_Scattering(K, N, Dp, F, D_RDG, K_RDG, Formula=1):
    try:

        if Formula == 1:  # Sorensen 2001
            Monomer_Total = (8 / 3) * pi * (K ** 4) * ((Dp / 2) ** 6) * F
            Rg = ((N / K_RDG) ** (1 / D_RDG)) * (Dp / 2)
            G = (1 + ((4 / (3 * D_RDG)) * ((K * Rg) ** 2))) ** (-0.5 * D_RDG)

            Aggregate_Total = (N ** 2) * Monomer_Total * G
            logging.debug(f"RDG_Total_Scattering,Monomer_Total={Monomer_Total},Rg={Rg},G={G},Aggregate_Total={Aggregate_Total}_Sorensen(2001)")
        if Formula == 2:  # Yang 2005
            Monomer_Total = (8 / 3) * pi * (K ** 4) * ((Dp / 2) ** 6) * F
            Rg = ((N / K_RDG) ** (1 / D_RDG)) * (Dp / 2)
            Betha = 3 * D_RDG / (8 * (K * Rg) ** 2)
            if Betha >= 1:
                G = 1 - (2 * ((K * Rg) ** 2) / 3)
            elif Betha < 1:
                G = ((Betha / 2) * (3 - 3 * Betha + 2 * (Betha ** 2))) - ((((K * Rg * Betha) ** 2) / 3) * (3 - 4 * Betha + 3 * (Betha ** 2))) + (
                        ((2 * K * Rg) ** (-1 * D_RDG)) * ((3 / (2 - D_RDG)) - (12 / ((6 - D_RDG) * (4 - D_RDG))) - (3 * (Betha ** (1 - D_RDG / 2)) * ((1 / (2 - D_RDG)) - (2 * Betha / (4 - D_RDG)) + (2 * (Betha ** 2) / (6 - D_RDG))))))

            Aggregate_Total = (N ** 2) * Monomer_Total * G
            logging.debug(f"RDG_Total_Scattering,Monomer_Total={Monomer_Total},Rg={Rg}, Betha={Betha},G={G},Aggregate_Total={Aggregate_Total}_Yang(2005)")
        return Aggregate_Total

    except Exception as e:
        logging.exception(e)
        raise


def Scattering_Wave_Vector(WaveLength_meter, Theta_radian):
    try:
        A = (4 * pi / WaveLength_meter) * sin(Theta_radian / 2)
        return A

    except Exception as e:
        logging.exception(e)
        raise


def Diff_Integral_Phi(Scattered_Theta, Phi_Radian, Theta_Rad, Theta_Diff, Phi_Diff, Formula=1):
    try:

        A = 0
        for phi in range(len(Phi_Radian)):
            A += Scattered_Theta * sin(Theta_Rad) * Theta_Diff * Phi_Diff * (1 - (((cos(Phi_Radian[phi])) * sin(Theta_Rad)) ** 2))
        return A

    except Exception as e:
        logging.exception(e)
        raise


def RDG_Structure_Factor(q, Rg, D_RDG, K_RDG, Formula=1, C=1):
    try:

        if Formula == 1:  # Sorensen 2001
            C = 1.35 / K_RDG
            Check = q * Rg
            if Check <= 1:
                A = 1 - ((Check ** 2) / 3)  # Guinier equation (Guinier 1939;Guinier et al. 1955; Teixeira 1986).
            elif Check > 1:
                A = C * (Check) ** (-1 * D_RDG)
        elif Formula == 2:  # Yang 2005
            A = (1 + ((8 * (q * Rg) ** 2) / (3 * D_RDG)) + (q * Rg) ** 8) ** (-1 * D_RDG / 8)

        return A

    except Exception as e:
        logging.exception(e)
        raise


def RDG_Def_Scattering(K, N, Dp, q, F, D_RDG, K_RDG, Formula=1):
    try:

        Monomer_Differential = F * ((Dp / 2) ** 6) * (K ** 4)
        Rg = ((N / K_RDG) ** (1 / D_RDG)) * (Dp / 2)
        S_q = RDG_Structure_Factor(q=q, Rg=Rg, D_RDG=D_RDG, K_RDG=K_RDG, Formula=Formula)
        Aggregate_Differential = (N ** 2) * Monomer_Differential * S_q

        logging.debug(f"RDG_Def_Scattering,Monomer_Differential={Monomer_Differential},Rg={Rg},S_q={S_q},Aggregate_Differential={Aggregate_Differential}")
        return Aggregate_Differential

    except Exception as e:
        logging.exception(e)
        raise


def RDG_Absorption(K, N, Dp, E):
    try:
        A = -1 * N * 4 * pi * K * E * ((Dp / 2) ** 3)
        return A

    except Exception as e:
        logging.exception(e)
        raise


def DTEM_FromEffectiveDensity_D_Alpha(Primary_D_Alpha, Eff_dm):
    try:
        K = (2 * Primary_D_Alpha - Eff_dm) / (2 * Primary_D_Alpha - 3)
        return K

    except Exception as e:
        logging.exception(e)
        raise


def Primary_dp100nm(Eff_rho_100nm, Primary_K_Alpha, Soot_Material_Density, Primary_D_Alpha):
    try:
        K = ((Eff_rho_100nm / (Primary_K_Alpha * Soot_Material_Density)) ** (1 / (3 - 2 * Primary_D_Alpha))) * (100 * 10 ** (-9))
        return K

    except Exception as e:
        logging.exception(e)
        raise


def Effective_Density_K_FromRho100nm(Eff_dm, Eff_rho_100nm):
    try:
        da = 100 * (10 ** -9)

        K = Eff_rho_100nm / (da ** (Eff_dm - 3))
        return K

    except Exception as e:
        logging.exception(e)
        raise
