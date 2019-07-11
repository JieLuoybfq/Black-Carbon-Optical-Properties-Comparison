# Source: https://github.com/keyhanB/FlareNet-Scattering-Absorption-Calculator

import AppFunctions as FN
import os
import numpy as np
from math import pi
from math import log
from ConfigParserM import logging

####################### Aggregate Distribution Setting
def RDG(Array):
    try:

        ####################### Plot options
        Sample_Diameter_Min = 50  # min diameter in nm for plotting
        Sample_Diameter_Max = 2100  # max diameter in nm for plotting
        Primary_Y_Min = 5  # in nm for plotting
        Primary_Y_Max = 50  # in nm for plotting

        ########### Labels

        X_LabelA1 = 'Mobility-equivalent Diameter (nm)'
        Y_LabelA1 = 'dN/dLogDm, (#/cm' + "$^{}$".format(3) + ')'
        Y_LabelA2 = "Primary Particle Diameter-dp(nm)"
        Y_LabelA3 = 'Absorption Cross Section(m' + "$^{}$".format(2) + '/dLogDm' + ')'
        Y_LabelA4 = "Absorption Efficiency"
        Y_LabelA5 = "Scattering Efficiency"
        Y_LabelA6 = 'Scattering Cross Section(m' + "$^{}$".format(2) + '/dLogDm' + ')'
        Y_LabelA7 = "Sample Mass(kg/dLogDm)"
        Y_LabelA8 = 'Effective Density (kg/m' + "$^{}$".format(3) + ')'
        Y_LabelA9 = "SSA"

        Y_LabelB1 = "Aggregate Mass Average1"
        Y_LabelB2 = "Aggregate Mass Average2"
        Y_LabelB3 = 'MAC(m' + "$^{}$".format(2) + '/g' + ')'
        Y_LabelB4 = 'MSC(m' + "$^{}$".format(2) + '/g' + ')'
        logging.info("Labels Created!")

        ####################### Saving Parameters

        Save_Situation = []
        Save_Diameter = {}
        Save_LogN_Sample_SizeD = {}
        Save_Soot_Primary_Diameter_Median_Nano = {}
        Save_Absorption_Cross_Section_Sample_dlnDp = {}
        Save_Primary_Diameter_Bank = {}
        Save_Primary_Number_Bank = {}
        Save_Primary_Probability_Bank = {}
        Save_Absorption_Efficiency_Sample = {}
        Save_Differential_Scattering_Cross_Section_Full_dlnDp = {}
        Save_Scattering_Efficiency_Sample = {}
        Save_SSA_Distribution = {}
        Save_SSA = {}
        Save_Scattering_Cross_Section_Total_Distribution_dlnDp = {}
        Save_Effective_Density = {}
        Save_Mass_Sample_1 = {}
        Save_Mass_Sample_2 = {}
        Save_MAC = {}
        Save_MSC_Diff = {}
        Save_MAC_Distribution = {}
        Save_MSC_Distribution = {}
        Save_MSC_Total = {}
        Counter = 0
        ####################### Main Loops

        for i1 in range(Sample_LogN_D_Median_Bins):  # Sample Median Loop
            for i2 in range(Sample_LogN_Sigma_Bins):  # Sample Sigma Loop
                for i3 in range(Eff_dm_Bins):  # Effective Density Dm Loop
                    for i4 in range(Eff_rho_100nm_Bins):  # Effective Density for 100nm aggregate Loop
                        for i6 in range(Primary_Sigma_da_CTE_Bins):  # Primary Diameter Sigma Loop

                            Sample_D_Median = Sample_D_Median_List[i1]
                            Sample_Sigma = Sample_Sigma_List[i2]
                            Eff_dm = Eff_dm_List[i3]
                            Eff_rho_100nm = Eff_rho_100nm_List[i4]
                            Primary_Sigma_da_CTE = Primary_Sigma_da_CTE_List[i6]

                            Eff_k = FN.Effective_Density_K_FromRho100nm(Eff_dm, Eff_rho_100nm)
                            Primary_D_TEM = FN.DTEM_FromEffectiveDensity_D_Alpha(Primary_D_Alpha, Eff_dm)
                            Primary_Diameter_100nm = FN.Primary_dp100nm(Eff_rho_100nm, Primary_K_Alpha, Soot_Material_Density, Primary_D_Alpha)

                            Counter += 1

                            Situation = f"S_Median Diameter={round(Sample_D_Median, 1)} (nm)- S_Sigma={round(Sample_Sigma, 1)}- Dm={round(Eff_dm, 2)}- k={round(Eff_k, 2)}- rho_100nm={round(Eff_rho_100nm, 1)}- D_TEM={round(Primary_D_TEM, 2)}- P100nm_Diam={round(Primary_Diameter_100nm * 10 ** 9, 1)} (nm)- P_Sigma={round(Primary_Sigma_da_CTE, 2)}"
                            Save_Situation.append(Situation)
                            Graph_Folder_Situation = Graph_Folder + "/" + Situation

                            Diameter_Meter, Diameter_Log, Diameter_Nano = FN.Bins_LogN_Distributed(Median_Diameter=Sample_D_Median * (10 ** (-9)), Sigma_G=Sample_Sigma, Sigma_G_Bound=Sample_Sigma_Bound, Total_Number_Bins=Sample_Sigma_Bins)

                            LogN_Sample_SizeDistribution_Plain = []
                            LogN_Sample_SizeDistribution = []

                            Save_Diameter[Situation] = Diameter_Nano[:-1]

                            SUM = 0
                            SUM1 = 0
                            logging.info(f"Sample's particle bins generated:{Situation}:Sigma Bound:{Sample_Sigma_Bound}: {Diameter_Nano}")

                            ########## Number Concentration

                            for k in range(Sample_Sigma_Bins - 1):
                                LogN_Sample_SizeDistribution_Plain.append(FN.LogN_Distribution(Median=Sample_D_Median, SigmaG=Sample_Sigma, Dp2=Diameter_Nano[k + 1], Dp1=Diameter_Nano[k]))
                                LogN_Sample_SizeDistribution.append(LogN_Sample_SizeDistribution_Plain[k] * Sample_Total_Number_Concentration)
                                SUM += LogN_Sample_SizeDistribution_Plain[k]
                                SUM1 += LogN_Sample_SizeDistribution[k]

                            logging.info(f"Sample's particle bins populated:{Situation}:Sigma Bound:{Sample_Sigma_Bound}: {LogN_Sample_SizeDistribution_Plain}")
                            logging.info(f"LogNormal Check:{Situation}:Sigma Bound:{Sample_Sigma_Bound}: {SUM}: {SUM1}")

                            Save_LogN_Sample_SizeD[Situation] = LogN_Sample_SizeDistribution
                            if Figure_Enable == 1:
                                Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder_Situation, FileName=str(Counter) + "_" + "NumberConcentration", Extension="jpg")
                                FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address=Address, X_Array=Diameter_Nano[:-1], Y_array=LogN_Sample_SizeDistribution, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1, Y_Legend="Number Concentration,\nTotal Conc.= " + str(round(Sample_Total_Number_Concentration, 1)), Y_label1=Y_LabelA1,
                                                                       Plot_Title=Situation)

                            ########## Primary Particle Size

                            Soot_Primary_Diameter_Median_meter = []
                            Soot_Primary_Diameter_Median_Nano = []
                            for k in range(Sample_Sigma_Bins - 1):
                                Soot_Primary_Diameter_Median_meter.append(FN.Primary_Particle_Size_meter(dm_meter=Diameter_Meter[k], dp100_meter=Primary_Diameter_100nm, Dtem=Primary_D_TEM))
                                Soot_Primary_Diameter_Median_Nano.append(Soot_Primary_Diameter_Median_meter[k] * 10 ** 9)
                            logging.info(f"Primary Particle Diameter(m):{Situation}:Sigma Bound:{Sample_Sigma_Bound}: {Soot_Primary_Diameter_Median_Nano}")

                            Save_Soot_Primary_Diameter_Median_Nano[Situation] = Soot_Primary_Diameter_Median_Nano
                            if Figure_Enable == 1:
                                Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder_Situation, FileName=str(Counter) + "_" + "PrimaryPS", Extension="jpg")
                                FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address=Address, X_Array=Diameter_Nano[:-1], Y_array=Soot_Primary_Diameter_Median_Nano, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, Y_Min=Primary_Y_Min, Y_Max=Primary_Y_Max, X_Label=X_LabelA1, Y_Legend="Primary Particle Size", Y_label1=Y_LabelA2, Plot_Title=Situation)

                            ########## Finding Primary Particle Number and Size

                            Primary_Diameter_Bank = []
                            Primary_Number_Bank = []
                            Primary_Probability_Bank = []
                            Aggregate_Mass_Bank = []

                            SUM2 = []
                            for k in range(Sample_Sigma_Bins - 1):
                                SUM_Temp = 0
                                Primary_Diameter, Primary_Number, Primary_Probability, Aggregate_Mass_array = FN.Primary_LogN_Generator_Mass(dp_Median_meter=Soot_Primary_Diameter_Median_meter[k], Sigma_g_Primary=Primary_Sigma_da_CTE, Sigma_g_Number=Primary_Sigma_da_CTE_Bound, Number_Points=Primary_Sigma_da_CTE_Nt, da_meter=Diameter_Meter[k],
                                                                                                                                             Eff_dm=Eff_dm, Eff_k=Eff_k, rho_cte=Soot_Material_Density)
                                Primary_Diameter_Bank.append(Primary_Diameter)
                                Primary_Number_Bank.append(Primary_Number)
                                Primary_Probability_Bank.append(Primary_Probability)
                                Aggregate_Mass_Bank.append(Aggregate_Mass_array)

                                for m in range(len(Primary_Probability)):
                                    SUM_Temp += Primary_Probability[m]
                                SUM2.append(SUM_Temp)
                                logging.debug(f"Primary Particle Diameter(m):{Situation}:Sigma Bound:{Sample_Sigma_Bound}: {Primary_Diameter}")
                                logging.debug(f"Primary Particle Number:{Situation}:Sigma Bound:{Sample_Sigma_Bound}: {Primary_Number}")
                            # FN.Fig_Plot_3D_Show_XCte(Diameter_Nano[:-1], Primary_Diameter_Bank, Primary_Probability_Bank)
                            logging.info(f"Primary Particle Probability Check:{Situation}:Sigma Bound:{Sample_Sigma_Bound}: {SUM2}")

                            Save_Primary_Diameter_Bank[Situation] = Primary_Diameter_Bank
                            Save_Primary_Number_Bank[Situation] = Primary_Number_Bank
                            Save_Primary_Probability_Bank[Situation] = Primary_Probability_Bank

                            ########### Absorption RDG

                            Absorption_Cross_Section_Sample = []
                            Absorption_Cross_Section = 0
                            Absorption_Cross_Section_Sample_dlnDp = []

                            for k in range(Sample_Sigma_Bins - 1):
                                Absorption_Cross_Section_Specific_Particle_Size = 0
                                for p in range(Primary_Sigma_da_CTE_Nt - 1):
                                    Absorption_Cross_Section_Agg = FN.RDG_Absorption(K=Wave_Number, N=Primary_Number_Bank[k][p], Dp=Primary_Diameter_Bank[k][p], E=Soot_EM)  # within  aggregate
                                    Absorption_Cross_Section_Specific_Particle_Size += Absorption_Cross_Section_Agg * Primary_Probability_Bank[k][p]

                                ABS = Absorption_Cross_Section_Specific_Particle_Size * LogN_Sample_SizeDistribution[k]
                                # Absorption_Cross_Section_Sample_dlnDp.append(ABS / (log(Diameter_Nano[k + 1]) - log(Diameter_Nano[k])))
                                Absorption_Cross_Section_Sample_dlnDp.append(ABS / (LogN_Sample_SizeDistribution[k]))
                                Absorption_Cross_Section_Sample.append(ABS)
                                Absorption_Cross_Section += ABS

                            Save_Absorption_Cross_Section_Sample_dlnDp[Situation] = Absorption_Cross_Section_Sample_dlnDp
                            if Figure_Enable == 1:
                                Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder_Situation, FileName=str(Counter) + "_" + "ABSCross", Extension="jpg")
                                FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address=Address, X_Array=Diameter_Nano[:-1], Y_array=Absorption_Cross_Section_Sample_dlnDp, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1,
                                                                       Y_Legend="Absorption Cross Section=\n" + str(round(Absorption_Cross_Section * 10 ** 12, 2)) + " um" + "$^{}$".format(2), Y_label1=Y_LabelA3,
                                                                       Plot_Title=Situation)
                            ########### Absorption EFF

                            Absorption_Efficiency_Sample = []
                            # Absorption_Efficiency_Sample_dlnDp = []

                            for k in range(Sample_Sigma_Bins - 1):
                                Absorption_Efficiency_Sample.append(FN.Abs_Scatt_Eff(CrossSection=Absorption_Cross_Section_Sample[k] / LogN_Sample_SizeDistribution[k], dm=Diameter_Meter[k]))
                                # Absorption_Efficiency_Sample_dlnDp.append(Absorption_Efficiency_Sample[k])

                            Save_Absorption_Efficiency_Sample[Situation] = Absorption_Efficiency_Sample
                            if Figure_Enable == 1:
                                Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder_Situation, FileName=str(Counter) + "_" + "ABSEff", Extension="jpg")
                                FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address=Address, X_Array=Diameter_Nano[:-1], Y_array=Absorption_Efficiency_Sample, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1, Y_Legend="Absorption Efficiency", Y_label1=Y_LabelA4, Plot_Title=Situation)

                            ########### Scattering RDG
                            ########### Differential

                            Differential_Scattering_Cross_Section_Full = []
                            Differential_Scattering_Cross_Section_Full_dlnDp = []
                            Scattering_Cross_Section_Diff_Tot = 0
                            qDp_Full = []
                            for k in range(Sample_Sigma_Bins - 1):
                                qDp_Temp = []
                                Differential_Scattering_Cross_Section_T = 0
                                for t in range(Theta_Number):
                                    q = FN.Scattering_Wave_Vector(WaveLength_meter=Wave_Length, Theta_radian=Theta_Radian[t])
                                    qDp_Temp.append(q * Soot_Primary_Diameter_Median_meter[k])
                                    Differential_Scattering_Cross_Section_Specific_Particle_Size = 0
                                    for p in range(Primary_Sigma_da_CTE_Nt - 1):
                                        Differential_Scattering_Cross_Section_Agg = FN.RDG_Def_Scattering(K=Wave_Number, N=Primary_Number_Bank[k][p], Dp=Primary_Diameter_Bank[k][p], q=q, F=Soot_FM, D_RDG=Soot_Fractal_D_mc_RDG, K_RDG=Soot_Prefactor_k_mc_RDG, Formula=2)
                                        Differential_Scattering_Cross_Section_Specific_Particle_Size += Differential_Scattering_Cross_Section_Agg * Primary_Probability_Bank[k][p]

                                    Differential_Scattering_Cross_Section_T += FN.Diff_Integral_Phi(Differential_Scattering_Cross_Section_Specific_Particle_Size, Phi_Radian, Theta_Radian[t], Theta_Diff, Phi_Diff)

                                qDp_Full.append(qDp_Temp)
                                Diff_Scatter = Differential_Scattering_Cross_Section_T * LogN_Sample_SizeDistribution[k]
                                Differential_Scattering_Cross_Section_Full.append(Diff_Scatter)
                                # Differential_Scattering_Cross_Section_Full_dlnDp.append(Diff_Scatter / (log(Diameter_Nano[k + 1]) - log(Diameter_Nano[k])))
                                Differential_Scattering_Cross_Section_Full_dlnDp.append(Diff_Scatter / (LogN_Sample_SizeDistribution[k]))
                                Scattering_Cross_Section_Diff_Tot += Diff_Scatter

                            Save_Differential_Scattering_Cross_Section_Full_dlnDp[Situation] = Differential_Scattering_Cross_Section_Full_dlnDp
                            if Figure_Enable == 1:
                                Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder_Situation, FileName=str(Counter) + "_" + "ScatterDiffCross", Extension="jpg")
                                FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address=Address, X_Array=Diameter_Nano[:-1], Y_array=Differential_Scattering_Cross_Section_Full_dlnDp, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1,
                                                                       Y_Legend="Total Scattering Cross Section (Differential)=\n" + str(round(Scattering_Cross_Section_Diff_Tot * 10 ** 12, 2)) + " um" + "$^{}$".format(2), Y_label1=Y_LabelA6, Plot_Title=Situation)

                            ########### Scattering EFF

                            Scattering_Efficiency_Sample = []
                            # Scattering_Efficiency_Sample_dlnDp = []

                            for k in range(Sample_Sigma_Bins - 1):
                                Scattering_Efficiency_Sample.append(FN.Abs_Scatt_Eff(CrossSection=Differential_Scattering_Cross_Section_Full[k] / LogN_Sample_SizeDistribution[k], dm=Diameter_Meter[k]))
                                # Scattering_Efficiency_Sample_dlnDp.append(Scattering_Efficiency_Sample[k])

                            Save_Scattering_Efficiency_Sample[Situation] = Scattering_Efficiency_Sample
                            if Figure_Enable == 1:
                                Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder_Situation, FileName=str(Counter) + "_" + "Scattering_Eff", Extension="jpg")
                                FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address=Address, X_Array=Diameter_Nano[:-1], Y_array=Scattering_Efficiency_Sample, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1, Y_Legend="Scattering Efficiency", Y_label1=Y_LabelA5, Plot_Title=Situation)

                            ########### SSA Distribution

                            SSA = []

                            for k in range(Sample_Sigma_Bins - 1):
                                SSA.append(FN.SSA_Calculator(Scattering_CS=Differential_Scattering_Cross_Section_Full[k], Absorption_CS=Absorption_Cross_Section_Sample[k]))
                                # Scattering_Efficiency_Sample_dlnDp.append(Scattering_Efficiency_Sample[k])

                            Save_SSA_Distribution[Situation] = SSA
                            if Figure_Enable == 1:
                                Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder_Situation, FileName=str(Counter) + "_" + "SSA", Extension="jpg")
                                FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address=Address, X_Array=Diameter_Nano[:-1], Y_array=SSA, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1, Y_Legend="SSA", Y_label1=Y_LabelA9, Plot_Title=Situation)

                            ########### Total Scattering

                            Scattering_Cross_Section_Total_Distribution = []
                            Scattering_Cross_Section_Total_Distribution_dlnDp = []
                            Scattering_Cross_Section_Total = 0
                            for k in range(Sample_Sigma_Bins - 1):

                                Scattering_Cross_Section_Total_Specific_Particle_Size = 0
                                for p in range(Primary_Sigma_da_CTE_Nt - 1):
                                    Scattering_Cross_Section_Total_Agg = FN.RDG_Total_Scattering(K=Wave_Number, N=Primary_Number_Bank[k][p], Dp=Primary_Diameter_Bank[k][p], F=Soot_FM, D_RDG=Soot_Fractal_D_mc_RDG, K_RDG=Soot_Prefactor_k_mc_RDG, Formula=2)
                                    Scattering_Cross_Section_Total_Specific_Particle_Size += Scattering_Cross_Section_Total_Agg * Primary_Probability_Bank[k][p]

                                Total_Scatter = Scattering_Cross_Section_Total_Specific_Particle_Size * LogN_Sample_SizeDistribution[k]
                                Scattering_Cross_Section_Total_Distribution.append(Total_Scatter)
                                Scattering_Cross_Section_Total_Distribution_dlnDp.append(Total_Scatter / (log(Diameter_Nano[k + 1]) - log(Diameter_Nano[k])))
                                Scattering_Cross_Section_Total += Total_Scatter

                            Save_Scattering_Cross_Section_Total_Distribution_dlnDp[Situation] = Scattering_Cross_Section_Total_Distribution_dlnDp
                            if Figure_Enable == 1:
                                Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder_Situation, FileName=str(Counter) + "_" + "ScatterTotCross", Extension="jpg")
                                FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address=Address, X_Array=Diameter_Nano[:-1], Y_array=Scattering_Cross_Section_Total_Distribution_dlnDp, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1,
                                                                       Y_Legend="Scattering Cross Section (Total)=\n" + str(round(Scattering_Cross_Section_Total * 10 ** 12, 2)) + " um" + "$^{}$".format(2),
                                                                       Y_label1=Y_LabelA6, Plot_Title=Situation)

                            ########### MAC_MSC

                            Mass_Sample_1 = []
                            Mass_Sample_2 = []
                            Density = []
                            Total_Mass_1 = 0
                            Total_Mass_2 = 0

                            Aggregate_Mass_Ave1 = []
                            Aggregate_Mass_Ave2 = []
                            MAC_Distribution = []
                            MSC_Distribution = []

                            for k in range(Sample_Sigma_Bins - 1):
                                Density.append(FN.Effective_Density(K=Eff_k, Dm=Eff_dm, da=Diameter_Meter[k]))
                                Aggregate_Mass_Ave1.append(FN.Mass_Calc(rho=Density[k], da=Diameter_Meter[k]))
                                Mass_Aggregate_Ave = 0
                                for p in range(Primary_Sigma_da_CTE_Nt - 1):
                                    Mass_Aggregate_Ave += Aggregate_Mass_Bank[k][p] * Primary_Probability_Bank[k][p]
                                Aggregate_Mass_Ave2.append(Mass_Aggregate_Ave)

                                MAC_Distribution.append(Absorption_Cross_Section_Sample[k] / (Aggregate_Mass_Ave2[k] * LogN_Sample_SizeDistribution[k] * 1000))
                                MSC_Distribution.append(Differential_Scattering_Cross_Section_Full[k] / (Aggregate_Mass_Ave2[k] * LogN_Sample_SizeDistribution[k] * 1000))

                                Mass_Effective1 = Aggregate_Mass_Ave1[k] * LogN_Sample_SizeDistribution[k]
                                Mass_Effective2 = Aggregate_Mass_Ave2[k] * LogN_Sample_SizeDistribution[k]

                                Mass_Sample_1.append(Mass_Effective1 / (log(Diameter_Nano[k + 1]) - log(Diameter_Nano[k])))
                                Mass_Sample_2.append(Mass_Effective2 / (log(Diameter_Nano[k + 1]) - log(Diameter_Nano[k])))
                                Total_Mass_1 += Mass_Effective1
                                Total_Mass_2 += Mass_Effective2

                            Total_Mass_1 = Total_Mass_1 * 1000  # kg to g
                            Total_Mass_2 = Total_Mass_2 * 1000  # kg to g

                            MAC_1 = Absorption_Cross_Section / Total_Mass_1
                            MSC_1 = Scattering_Cross_Section_Total / Total_Mass_1
                            MSC_Diff_1 = Scattering_Cross_Section_Diff_Tot / Total_Mass_1
                            SSA = Scattering_Cross_Section_Diff_Tot / (Scattering_Cross_Section_Diff_Tot + Absorption_Cross_Section)

                            MAC_2 = Absorption_Cross_Section / Total_Mass_2
                            MSC_2 = Scattering_Cross_Section_Total / Total_Mass_2
                            MSC_Diff_2 = Scattering_Cross_Section_Diff_Tot / Total_Mass_2

                            Save_Effective_Density[Situation] = Density
                            Save_Mass_Sample_1[Situation] = Mass_Sample_1
                            Save_Mass_Sample_2[Situation] = Mass_Sample_2
                            Save_MAC_Distribution[Situation] = MAC_Distribution
                            Save_MSC_Distribution[Situation] = MSC_Distribution
                            Save_MAC[Situation] = str(MAC_1) + ", " + str(MAC_2)
                            Save_MSC_Total[Situation] = str(MSC_1) + ", " + str(MSC_2)
                            Save_MSC_Diff[Situation] = str(MSC_Diff_1) + ", " + str(MSC_Diff_2)
                            Save_SSA[Situation] = SSA

                            if Figure_Enable == 1:
                                Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder_Situation, FileName=str(Counter) + "_" + "EffectiveDensity", Extension="jpg")
                                FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address=Address, X_Array=Diameter_Nano[:-1], Y_array=Density, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1, Y_label1=Y_LabelA8, Plot_Title=Situation)

                                Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder_Situation, FileName=str(Counter) + "_" + "Aggregate_Mass_Ave1", Extension="jpg")
                                FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address=Address, X_Array=Diameter_Nano[:-1], Y_array=Aggregate_Mass_Ave1, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1, Y_label1=Y_LabelB1, Plot_Title=Situation)

                                Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder_Situation, FileName=str(Counter) + "_" + "Aggregate_Mass_Ave2", Extension="jpg")
                                FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address=Address, X_Array=Diameter_Nano[:-1], Y_array=Aggregate_Mass_Ave2, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1, Y_label1=Y_LabelB2, Plot_Title=Situation)

                                Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder_Situation, FileName=str(Counter) + "_" + "MAC_Distribution", Extension="jpg")
                                FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address=Address, X_Array=Diameter_Nano[:-1], Y_array=MAC_Distribution, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1, Y_label1=Y_LabelB3, Plot_Title=Situation)

                                Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder_Situation, FileName=str(Counter) + "_" + "MSC_Distribution", Extension="jpg")
                                FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address=Address, X_Array=Diameter_Nano[:-1], Y_array=MSC_Distribution, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1, Y_label1=Y_LabelB4, Plot_Title=Situation)

                                Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder_Situation, FileName=str(Counter) + "_" + "MassDistribution1", Extension="jpg")
                                FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address=Address, X_Array=Diameter_Nano[:-1], Y_array=Mass_Sample_1, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1, Y_Legend="MAC=" + str(round(MAC_1, 2)) + "\n" + "MSC=" + str(round(MSC_1, 2)) + "\n" + "MSC-Diff=" + str(round(MSC_Diff_1, 2)),
                                                                       Y_label1=Y_LabelA7,
                                                                       Plot_Title=Situation)

                                Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder_Situation, FileName=str(Counter) + "_" + "MassDistribution2", Extension="jpg")
                                FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address=Address, X_Array=Diameter_Nano[:-1], Y_array=Mass_Sample_2, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1,
                                                                       Y_Legend="MAC=" + str(round(MAC_2, 2)) + "\n" + "MSC=" + str(round(MSC_2, 2)) + "\n" + "MSC-Diff=" + str(round(MSC_Diff_2, 2)),
                                                                       Y_label1=Y_LabelA7, Plot_Title=Situation)

        ################# SAVING

        Address = FN.File_Pointer(Main=script_dir, FolderName=Results_Folder, FileName="Sample_Agggregate_Size", Extension="csv")
        FN.Dictionary_ToCSV(Address, Save_Diameter)

        Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="Sample_SizeDistribution", Extension="jpg")
        FN.Fig_Plot_Save_1Lines_X_Log_Y_Dictionary_Linear(Address=Address, Identifier=Save_Situation, X_Array=Save_Diameter, Y_array=Save_LogN_Sample_SizeD, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1, Y_label1=Y_LabelA1, Plot_Title="Number Concentration")
        Address = FN.File_Pointer(Main=script_dir, FolderName=Results_Folder, FileName="Sample_SizeDistribution", Extension="csv")
        FN.Dictionary_ToCSV(Address, Save_LogN_Sample_SizeD)

        Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="Soot_Primary_Diameter_Median", Extension="jpg")
        FN.Fig_Plot_Save_1Lines_X_Log_Y_Dictionary_Linear(Address=Address, Identifier=Save_Situation, X_Array=Save_Diameter, Y_array=Save_Soot_Primary_Diameter_Median_Nano, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1, Y_label1=Y_LabelA2, Plot_Title="Primary Particle Diameter")
        Address = FN.File_Pointer(Main=script_dir, FolderName=Results_Folder, FileName="Soot_Primary_Diameter_Median", Extension="csv")
        FN.Dictionary_ToCSV(Address, Save_Soot_Primary_Diameter_Median_Nano)

        Address = FN.File_Pointer(Main=script_dir, FolderName=Results_Folder, FileName="Primary_Diameter_Distribution", Extension="csv")
        FN.Dictionary_ToCSV(Address, Save_Primary_Diameter_Bank)
        Address = FN.File_Pointer(Main=script_dir, FolderName=Results_Folder, FileName="Primary_Number_Distribution", Extension="csv")
        FN.Dictionary_ToCSV(Address, Save_Primary_Number_Bank)
        Address = FN.File_Pointer(Main=script_dir, FolderName=Results_Folder, FileName="Primary_Probability_Distribution", Extension="csv")
        FN.Dictionary_ToCSV(Address, Save_Primary_Probability_Bank)

        Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="Absorption_Cross_Section_Sample", Extension="jpg")
        FN.Fig_Plot_Save_1Lines_X_Log_Y_Dictionary_Linear(Address=Address, Identifier=Save_Situation, X_Array=Save_Diameter, Y_array=Save_Absorption_Cross_Section_Sample_dlnDp, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1, Y_label1=Y_LabelA3, Plot_Title="Total Absorption Cross Section")
        Address = FN.File_Pointer(Main=script_dir, FolderName=Results_Folder, FileName="Absorption_Cross_Section_Sample", Extension="csv")
        FN.Dictionary_ToCSV(Address, Save_Absorption_Cross_Section_Sample_dlnDp)

        Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="Absorption_Efficiency_Sample", Extension="jpg")
        FN.Fig_Plot_Save_1Lines_X_Log_Y_Dictionary_Linear(Address=Address, Identifier=Save_Situation, X_Array=Save_Diameter, Y_array=Save_Absorption_Efficiency_Sample, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1, Y_label1=Y_LabelA4, Plot_Title="Absorption Efficiency")
        Address = FN.File_Pointer(Main=script_dir, FolderName=Results_Folder, FileName="Absorption_Efficiency_Sample", Extension="csv")
        FN.Dictionary_ToCSV(Address, Save_Absorption_Efficiency_Sample)

        Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="Differential_Scattering_Cross_Section", Extension="jpg")
        FN.Fig_Plot_Save_1Lines_X_Log_Y_Dictionary_Linear(Address=Address, Identifier=Save_Situation, X_Array=Save_Diameter, Y_array=Save_Differential_Scattering_Cross_Section_Full_dlnDp, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1, Y_label1=Y_LabelA6, Plot_Title="Scattering Cross Section_ Diff")
        Address = FN.File_Pointer(Main=script_dir, FolderName=Results_Folder, FileName="Differential_Scattering_Cross_Section", Extension="csv")
        FN.Dictionary_ToCSV(Address, Save_Differential_Scattering_Cross_Section_Full_dlnDp)

        Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="Scattering_Efficiency", Extension="jpg")
        FN.Fig_Plot_Save_1Lines_X_Log_Y_Dictionary_Linear(Address=Address, Identifier=Save_Situation, X_Array=Save_Diameter, Y_array=Save_Scattering_Efficiency_Sample, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1, Y_label1=Y_LabelA5, Plot_Title="Scattering Efficiency")
        Address = FN.File_Pointer(Main=script_dir, FolderName=Results_Folder, FileName="Scattering_Efficiency", Extension="csv")
        FN.Dictionary_ToCSV(Address, Save_Scattering_Efficiency_Sample)

        Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="SSA", Extension="jpg")
        FN.Fig_Plot_Save_1Lines_X_Log_Y_Dictionary_Linear(Address=Address, Identifier=Save_Situation, X_Array=Save_Diameter, Y_array=Save_SSA_Distribution, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1, Y_label1=Y_LabelA9, Plot_Title="SSA Distribution")
        Address = FN.File_Pointer(Main=script_dir, FolderName=Results_Folder, FileName="SSA", Extension="csv")
        FN.Dictionary_ToCSV(Address, Save_SSA_Distribution)

        Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="Scattering_Cross_Section_Total", Extension="jpg")
        FN.Fig_Plot_Save_1Lines_X_Log_Y_Dictionary_Linear(Address=Address, Identifier=Save_Situation, X_Array=Save_Diameter, Y_array=Save_Scattering_Cross_Section_Total_Distribution_dlnDp, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1, Y_label1=Y_LabelA6, Plot_Title="Scattering Cross Section_ Total")
        Address = FN.File_Pointer(Main=script_dir, FolderName=Results_Folder, FileName="Scattering_Cross_Section_Total", Extension="csv")
        FN.Dictionary_ToCSV(Address, Save_Scattering_Cross_Section_Total_Distribution_dlnDp)

        Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="Effective_Density", Extension="jpg")
        FN.Fig_Plot_Save_1Lines_X_Log_Y_Dictionary_Linear(Address=Address, Identifier=Save_Situation, X_Array=Save_Diameter, Y_array=Save_Effective_Density, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1, Y_label1=Y_LabelA8, Plot_Title="Effective Density")
        Address = FN.File_Pointer(Main=script_dir, FolderName=Results_Folder, FileName="Effective_Density", Extension="csv")
        FN.Dictionary_ToCSV(Address, Save_Effective_Density)

        Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="Mass_Sample_Effective_Density", Extension="jpg")
        FN.Fig_Plot_Save_1Lines_X_Log_Y_Dictionary_Linear(Address=Address, Identifier=Save_Situation, X_Array=Save_Diameter, Y_array=Save_Mass_Sample_1, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1, Y_label1=Y_LabelA7, Plot_Title="Mass Distribution using Effective Density")
        Address = FN.File_Pointer(Main=script_dir, FolderName=Results_Folder, FileName="Mass_Sample_Effective_Density", Extension="csv")
        FN.Dictionary_ToCSV(Address, Save_Mass_Sample_1)

        Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="MAC_Distribution", Extension="jpg")
        FN.Fig_Plot_Save_1Lines_X_Log_Y_Dictionary_Linear(Address=Address, Identifier=Save_Situation, X_Array=Save_Diameter, Y_array=Save_MAC_Distribution, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1, Y_label1=Y_LabelB3, Plot_Title="MAC Distribution")
        Address = FN.File_Pointer(Main=script_dir, FolderName=Results_Folder, FileName="MAC_Distribution", Extension="csv")
        FN.Dictionary_ToCSV(Address, Save_MAC_Distribution)

        Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="MSC_Distribution", Extension="jpg")
        FN.Fig_Plot_Save_1Lines_X_Log_Y_Dictionary_Linear(Address=Address, Identifier=Save_Situation, X_Array=Save_Diameter, Y_array=Save_MSC_Distribution, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1, Y_label1=Y_LabelB4, Plot_Title="MSC Distribution")
        Address = FN.File_Pointer(Main=script_dir, FolderName=Results_Folder, FileName="MSC_Distribution", Extension="csv")
        FN.Dictionary_ToCSV(Address, Save_MSC_Distribution)

        Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="Mass_Sample_Rho_Cte", Extension="jpg")
        FN.Fig_Plot_Save_1Lines_X_Log_Y_Dictionary_Linear(Address=Address, Identifier=Save_Situation, X_Array=Save_Diameter, Y_array=Save_Mass_Sample_2, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_LabelA1, Y_label1=Y_LabelA7, Plot_Title="Mass Distribution using Constant Density")
        Address = FN.File_Pointer(Main=script_dir, FolderName=Results_Folder, FileName="Mass_Sample_Rho_Cte", Extension="csv")
        FN.Dictionary_ToCSV(Address, Save_Mass_Sample_2)

        Address = FN.File_Pointer(Main=script_dir, FolderName=Results_Folder, FileName="MAC", Extension="csv")
        FN.Dictionary_ToCSV(Address, Save_MAC)
        Address = FN.File_Pointer(Main=script_dir, FolderName=Results_Folder, FileName="MSC_Total", Extension="csv")
        FN.Dictionary_ToCSV(Address, Save_MSC_Total)
        Address = FN.File_Pointer(Main=script_dir, FolderName=Results_Folder, FileName="MSC_Diff", Extension="csv")
        FN.Dictionary_ToCSV(Address, Save_MSC_Diff)


    except Exception as e:
        logging.exception(e)
        raise


Sample_Total_Number_Concentration = 5000  # Sample total concentration #/cm^3
Sample_LogN_D_Median_Min = 150  # Smallest median diameter of computation (nm)
Sample_LogN_D_Median_Max = 350  # Largest median diameter of computation (nm)
Sample_LogN_D_Median_Bins = 1  # Number of the Steps

Sample_LogN_Sigma_Min = 1.3  # Smallest Sample Sigma G
Sample_LogN_Sigma_Max = 1.9  # Largest Sample Sigma G
Sample_LogN_Sigma_Bins = 2  # Number of the Steps

####################### Aggregate Distribution bounds

Sample_Sigma_Bound = 3  # Number of Sigma G to cover
Sample_Sigma_Bins = 124  # Number of bins

####################### Effective Density

Eff_dm_Min = 2.56  # Smallest mass mobility exponent for the result to be in kg/m^3
Eff_dm_Max = 2.8  # Largest mass mobility exponent for the result to be in kg/m^3
Eff_dm_Bins = 1  # Number of bins

####################### Effective Density for 100nm aggregate

Eff_rho_100nm_Min = 502  # Smallest effective Density of 100nm aggregate in kg/m^3
Eff_rho_100nm_Max = 650  # Largest effective Density of 100nm aggregate in kg/m^3
Eff_rho_100nm_Bins = 1  # Number of bins

Soot_Material_Density = 1800  # in kg/m^3

####################### Prefactor and projected area exponent of the aggregate

Primary_D_Alpha = 1.1
Primary_K_Alpha = 1.13

#######################

Primary_Sigma_da_CTE_Min = 1  # Smallest Sigma G around specific aggregate mobility diameter
Primary_Sigma_da_CTE_Max = 1.4  # Largest Sigma G around specific aggregate mobility diameter
Primary_Sigma_da_CTE_Bins = 2  # Number of bins

####################### Primary Particle size distrubution options

Primary_Sigma_da_CTE_Bound = 3  # Number of Sigma G to cover
Primary_Sigma_da_CTE_Nt = 75 + 1  # Number of bins

####################### Obsolete Varibales

# Soot_Primary_Agg_dp_sigma = 1  # Geometric std of primary particles within aggregate
# Primary_Betha = 0.9
# Soot_Fractal_D_mc = 0.52  # Fractal Properties for continuum regime
# Soot_Prefactor_k_mc = 0.85 + (0.03 * Soot_Primary_Agg_dp_sigma ** (4.4))  # Fractal continuum regime
# Soot_Fractal_D_mc_RDG = 1 / Soot_Fractal_D_mc
# Soot_Prefactor_k_mc_RDG = (Primary_Betha / Soot_Prefactor_k_mc) ** (1 / Soot_Fractal_D_mc)

####################### RDG Variables

Soot_Fractal_D_mc_RDG = 1.78
Soot_Prefactor_k_mc_RDG = 1.3

####################### Refractive index and Wave Parameters

Soot_Refractive_Index = 2 - 1j  # Complex refractive index
Wave_Length = 870 * 1e-9  # nm
Wave_Number = 2 * pi / Wave_Length  # k
Soot_Complex = (((Soot_Refractive_Index ** 2) - 1) / ((Soot_Refractive_Index ** 2) + 2))
Soot_FM = (abs(Soot_Complex)) ** 2
Soot_EM = Soot_Complex.imag

####################### Scattering Parameters

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
# alpha=2*pi*radius/lambda

####################### Detail Figures

Figure_Enable = 1  # 0 for disable and 1 for enable

####################### List Generation

Sample_Sigma_List = np.linspace(Sample_LogN_Sigma_Min, Sample_LogN_Sigma_Max, num=Sample_LogN_Sigma_Bins)
Sample_D_Median_List = np.linspace(Sample_LogN_D_Median_Min, Sample_LogN_D_Median_Max, num=Sample_LogN_D_Median_Bins)
Eff_dm_List = np.linspace(Eff_dm_Min, Eff_dm_Max, num=Eff_dm_Bins)
Eff_rho_100nm_List = np.linspace(Eff_rho_100nm_Min, Eff_rho_100nm_Max, num=Eff_rho_100nm_Bins)
Primary_Sigma_da_CTE_List = np.linspace(Primary_Sigma_da_CTE_Min, Primary_Sigma_da_CTE_Max, num=Primary_Sigma_da_CTE_Bins)
"""
