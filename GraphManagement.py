from ConfigParserM import logging
import GeneralFunctions as GF
import numpy as np
from matplotlib import rcParams
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.ticker import FormatStrFormatter

####### Plotting Parameters
rcParams['mathtext.fontset'] = 'stix'
rcParams['font.family'] = 'STIXGeneral'


##############

class GraphTools:
    def __init__(self, FolderInfo):
        try:
            self.__CRSGraphs = True

            self.__DmCTE = True
            self.__rhoEff100nmCTE = True
            self.__sigmaMobCTE = True

            self.__folderNameGraph = FolderInfo['FOLDER_NAME_GRAPH']
            self.__folderNameData = FolderInfo['FOLDER_NAME_DATA']
            self.__xLabelPlotCrossSections = "Mobility-equivalent Diameter (nm)"
            self.__yLabelPlotCrossSections = "Cross Section (µm" + "$^{}$".format(2) + ")"

            self.__plotTitle_Size = 12
            self.__labelFont_size = 12
            self.__figureDPI = 700
            self.__markerSize = 3.5
            self.__alphaMainLine = 0.75
            self.__lineColor = ['red', 'blue', 'green']
            self.__lineStyle = ['-', '--', ':']
            self.__markerStyle = ["o", "X", "^"]
            self.__lineWidth = [1.5, 1, 0.5]
            self.__mobilityDiamLimits = [75, 840]

        except Exception as e:
            logging.exception(e)
            raise

    def RDG_TMatrixComparisonGraphs(self):
        try:
            self.MainDF = pd.read_csv(GF.getAddressTo(FolderName=self.__folderNameData, FileName='beacon', Extension='csv'))
            # fileNames = MainDF.loc[MainDF['AA_Plot'] == 1, 'AA_FileName']
            # self.fullDataInfoDF = MainDF.loc[MainDF['AA_Plot'] == 1, 'AA_Plot':]
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
            a = 3
        except Exception as e:
            logging.exception(e)
            raise

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
                Address = GF.getAddressTo(FolderName=self.__folderNameGraph + "\Cross Sections", FileName=F1, Extension="jpg")
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
            tlName = "RDG Absorption Cross Section"
            self.PlotCrossSectionsIntermediate(clName=clName, tlName=tlName)
            ######################## Scattering RDG
            clName = 'SCA_RDG'
            tlName = "RDG Scattering Cross Section"
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
