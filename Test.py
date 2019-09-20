import os
from pathlib import Path
from ConfigParserM import logging
import subprocess
from matplotlib import pyplot as plt


class SVG_WMF_Plotting:

    def __init__(self):
        self.__folderNameGraph = 'Graphs'
        self.__WMF_SVGSaving = True
        self.__inkScapePath = "C://Program Files//inkscape//inkscape.exe"
        self.__figureDPI = 500

    def getRootDirectory(self):
        try:
            return Path(os.path.dirname(os.path.realpath('__file__')))

        except Exception as e:
            logging.exception(e)
            raise

    def getAddressTo(self, Main=None, FolderName=None, FileName=None, Extension=None):
        try:
            if Main is None:
                Main = self.getRootDirectory()
            if FolderName:
                Path1 = Path(Main) / Path(FolderName)
            else:
                Path1 = Path(Main)

            if not os.path.exists(Path1):
                os.makedirs(Path1)
            if FileName:
                if Extension:
                    File_Address = Path1 / Path(FileName + "." + Extension)
                else:
                    File_Address = Path1 / Path(FileName)
            else:
                File_Address = Path1
            return File_Address

        except Exception as e:
            logging.exception(e)
            raise

    def TestPlot(self):
        try:

            fig, ax1 = plt.subplots()
            x = [1, 2]
            y = [1, 2]
            F1 = 'test'
            ax1.plot(x, y)
            self.SaveAndClosePlot(folderName=self.__folderNameGraph, fileName=F1)


        except Exception as e:
            logging.exception(e)
            raise

    def SaveAndClosePlot(self, folderName, fileName):
        try:
            Address = self.getAddressTo(FolderName=self.__folderNameGraph + f"\{folderName}", FileName=fileName, Extension="jpg")
            plt.savefig(Address, format='jpg', dpi=self.__figureDPI, bbox_inches='tight')

            if self.__WMF_SVGSaving:
                Address = self.getAddressTo(FolderName=self.__folderNameGraph + f"\{folderName}", FileName=fileName, Extension="svg")
                plt.savefig(Address, format='svg', dpi=self.__figureDPI, bbox_inches='tight')
                # add removing SVG if needed

                AddressWMF = self.getAddressTo(FolderName=self.__folderNameGraph + f"\{folderName}", FileName=fileName, Extension="wmf")
                subprocess.call([self.__inkScapePath, str(Address.resolve()), '--export-wmf', str(AddressWMF.resolve())])

            plt.clf()
            plt.close()
        except Exception as e:
            logging.exception(e)
            raise
