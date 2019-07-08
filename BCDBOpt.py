import subprocess
import blackbox as bb
import sys
import configparser

InputArray = sys.argv[1:]

Df_Min = float(InputArray[0])
Df_Max = float(InputArray[1])
GlobalSteps = int(InputArray[2])
LocalSteps = int(InputArray[3])
Batches = int(InputArray[4])
Number = int(InputArray[5])

config = configparser.ConfigParser()
config.optionxform = str

config.read('config/BC_DB_Execute.cnf')
configFull = {}
for section in config.sections():
    configFull[section] = {}
    for key, val in config.items(section):
        configFull[section][key] = val

BC_DB_Limits = configFull['BC_DB_Limits']

Real_RI_Max = float(BC_DB_Limits['RI_Real_Max'])
Real_RI_Min = float(BC_DB_Limits['RI_Real_Min'])
Imag_RI_Max = float(BC_DB_Limits['RI_Imag_Max'])
Imag_RI_Min = float(BC_DB_Limits['RI_Imag_Min'])
Np_Max = float(BC_DB_Limits['Np_Max'])
Np_Min = float(BC_DB_Limits['Np_Min'])
dp_Max = float(BC_DB_Limits['dp_Max'])
dp_Min = float(BC_DB_Limits['dp_Min'])
Wavelength_Max = float(BC_DB_Limits['Wavelength_Max'])
Wavelength_Min = float(BC_DB_Limits['Wavelength_Min'])


def BC_opt_CoreV1(ParametersList):
    command = "python BCDBExecuteC" + str(Number) + ".py " + str(ParametersList[0]) + " " + str(ParametersList[1]) + " " + str(ParametersList[2]) + " " + str(ParametersList[3]) + " " + str(ParametersList[4]) + " " + str(ParametersList[5])
    out = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding='utf-8')
    A = out.stdout
    return float(A)


def OptimizeBC(Constraint, GlobalSteps, LocalSteps, Batches):
    bb.search(f=BC_opt_CoreV1,  # given function
              box=Constraint,  # [[-10., 10.], [-10., 10.]],  # range of values for each parameter (2D case)
              n=GlobalSteps,  # 20,  # number of function calls on initial stage (global search)
              m=LocalSteps,  # 20,  # number of function calls on subsequent stage (local search)
              batch=Batches,  # 4,  # number of calls that will be evaluated in parallel
              resfile="Optimization/BC_Opt_C" + str(Number) + ".csv")  # text file where results will be saved


if __name__ == '__main__':
    OptimizeBC(Constraint=[[Df_Min, Df_Max], [Real_RI_Min, Real_RI_Max], [Imag_RI_Min, Imag_RI_Max], [Np_Min, Np_Max], [dp_Min, dp_Max], [Wavelength_Min, Wavelength_Max]], GlobalSteps=GlobalSteps, LocalSteps=LocalSteps, Batches=Batches)
