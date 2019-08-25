'''
TMatrix_Table_Input_Headers1 = ['Df', 'kf', 'R_RI', 'I_RI', 'WaveL', 'dp', 'Np', 'Version']
TMatrix_Table_Output_Headers1 = ['TMatrix_ABS_CRS', 'TMatrix_SCA_CRS', 'MLink']
Array1 = [[2.3, 1.2, 1.6, 0.6, 860, 33, 50, 0.1], [2.4, 1.2, 1.6, 0.6, 860, 33, 50, 0.1], [2.6, 1.2, 1.6, 0.6, 860, 33, 50, 0.1]]
Array2 = [[0.056, 0.666, 20], [0.053, 0.666, 21], [0.051, 0.666, 22]]
DB.insertArrayIntoTable(INFO=DB_Info, TableName=TMatrix_Table_Name, NameArray=TMatrix_Table_Input_Headers1, Array=Array1)
DB.insertArrayIntoTable(INFO=DB_Info, TableName=TMatrix_Table_Name + "_out", NameArray=TMatrix_Table_Output_Headers1, Array=Array2)
DB.showAllTablesInDBSummary(DB_Info)
'''
'''
    TMatrix_Main_Input_Array.append([2.3, 1.2, 1.6, 0.6, 860, 33, 50, 0.1])
    TMatrix_Main_Input_Array.append([2.4, 1.2, 1.6, 0.6, 860, 33, 50, 0.1])
    TMatrix_Main_Input_Array.append([2.6, 1.2, 1.6, 0.6, 860, 33, 50, 0.1])
    TMatrix_Main_Input_Array.append([2.1, 1.2, 1.6, 0.6, 860, 33, 50, 0.1])
    '''
import DBManagement
import matplotlib.pyplot as plt

S1 = DBManagement.MySQLManagement()

def corrdot(*args, **kwargs):
    corr_r = args[0].corr(args[1], 'pearson')
    corr_text = f"{corr_r:2.2f}".replace("0.", ".")
    ax = plt.gca()
    ax.set_axis_off()
    marker_size = abs(corr_r) * 10000
    ax.scatter([.5], [.5], marker_size, [corr_r], alpha=0.6, cmap="coolwarm",
               vmin=-1, vmax=1, transform=ax.transAxes)
    font_size = abs(corr_r) * 40 + 5
    ax.annotate(corr_text, [.5, .5, ], xycoords="axes fraction",
                ha='center', va='center', fontsize=font_size)


sns.set(style='white', font_scale=1.6)
iris = sns.load_dataset('iris')
g = sns.PairGrid(iris, aspect=1.4, diag_sharey=False)
g.map_lower(sns.regplot, ci=False, line_kws={'color': 'black'})
g.map_diag(sns.distplot, kde_kws={'color': 'black'})
g.map_upper(corrdot)
plt.show()
