import os

import matplotlib.pyplot as plt
import pandas as pd

if __name__ == '__main__':
    os.chdir('/Users/u2176312/OneDrive - University of Warwick/Thesis/Paper Draft/')
    Substrates = pd.read_excel('/Users/u2176312/OneDrive - University of Warwick/Thesis/Paper '
                               'Draft/organicPbymolybdatePsip.xlsx', sheet_name='ToPlot2', header=0, index_col=0)
    x_value = list(Substrates.index)
    plt.figure(figsize=(10, 5))
    plt.rcParams["axes.labelweight"] = "bold"

    shapedict = {
        'PE': ['-', 'o'],
        'PC': ['--', 'o'],
        'G3P': ['-.', 'o'],
        'G1P': [':', 'o'],
        'PNPP': ['--', '*']
    }
    for column in Substrates.columns:
        if '_std' in str(column):
            continue
        else:
            y = list(Substrates[column])
            error = list(Substrates[column + '_std'])
            plt.errorbar(x_value, y, yerr=error, capsize=2, capthick=1,
                         fmt=str(shapedict[str(column)][0]) + str(shapedict[str(column)][1]),
                         markerfacecolor='none',
                         color='black',
                         markeredgecolor='black',
                         label=str(column))
    plt.xlabel('Concentration (μM)')
    plt.ylabel('pmols/h of Pi released')
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig('Substrate_paper_Psip1_Andrew_2.pdf', dpi=300)
    plt.show()
