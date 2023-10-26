import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

if __name__ == '__main__':
    os.chdir('/Users/u2176312/OneDrive - University of Warwick/Thesis/Paper Draft/')
    Substrates = pd.read_excel('/Users/u2176312/OneDrive - University of Warwick/Thesis/Paper '
                               'Draft/organicPbymolybdatePsip.xlsx', sheet_name='ToPlot2', header=0, index_col=0)
    x_value = list(Substrates.index)
    mpl.rcParams['font.size'] = 15
    mpl.rcParams['font.family'] = 'Arial'
    fig, ax = plt.subplots(figsize=(11, 11), dpi=450)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    shapedict = {
        'PE': ['-', 'D', 'none'],
        'PC': ['--', 'D', 'black'],
        'G3P': ['-.', 'o', 'none'],
        'G1P': [':', 'o', 'black'],
        'PNPP': ['--', 's', 'black']
    }
    for column in Substrates.columns:
        if '_std' in str(column):
            continue
        else:
            y = list(Substrates[column])
            error = list(Substrates[column + '_std'])
            plt.errorbar(x_value, y, yerr=error, capsize=2, capthick=1,
                         fmt=str(shapedict[str(column)][0]) + str(shapedict[str(column)][1]),
                         # markerfacecolor=shapedict[str(column)][2],
                         # color='black',
                         # markeredgecolor='black',
                         label=str(column),
                         markersize=8)
    plt.xlabel('Concentration (μM)')
    plt.ylabel('pmol Pi released h$^{-1}$')
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig('Substrate_paper_Psip1_Andrew_coloured.pdf', dpi=600)
    plt.show()
