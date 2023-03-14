#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import snakemake

if __name__ == '__main__':
    TableMetals = pd.read_excel(str(snakemake.input[0]), sheet_name='Sheet2', header=0, index_col=0)

    x_value = list(TableMetals.index)
    plt.figure()
    plt.rcParams["axes.labelweight"] = "bold"
    colordict = {
        'Calcium': 'blue',
        'Iron': 'red',
        'Cobalt': 'green',
        'Magnesium': 'purple',
        'Manganese': 'orange'
    }

    for column in TableMetals.columns:
        if '_std' in str(column):
            continue
        if 'Calcium +' in str(column):
            y = list(TableMetals[column])
            error = list(TableMetals[column + '_std'])
            plt.errorbar(x_value, y, yerr=error, capsize=2, capthick=1, fmt='-o', markerfacecolor='none',
                         color=colordict[str(column).replace('Calcium + ', '')],
                         markeredgecolor=colordict[str(column).replace('Calcium + ', '')],
                         label=str(column))
        else:
            y = list(TableMetals[column])
            error = list(TableMetals[column + '_std'])
            plt.errorbar(x_value, y, yerr=error, capsize=2, capthick=1, fmt='-o',
                         color=colordict[str(column).replace('Calcium + ', '')],
                         label=str(column))
    plt.xlabel('Time (mins)')
    plt.ylabel('Abs(405nm)')
    plt.legend(loc='upper left')
    plt.tight_layout()
    plt.savefig(str(snakemake.output[0]), dpi=300)
    plt.show()