#!/usr/bin/env python

import os

import matplotlib.pyplot as plt
import pandas as pd

if __name__ == '__main__':
    os.chdir('/Users/u2176312/OneDrive - University of Warwick/Thesis/Paper Draft/')
    TableMetals = pd.read_excel('/Users/u2176312/OneDrive - University of Warwick/Thesis/Paper '
                                'Draft/Andrew_PsipMetals.xlsx', sheet_name='Sheet2', header=0, index_col=0)

    x_value = list(TableMetals.index)
    fig, ax = plt.subplots()
    #plt.rcParams["axes.labelweight"] = "bold"
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
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
            ax.errorbar(x_value, y, yerr=error, capsize=2, capthick=1, fmt='-o', markerfacecolor='none',
                        color=colordict[str(column).replace('Calcium + ', '')],
                        markeredgecolor=colordict[str(column).replace('Calcium + ', '')],
                        label=str(column))
        else:
            y = list(TableMetals[column])
            error = list(TableMetals[column + '_std'])
            ax.errorbar(x_value, y, yerr=error, capsize=2, capthick=1, fmt='-o',
                        color=colordict[str(column).replace('Calcium + ', '')],
                        label=str(column))
    plt.xlabel('Time (mins)')
    plt.ylabel('Abs(405nm)')
    plt.legend(loc='upper left')
    plt.tight_layout()
    plt.savefig('Metals_paper_Psip1_Andrew.png', dpi=300)
    plt.show()

    x_value = list(TableMetals.index)
    fig, ax = plt.subplots()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #plt.rcParams["axes.labelweight"] = "bold"
    shapedict = {
        'Calcium': 's',
        'Iron': 'o',
        'Cobalt': 'X',
        'Magnesium': 'p',
        'Manganese': 'd'
    }
    for column in TableMetals.columns:
        if '_std' in str(column):
            continue
        if 'Calcium +' in str(column):
            y = list(TableMetals[column])
            error = list(TableMetals[column + '_std'])
            ax.errorbar(x_value, y, yerr=error, capsize=2, capthick=1,
                        fmt='--' + str(shapedict[str(column).replace('Calcium + ', '')]),
                        markerfacecolor='none',
                        color='black',
                        markeredgecolor='black',
                        label=str(column))
        else:
            y = list(TableMetals[column])
            error = list(TableMetals[column + '_std'])
            ax.errorbar(x_value, y, yerr=error, capsize=2, capthick=1,
                        fmt='-' + str(shapedict[str(column).replace('Calcium + ', '')]),
                        markerfacecolor='none',
                        color='black',
                        markeredgecolor='black',
                        label=str(column))
    plt.xlabel('Time (mins)')
    plt.ylabel('Abs(405nm)')
    # plt.legend(loc='upper left')
    plt.tight_layout()
    plt.savefig('Metals_paper_Psip1_Andrew_2.png', dpi=300)
    plt.show()

    os.chdir('/Users/u2176312/OneDrive - University of Warwick/Thesis/Paper Draft/')
    TableMetals = pd.read_excel('/Users/u2176312/OneDrive - University of Warwick/Thesis/Paper '
                                'Draft/Andrew_PsipMetals.xlsx', sheet_name='Sheet3', header=0, index_col=0)
    x_value = list(TableMetals.index)
    fig, ax = plt.subplots(dpi=450)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #plt.rcParams["axes.labelweight"] = "bold"
    shapedict = {
        'Calcium': 's',
        'Iron': 'o'
    }
    for column in TableMetals.columns:
        if '_std' in str(column):
            continue
        if 'Calcium' in str(column):
            y = list(TableMetals[column])
            error = list(TableMetals[column + '_std'])
            ax.errorbar(x_value, y, yerr=error, capsize=2, capthick=1,
                        fmt='--' + str(shapedict[str(column).replace('Calcium + ', '')]),
                        markerfacecolor='black',
                        color='black',
                        markeredgecolor='black',
                        label=str(column))
        else:
            y = list(TableMetals[column])
            error = list(TableMetals[column + '_std'])
            ax.errorbar(x_value, y, yerr=error, capsize=2, capthick=1,
                        fmt='-' + str(shapedict[str(column).replace('Calcium + ', '')]),
                        markerfacecolor='none',
                        color='black',
                        markeredgecolor='black',
                        label=str(column))
    plt.xlabel('Time (mins)')
    plt.ylabel('Abs(405nm)')
    ax.legend(loc='upper left', fontsize=15)
    plt.tight_layout()
    plt.savefig('Metals_paper_Psip1_Andrew_3.pdf', dpi=450)
    plt.show()
