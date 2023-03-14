import os
import statistics

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import optimize
from sklearn.metrics import r2_score

os.chdir('/Users/u1866168/Documents/OneDrive - University of '
         'Warwick/Experiments/Analysis_Genes/Psip1/Activity_Experiments/')
Generaldf = pd.read_excel('Phosphoethanolamine.xlsx', sheet_name='DataTransformed', index_col=0)

#%% Constant concentration of PNPP

PNPP_80uM = Generaldf[Generaldf.filter(like=str('200uM_PNPP')).columns]

index = np.array(PNPP_80uM.index)  # Generate the X values (x-axis)
amounts = [40, 20, 10, 5]

slopes = {}
r_squares = {}
for row in PNPP_80uM:  # This will fit a curve using the x-axis for each column
    y = np.array(PNPP_80uM[row])
    idx = np.isfinite(index) & np.isfinite(y)  # Test to prevent NaNs
    ab = np.polyfit(index[idx], y[idx], 1)  # Combine those elements with no NaNs
    predict = np.poly1d(ab)
    r_score = r2_score(y[idx], predict(index[idx]))
    if str(row).split("_")[2].replace('uM', '') in slopes:
        slopes[str(row).split("_")[2].replace('uM', '')].update({str(row).rsplit("_")[3]: ab[0]})
        r_squares[str(row).split("_")[2].replace('uM', '')].update({str(row).rsplit("_")[3]: r_score})
    else:
        slopes.update({str(row).split("_")[2].replace('uM', ''): {str(row).rsplit("_")[3]: ab[0]}})
        r_squares.update({str(row).split("_")[2].replace('uM', ''): {str(row).rsplit("_")[3]: r_score}})

velocity = np.array([list(slopes[str(x)].values()) for x in amounts])
means = [statistics.mean(x) for x in velocity]
error = [statistics.stdev(x) for x in velocity]


def michaelis_menten(x, Vm, Km):
    return Vm * (x / (Km + x))


def competitive_inhibitor_michaelis_menten(x2, Vm, Km2):
    return (Vm * x2) / (2.5 * (1 + (40 / Km2)) + 10)


p0 = [0.0001, 0.0001]
params, cv = optimize.curve_fit(competitive_inhibitor_michaelis_menten, amounts, means, absolute_sigma=True, p0=p0)
vmax, km = params

plt.rcParams["axes.labelweight"] = "bold"
plt.errorbar(amounts, means, yerr=error, fmt='o-', capsize=2, capthick=1)

plt.show()

#%% Constant concentration of inhibitor


PE = Generaldf[Generaldf.filter(like=str('40uM_PE')).columns]

index = np.array(PE.index)  # Generate the X values (x-axis)
amounts = [80, 40, 20, 10,  5]

slopes = {}
r_squares = {}
for row in PE:  # This will fit a curve using the x-axis for each column
    y = np.array(PE[row])
    idx = np.isfinite(index) & np.isfinite(y)  # Test to prevent NaNs
    ab = np.polyfit(index[idx], y[idx], 1)  # Combine those elements with no NaNs
    predict = np.poly1d(ab)
    r_score = r2_score(y[idx], predict(index[idx]))
    if str(row).split("_")[0].replace('uM', '') in slopes:
        slopes[str(row).split("_")[0].replace('uM', '')].update({str(row).rsplit("_")[3]: ab[0]})
        r_squares[str(row).split("_")[0].replace('uM', '')].update({str(row).rsplit("_")[3]: r_score})
    else:
        slopes.update({str(row).split("_")[0].replace('uM', ''): {str(row).rsplit("_")[3]: ab[0]}})
        r_squares.update({str(row).split("_")[0].replace('uM', ''): {str(row).rsplit("_")[3]: r_score}})

velocity = np.array([list(slopes[str(x)].values()) for x in amounts])
means = [statistics.mean(x) for x in velocity]
error = [statistics.stdev(x) for x in velocity]

plt.rcParams["axes.labelweight"] = "bold"
plt.errorbar(amounts, means, yerr=error, fmt='o-', capsize=2, capthick=1)

plt.show()

p0 = [0.0001, 0.0001]
params, cv = optimize.curve_fit(michaelis_menten, amounts, means, absolute_sigma=True, p0=p0)
vmax, km = params
