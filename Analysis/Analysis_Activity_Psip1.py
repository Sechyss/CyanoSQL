import os
import statistics

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scikit_posthocs as sp
from scipy import stats, optimize

from CyanoScripts.MyPackage.ActivityAnalysisClass import MetalActivityData, plot_curve_errorbar, \
    plot_abs_time, r_coefficient_values, michaelis_menten, transform_substratedf_productdf
from CyanoScripts.MyPackage.Test_Km import plot_lineweaverburk

os.chdir('/Users/u1866168/Documents/OneDrive - University of '
         'Warwick/Experiments/Analysis_Genes/Psip1/Activity_Experiments/')

# %% Analysis of Bis-pNPP
bpnppdata = pd.read_excel('Psip1_experiment_pag41_20211125.xlsx', sheet_name='Corrected', index_col='Time')
experiment = MetalActivityData(bpnppdata)
curves, rvalues = experiment.fit_curve()
bisPNPP = []
pNpp = []

for i in curves.keys():
    bisPNPP.append(curves[i]['BpNPP'])
    pNpp.append(curves[i]['pNPP'])

abs_mg_min_pNPP = [x / (1.5*1000/2) for x in pNpp]
abs_mg_min_BpNPP = [x / (1.5*1000/2) for x in bisPNPP]

plt.figure(figsize=(8, 7))
plt.rcParams["axes.labelweight"] = "bold"
list_to_plot = ['Mg', 'Mn', 'Fe', 'Ca']

X_axis = np.arange(len(list_to_plot))

plt.bar(X_axis - 0.2, abs_mg_min_BpNPP, 0.4, label='BpNPP', color='g')
plt.bar(X_axis + 0.2, abs_mg_min_pNPP, 0.4, label='pNPP')

plt.xticks(X_axis, list_to_plot)
plt.xlabel("Metals (0.1 mM)")
plt.ylabel("ΔAbs(405nm)/min/mg of protein")
plt.legend()
plt.show()

plt.show()

# %% Plotting metals abs and mins

trial = pd.read_excel('Psip1_experiment_pag47_table.xlsx', sheet_name='MeanReplicatesCorr', index_col=0)

plt.figure(figsize=(8, 7))

plot_abs_time('Ca', trial)
# %% Loading the data for metals

plate1 = pd.read_excel('Psip1_experiment_pag47_table.xlsx', sheet_name='Plate1_mancorrected', index_col=0)
plate2 = pd.read_excel('Psip1_experiment_pag47_table.xlsx', sheet_name='Plate2_mancorrected', index_col=0)
experiment1 = MetalActivityData(plate1)
experiment2 = MetalActivityData(plate2)
curves1, rvalues1 = experiment1.fit_curve()
curves2, rvalues2 = experiment2.fit_curve()

# %% Creation of elements to plot
fig = plt.figure(figsize=(8, 7))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

metal1 = np.array(list(curves1['Zn'].values()))
r_scoreZn = np.array(list(rvalues1['Zn'].values()))
metal2 = np.array(list(curves2['Zn'].values()))
r_scoreZn2 = np.array(list(rvalues2['Zn'].values()))
meanmetal = []
varmetal = []

# %% Once you have removed the outliers run the plotting
r_coefficient_values(metal1, metal2, meanmetal, varmetal)
x = list(curves1['Zn'].keys())
xint = [float(item) for item in x]

ax1.errorbar(xint, meanmetal, yerr=varmetal, fmt='-o', color='c', capsize=2, capthick=1)

metal1 = np.array(list(curves1['Fe'].values()))
r_scoreFe = np.array(list(rvalues1['Fe'].values()))
metal2 = np.array(list(curves2['Fe'].values()))
r_scoreFe2 = np.array(list(rvalues2['Fe'].values()))
meanmetal = []
varmetal = []

# %% Check point to remove the Fe
r_coefficient_values(metal1, metal2, meanmetal, varmetal)
x = list(curves1['Fe'].keys())
xint = [float(item) for item in x]

ax1.errorbar(xint, meanmetal, yerr=varmetal, fmt='-o', color='r', capsize=2, capthick=1)

ax1.legend_2(['Zn', 'Fe'], loc='lower left', prop={'size': 18})

metal1 = np.array(list(curves1['Mg'].values()))
r_scoreMg = np.array(list(rvalues1['Mg'].values()))
metal2 = np.array(list(curves2['Mg'].values()))
r_scoreMg2 = np.array(list(rvalues2['Mg'].values()))
meanmetal = []
varmetal = []

metal1[4] = metal2[4]  # This will remove the outlier at that particular time point

# %% Check point to remove the Mg
r_coefficient_values(metal1, metal2, meanmetal, varmetal)
x = list(curves1['Mg'].keys())
xint = [float(item) for item in x]

ax2.errorbar(xint, meanmetal, yerr=varmetal, fmt='-o', color='y', capsize=2, capthick=1)

metal1 = np.array(list(curves1['Mn'].values()))
r_scoreMn = np.array(list(rvalues1['Mn'].values()))
metal2 = np.array(list(curves2['Mn'].values()))
r_scoreMn2 = np.array(list(rvalues2['Mn'].values()))
meanmetal = []
varmetal = []

metal2[4] = metal1[4]
metal2[3] = metal1[3]

# %% Check point to remove the Mn
r_coefficient_values(metal1, metal2, meanmetal, varmetal)
x = list(curves1['Mn'].keys())
xint = [float(item) for item in x]

ax2.errorbar(xint, meanmetal, yerr=varmetal, fmt='-o', color='m', capsize=2, capthick=1)

metal1 = np.array(list(curves1['Ca'].values()))
r_scoreCa = np.array(list(rvalues1['Ca'].values()))
metal2 = np.array(list(curves2['Ca'].values()))
r_scoreCa2 = np.array(list(rvalues2['Ca'].values()))
meanmetal = []
varmetal = []

metal2[1] = metal1[1]
metal2[2] = metal1[2]

# %% Check point to remove the Ca
r_coefficient_values(metal1, metal2, meanmetal, varmetal)
x = list(curves1['Ca'].keys())
xint = [float(item) for item in x]

ax2.errorbar(xint, meanmetal, yerr=varmetal, fmt='-o', color='b', capsize=2, capthick=1)

ax2.legend_2(['Mg', 'Mn', 'Ca'], loc='upper left', prop={'size': 18})

ax1.set_xlabel('Concentration (mM)')
ax1.set_ylabel('R coefficient')

plt.savefig('./Figures_22/R_coefficients_metals.pdf')
plt.show()

# %% Analysis of pH (first draft)

plate = pd.read_excel('Psip1_experiment_pag52_table.xlsx', sheet_name='Corrected', index_col=0)
experiment = MetalActivityData(plate)
curves = experiment.fit_curve()

pH = np.array([list(curves['pH6.8'].values()), list(curves['pH7.5'].values()), list(curves['pH7.8'].values()),
               list(curves['pH8.0'].values()), list(curves['pH8.8'].values())])
means = [statistics.mean(x) for x in pH]
error = [statistics.stdev(x) for x in pH]
x = [6.8, 7.5, 7.8, 8.0, 8.8]

fig = plt.figure(figsize=(8, 7))

plt.errorbar(x, means, yerr=error, fmt='-o', color='b', capsize=2, capthick=1)
plt.xlabel('pH values')
plt.ylabel('R coefficient')

plt.show()

# %% Analysis of pH final

trial = pd.read_excel('Psip1_experiment_pag53_table.xlsx', sheet_name='Study_noOut', index_col=0)
experiment = MetalActivityData(trial)
curves, rvalue = experiment.fit_curve()

pH = list([list(curves['pH6.8'].values()), list(curves['pH7.5'].values()), list(curves['pH8.8'].values()),
           list(curves['pH9.4'].values()), list(curves['pH9.8'].values()), list(curves['pH10.4'].values()),
           list(curves['pH11.2'].values())])

r_score = list([list(rvalue['pH6.8'].values()), list(rvalue['pH7.5'].values()), list(rvalue['pH8.8'].values()),
                list(rvalue['pH9.4'].values()), list(rvalue['pH9.8'].values()), list(rvalue['pH10.4'].values()),
                list(rvalue['pH11.2'].values())])
means = [statistics.mean(x) for x in pH]
error = [statistics.stdev(x) for x in pH]
x = [6.8, 7.5, 8.8, 9.4, 9.8, 10.4, 11.2]

plt.figure()
plt.rcParams["axes.labelweight"] = "bold"
plt.errorbar(x, means, yerr=error, fmt='-o', capsize=2, capthick=1, color='black')
plt.xlabel('pH values')
plt.ylabel('ΔAbs(405nm)/min/mg of protein')
plt.tight_layout()
plt.savefig('./Figures_22/R_coefficient_pH_paper.pdf', dpi = 300)
plt.show()

# %% Final metal comparison

metals = pd.read_excel('Psip1_experiment_pag56_table.xlsx', sheet_name='Raw_NoControl', index_col=0)
metals_NoCa = pd.read_excel('Psip1_experiment_pag56_table.xlsx', sheet_name='Raw_Controls', index_col=0)
plt.figure(figsize=(8, 7))
plt.rcParams["axes.labelweight"] = "bold"
list_to_plot = ['Mg', 'Mn', 'Zn', 'Fe', 'Ca']

for element in list_to_plot:
    plot_curve_errorbar(metals, element)

plt.xlabel('Time(mins)')
plt.ylabel('Abs 405nm')
plt.legend(list_to_plot, loc='upper left', prop={'size': 18})
#plt.savefig('./Figures_22/Raw_data_metals01mM.pdf')
plt.show()

experiment = MetalActivityData(metals)
curves, rvalue_fit = experiment.fit_curve()

experiment_noCa = MetalActivityData(metals_NoCa)
curves_noCa, rvalue_fit_noCa = experiment_noCa.fit_curve()

metal = np.array([list(curves[x].values()) for x in list_to_plot])

means = [statistics.mean(x) for x in metal]
error = [statistics.stdev(x) for x in metal]

noCavalues = [list(curves_noCa[x].values()) for x in list_to_plot]
y_noMetal = [x[0] for x in noCavalues]

x_pos = np.arange(len(list_to_plot))

figure, ax = plt.subplots()
plt.rcParams["axes.labelweight"] = "bold"
ax.bar(x_pos, means, yerr=error, align='center', ecolor='black', capsize=10, edgecolor='black', color='white')
ax.bar(x_pos, y_noMetal, align='center', color='black')
ax.set_ylabel('ΔAbs(405nm)/min/mg of protein')
ax.set_xticks(x_pos, weight='bold')
ax.set_xticklabels(list_to_plot)
ax.yaxis.grid(False)

# Save the figure and show
plt.tight_layout()
plt.savefig('./Figures_22/R_coefficient_metals01mM_paper.pdf', dpi=300)
plt.show()

# %% Statistics to see differences between groups

Mg = list(curves['Mg'].values())
Mn = list(curves['Mn'].values())
Zn = list(curves['Zn'].values())
Fe = list(curves['Fe'].values())
Ca = list(curves['Ca'].values())

print(stats.bartlett(Mg, Mn, Zn, Fe, Ca))  # To test their variance is the equal, this step is necessary before
# Kruskal-Wallis, test of Levenne cannot be used since data is not normally distributed
print(stats.kruskal(Mg, Mn, Zn, Fe, Ca))

sp.posthoc_dunn([Mg, Mn, Zn, Fe, Ca], p_adjust='bonferroni')

df = pd.DataFrame({'data': [0.0005044600401535704, 0.0006028494945581356, 0.0005855993096403793, 0.0008213606424571163,
                            0.0009164523968863366, 0.0007230680497340704, 0.00011314765242506491, 6.336074812440481e-05,
                            0.0010105667288929589, 0.0010374678595329496, 0.0010340116586242125, 0.000732181853404248,
                            0.0006763543024197808, 0.0007022207741890035],
                   'metal': ['Mg', 'Mg', 'Mg', 'Mn', 'Mn', 'Mn', 'Zn', 'Zn', 'Fe', 'Fe', 'Fe', 'Ca', 'Ca', 'Ca']})

sp.posthoc_dunn(df, group_col='metal', val_col='data', p_adjust='bonferroni')

# %% Analysis of kinetics (exploration)

subsrate = pd.read_excel('Psip1_experiment_pag61_table.xlsx', sheet_name='Study', index_col=0)
list_to_plot = [2, 2.5, 3, 3.5, 4, 4.5, 5, 10]
legend = ['2 μM', '2.5 μM', '3 μM', '3.5 μM', '4 μM', '4.5 μM', '5 μM', '10 μM']
plt.figure(figsize=(8, 7))
plt.rcParams["axes.labelweight"] = "bold"
for element in list_to_plot:
    plot_curve_errorbar(subsrate, element)

plt.xlabel('Time(mins)')
plt.ylabel('Abs 405nm')
plt.legend(legend, loc='upper left', prop={'size': 15})
#plt.savefig('./Figures_22/Raw_data_kinetics.pdf')
plt.show()

# --- Fit the curves using the transformed abs to nmoles of product

standard = pd.read_excel('Psip1_experiment_pag61_table.xlsx', sheet_name='StandardCurve', index_col=0)
pnp = np.array(standard.index)
abs_y = np.array(standard['Mean_abs'])

plt.errorbar(pnp, abs_y, yerr=[standard['std']], fmt='-', capsize=2, capthick=1)
plt.xlabel('PNP (nmoles)')
plt.ylabel('Abs 405nm')

slope, intercept = np.polyfit(pnp, abs_y, 1)
model = np.poly1d(np.polyfit(pnp, abs_y, 1))

plt.plot(pnp, model(pnp), "r--")
plt.legend(['y = ' + str("{:e}".format(slope)) + 'x + ' + str("{:e}".format(intercept)), 'Abs_PNP'], loc='upper left')
# plt.savefig('./Figures_22/StandardCurve_data_kinetics.pdf')
plt.show()

product_df = pd.DataFrame(index=subsrate.index, columns=list(subsrate))

for row in subsrate:
    product_df[row] = (subsrate[row] - intercept) / slope

product_df = product_df/0.002

experiment = MetalActivityData(product_df)
curves, rvalue_fit = experiment.fit_curve()
velocity = np.array([list(curves[str(x)].values()) for x in list_to_plot])  # mM/min

means = [statistics.mean(x) for x in velocity]
error = [statistics.stdev(x) for x in velocity]


# --- Estimation of Km and Vmax by fitting Michaelis-Menten curve into the data
p0 = [0.0001, 0.0001]
params, cv = optimize.curve_fit(michaelis_menten, list_to_plot, means, absolute_sigma=True, p0=p0)
vmax, km = params


print("Vmax = ", vmax)
print("Km = ", km)
print('Estimated variance (Vm, Km) = '+ str(cv[0,0]) + ', '+ str(cv[1,1]))
print('Estimated standard devitation (Vm, Km) = ', np.sqrt(np.diag(cv)))


plt.rcParams["axes.labelweight"] = "bold"
plt.errorbar(list_to_plot, means, yerr=error, fmt='o', capsize=2, capthick=1)
plt.plot(np.linspace(0, 15, 1000),michaelis_menten(np.linspace(0, 15, 1000), vmax, km),
         "g--", label='Fitted Michaelis-Menten equation')

plt.xlabel('PNPP (μM)')
plt.ylabel('Velocity (nmoles/min/mg)')
plt.tight_layout()
plt.legend(loc='lower right')
plt.savefig('./Figures_22/RateVelocitykinetics_MM_paper.pdf', dpi = 300)
plt.show()


# --- Estimation of Km and Vmax using the Lineweaver-Burk plot
substrate_conc = list_to_plot[2:]
v_values = means[2:]

plot_lineweaverburk(substrate_conc, v_values)
#plt.savefig('./Figures_22/LineweaverBurk_kinetics_MM.pdf', dpi = 300)

plt.show()

# %% Final analysis of kinetics (Combination of both analysis)

# --- Loading of dataframes with data normalized by the negative control
concentrations_1 = pd.read_excel('Psip1_experiment_pag61_table.xlsx', sheet_name='Test', index_col=0)
concentrations_2 = pd.read_excel('Psip1_experiment_pag63_2_table.xlsx', sheet_name='Test', index_col=0)
AdditionalCon = pd.read_excel('Psip1_experiment_pag61_table.xlsx', sheet_name='Test2', index_col=0)

# --- Loading of standard curve tables
standard_1 = pd.read_excel('Psip1_experiment_pag61_table.xlsx', sheet_name='StandardCurve', index_col=0)
standard_2 = pd.read_excel('Psip1_experiment_pag63_2_table.xlsx', sheet_name='StandardCurve', index_col=0)
standard_3 = pd.read_excel('Psip1_experiment_pag61_table.xlsx', sheet_name='StandardCurveTest', index_col=0)

# --- Transformation of df to productdf to fit the curves of Michaelis-Menten
product_df1 = transform_substratedf_productdf(concentrations_1, standard_1, 0.0020)
product_df2 = transform_substratedf_productdf(concentrations_2, standard_2, 0.0020)
AddtionalProductdf = transform_substratedf_productdf(AdditionalCon, standard_3, 0.0028)

Concat_df = pd.merge(product_df1, AddtionalProductdf, left_index=True, right_index=True)

# --- Specify the concentrations used for the assays
amounts_1 = [0, 2, 2.5, 3, 3.5, 4, 5, 10, 15]
amounts_2 = [0, 1.5, 1.8, 2.2, 2.4, 2.5, 3, 4, 5, 6, 10]

# --- Fitting of the curves and combination of data into one single df
experiment = MetalActivityData(Concat_df)
curves, rvalue_fit = experiment.fit_curve()

velocity = np.array([list(curves[str(x)].values()) for x in amounts_1])
means = [statistics.mean(x) for x in velocity]
error = [statistics.stdev(x) for x in velocity]

experiment_2 = MetalActivityData(product_df2)
curves_2, rvalue_fit_2 = experiment_2.fit_curve()

velocity_2 = np.array([list(curves_2[str(x)].values()) for x in amounts_2])
means_2 = [statistics.mean(x) for x in velocity_2]
error_2 = [statistics.stdev(x) for x in velocity_2]

# --- Combination of the points for Michaelis-Menten
for d in curves_2.keys():
    if d not in curves.keys():
        curves.update({d: curves_2[d]})
    else:
        for i in curves[d]:
            curves[d][i] = [curves[d][i],curves_2[d][i]]

x_pos = [0, 1.5, 1.8, 2, 2.2, 2.4, 2.5, 3, 3.5, 4, 5, 6, 10, 15]

velocityMM = []

for i in x_pos:
    collector = []
    for element in curves[str(i)].values():
        if type(element) == list:
            collector.append(element[0]); collector.append(element[1])
        else:
            collector.append(element)
    velocityMM.append(collector)

meansMM = [statistics.mean(x) for x in velocityMM]
errorMM = [statistics.stdev(x) for x in velocityMM]

# --- Fitting of the
p0 = [10, 10]
params, cv = optimize.curve_fit(michaelis_menten, amounts_1, means)
vmax, km = params
print("Vmax = ", vmax)
print("Km = ", km)
print('Estimated variance (Vm, Km) = '+ str(cv[0,0]) + ', '+ str(cv[1,1]))
print('Estimated standard devitation (Vm, Km) = ', np.sqrt(np.diag(cv)))

residuals = means-michaelis_menten(amounts_1, *params)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((means-np.mean(means))**2)
r_squared = 1 - (ss_res / ss_tot)

plt.rcParams["axes.labelweight"] = "bold"
plt.errorbar(amounts_1, means, yerr=error, fmt='o', capsize=2, capthick=1, color='black')
plt.errorbar(amounts_2, means_2, yerr=error_2, fmt='ko', capsize=2, capthick=1, markerfacecolor='none',
             markeredgecolor='black')
plt.plot(np.linspace(0, 15, 1000),michaelis_menten(np.linspace(0, 15, 1000), vmax, km),
        "k--", label='Fitted Michaelis-Menten equation')
plt.text(5, 1, u"R\u00b2 score: {:0.2f}".format(r_squared), style='italic')
plt.xlabel('PNPP (μM)')
plt.ylabel('Velocity (nmoles/min/mg)')

plt.tight_layout()
plt.legend(loc='lower right')
plt.savefig('./Figures_22/RateVelocitykinetics_MM_2_paper.pdf', dpi=300)
plt.show()

# --- Estimation of Km and Vmax using the Lineweaver-Burk plot
substrate_conc = amounts_1[2:]
v_values = means[2:]

plot_lineweaverburk(substrate_conc, v_values)
#plt.savefig('./Figures_22/LineweaverBurk_kinetics_MM_final.pdf')

plt.show()


#%% Screening of other substrates for Psip1

substrates = pd.read_excel('Psip1_experiment_pag65_table.xlsx', sheet_name='Data', index_col=0)
means = np.array(substrates['pmoles/min/mg protein']); error = np.array(substrates['Std_pmoles'])
x_pos = np.arange(len(np.array(substrates.index)))
figure, ax = plt.subplots()
plt.rcParams["axes.labelweight"] = "bold"
ax.bar(x_pos, means, yerr=error, align='center', ecolor='black', capsize=10, edgecolor='black', color='white')
plt.axhline(y=0,linewidth= 1, color='k')
ax.set_ylabel('pmoles/min/mg of protein')
ax.set_xticks(x_pos, weight='bold')
ax.set_xticklabels(np.array(substrates.index))
ax.tick_params(axis='x', rotation=45)
ax.yaxis.grid(False)
plt.tight_layout()
plt.savefig('./Figures_22/pmoles_malakitegreen_paper.pdf', dpi=300)
plt.show()