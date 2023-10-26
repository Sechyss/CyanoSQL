import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

data = pd.read_excel('/Users/u1866168/Desktop/VirulenceDataWT_MT_NP.xlsx', sheet_name='Sheet1', index_col='Time')
fig, ax = plt.subplots()
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams.update({'font.size': 22})
plt.figure(figsize=(10, 10))
plt.rcParams['axes.xmargin'] = 0

x = np.array(data.index)

y_WT = np.array(data['T4 Wild Type'])
y_mutant = np.array(data['T4 ∆tRNA::HygR mutant phage'])
y_nP = np.array(data['No Phage'])
st_wt = np.array(data['WT_SD'])
st_mt = np.array(data['MT_SD'])
st_np = np.array(data['NP_SD'])


plt.axvline(x[3], 0, y_nP.max(), color='#838093')
plt.axvline(x[8], 0, y_nP.max(), color='#B1C2BE')

a = plt.plot(x, y_nP, color='#70A0AF', label='No Phage')
plt.fill_between(x, y_nP - st_np, y_nP + st_np,
                 alpha=0.5, edgecolor='#70A0AF', facecolor='#8EA5AC')


b = plt.plot(x, y_mutant, color='#A0C1B9', label='T4 ∆tRNA::HygR mutant phage')
plt.fill_between(x, y_mutant - st_mt, y_mutant + st_mt,
                 alpha=0.5, edgecolor='#A0C1B9', facecolor='#B1C2BE')

c = plt.plot(x, y_WT, color='#706993', label='T4 Wild Type')
plt.fill_between(x, y_WT - st_wt, y_WT + st_wt,
                 alpha=0.5, edgecolor='#706993', facecolor='#838093')

plt.legend(loc='upper right')
plt.xlabel('Time (mins)')
plt.ylabel('Abs(650nm)')

start, end = ax.get_xlim()

ax.xaxis.set_ticks(np.arange(0, 2001, 200))
plt.tight_layout()
plt.savefig('/Users/u1866168/Desktop/GrowthDataWT_MT_noSTDev_2.png')
plt.show()
