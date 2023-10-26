import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv('/Users/u1866168/Documents/OneDrive - University of '
                 'Warwick/Experiments/Analysis_Genes/Psip1/Activity_Experiments/DataCorrelation.csv',
                 index_col='Time(mins)')

x = np.array(df.index)
FeC = df[df.filter(like='FeC').columns]
Mn = df[df.filter(like='Mn').columns]
Mg = df[df.filter(like='Mg').columns]
Zn = df[df.filter(like='Zn').columns]
Ca = df[df.filter(like='Ca').columns]

collector = {}

for row in df:
    y = np.array(df[row])
    idx = np.isfinite(x) & np.isfinite(y)
    ab = np.polyfit(x[idx], y[idx], 1)
    if str(row).rsplit("_", 1)[0] in collector:
        collector[str(row).rsplit("_", 1)[0]].update({str(row).rsplit("_", 1)[1]: ab[0]})
    else:
        collector.update({str(row).rsplit("_", 1)[0]: {str(row).rsplit("_", 1)[1]: ab[0]}})


MgX = np.array(list(str(x).replace('mM', '') for x in collector['Mg'].keys()), dtype='float')
MgY = np.array(list(collector['Mg'].values()))

MgX = np.delete(MgX, [1,2], None)
MgY = np.delete(MgY, [1,2], None)

plt.figure(figsize=(10, 12))
plt.subplot(211)
plt.plot(x, Mg['Mg_0mM'],x, Mg['Mg_0.1mM'], x, Mg['Mg_1mM'],
         x, Mg['Mg_10mM'], x, Mg['Mg_30mM'], x, Mg['Mg_60mM'], x, Mg['Mg_100mM'])
plt.xlabel('Time (mins)')
plt.ylabel('Abs 405nm')
plt.legend(list(Mg.columns), loc='upper left')

plt.subplot(212)
plt.plot(MgX, MgY, 'b-o')
plt.xlabel('Concentration (mM)')
plt.ylabel('R coefficient')

plt.show()

