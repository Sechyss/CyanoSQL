import os

import numpy as np
import pandas as pd

os.chdir('/Users/u1866168/Documents/OneDrive - University of Warwick/MHC_Scripts/')

cardionlipin_table = pd.read_excel('ValvanoData_Alberto.xlsx', index_col=0, sheet_name='Cardiolipin')
pg_table = pd.read_excel('ValvanoData_Alberto.xlsx', index_col=0, sheet_name='pg')
apg_table = pd.read_excel('ValvanoData_Alberto.xlsx', index_col=0, sheet_name='APG')
pe_table = pd.read_excel('ValvanoData_Alberto.xlsx', index_col=0, sheet_name='PE')


wt = np.array(cardionlipin_table['WT'])
wt_pmb = np.array(cardionlipin_table['WT+Pmb'])
lcoA = np.array(cardionlipin_table['lcoA'])
lcoA_pmb = np.array(cardionlipin_table['lcoA+Pmb'])
psrA = np.array(cardionlipin_table['psrA'])
psrA_pmb = np.array(cardionlipin_table['psrA+pmb'])
bcnA = np.array(cardionlipin_table['bcnA'])
bcnA_pmb = np.array(cardionlipin_table['bcnA+pmb'])
