# Upload of packages to run the script
import argparse
import textwrap
import warnings

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

warnings.simplefilter(action='ignore', category=FutureWarning)
# =============================================================================

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
Creation of Heatmap.
-------------------------------------------------------------
'''))
# Parse command line arguments
# -i csvfile -o OUTPUT

parser.add_argument("-i", metavar='file.csv', dest="csv", help="Input", type=str)
parser.add_argument("-o", metavar='file.png', dest="output", help="Output", type=str)

args = parser.parse_args()

table = pd.read_csv(args.csv, sep=',', index_col=0)
figure, ax = plt.subplots()
sns.heatmap(table, annot=False, cmap="vlag", xticklabels=True,
            yticklabels=True, ax=ax, fmt='d')

ax.xaxis.tick_top()
ax.tick_params(axis='x', which='major', labelsize=3)
ax.tick_params(axis='y', which='major', labelsize=1)
plt.savefig(args.output, dpi=300)
