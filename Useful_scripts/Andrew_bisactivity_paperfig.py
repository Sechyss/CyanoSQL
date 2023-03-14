import matplotlib.pyplot as plt
import pandas as pd

if __name__ == '__main__':
    Diester = pd.read_excel('/Users/u2176312/OneDrive - University of Warwick/Thesis/Paper '
                            'Draft/Diesterase_test.xlsx', sheet_name='Sheet2', header=0)

    x_values = [x for x in Diester.columns if '_std' not in x]
    error_columns = [x for x in Diester.columns if '_std' in x]
    means = Diester[x_values].values.tolist()
    means = [item for sublist in means for item in sublist]
    error = Diester[error_columns].values.tolist()
    error = [item for sublist in error for item in sublist]

    figure, ax = plt.subplots()
    plt.rcParams["axes.labelweight"] = "bold"
    x_values = [x.replace('F. johnsoniae', '$\it{F. johnsoniae}$') for x in x_values]
    ax.bar(x_values, means, yerr=error, align='center', ecolor='black', capsize=10, edgecolor='black', color='white')
    ax.set_ylabel('Abs (405nm)', weight='bold')
    ax.set_xticks(x_values)
    ax.tick_params(axis='x', rotation=45)
    ax.set_xticklabels(x_values)
    ax.yaxis.grid(False)

    ax.axhline(y=0, color='black', linestyle='-')

    # Save the figure and show
    plt.tight_layout()
    plt.savefig('/Users/u2176312/OneDrive - University of Warwick/Thesis/Paper Draft/Diesterase_test.pdf', dpi=300)
    plt.show()
