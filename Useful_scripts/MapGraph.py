import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
data = pd.read_excel('/Users/u2176312/Downloads/PhosphatasesAMT.xlsx', sheet_name='AMT23', index_col='Station')
gdf = gpd.GeoDataFrame(data, geometry=gpd.points_from_xy(data.Longitude, data.Latitude))

fig = plt.figure()
ax_map = fig.add_axes([0, 0, 1, 1])
world.plot(ax=ax_map, color='grey', alpha=0.6)


# function to create inset axes and plot bar chart on it
# this is good for 3 items bar chart
def build_pie(mapx, mapy, ax, width, vals: list, fcolors: list, size):
    ax_h = inset_axes(ax, width=width,
                      height=width,
                      loc=3,
                      bbox_to_anchor=(mapx, mapy),
                      bbox_transform=ax.transData,
                      borderpad=0,
                      axes_kwargs={'alpha': 0.35, 'visible': True})

    ax_h.pie(vals, colors=fcolors, radius=size)
    ax_h.axis('off')
    return ax_h


bar_width = 0.1  # inch
colors = ['#1b9e77', '#D95F02', '#7570B3']
for index, row in gdf.iterrows():
    x1, y1 = gdf.loc[index]['Longitude'], gdf.loc[index]['Latitude']  # get data coordinates for plotting
    bax = build_pie(x1, y1, ax_map, 0.2, vals=[gdf.loc[index]['phoD'], gdf.loc[index]['phoX'],
                                               gdf.loc[index]['psip1']],
                    fcolors=colors, size=gdf.loc[index]['SizePieChart'])

patch0 = mpatches.Patch(color=colors[0], label='$\mathit{phoD}$')
patch1 = mpatches.Patch(color=colors[1], label='$\mathit{phoX}$')
patch2 = mpatches.Patch(color=colors[2], label='$\mathit{psip1}$')
legend = ax_map.legend(handles=[patch0, patch1, patch2], loc=1)
legend.get_frame().set_alpha(None)
plt.savefig('/Users/u2176312/OneDrive - University of '
            'Warwick/Thesis/Paper Draft/PlotAMT23DATA_Phosphatases_paper.png', dpi=300)

plt.show()
