import geopandas as gpd
import pandas as pd
from matplotlib import pyplot as plt, font_manager
import matplotlib.patches as mpatches


world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
data = pd.read_excel('/Users/u2176312/Downloads/PhosphatasesAMT.xlsx', sheet_name='AMT22', index_col='Station')
data2 = pd.read_excel('/Users/u2176312/Downloads/PhosphatasesAMT.xlsx', sheet_name='AMT23', index_col='Station')

gdf = gpd.GeoDataFrame(data, geometry=gpd.points_from_xy(data.Longitude, data.Latitude))
gdf2 = gpd.GeoDataFrame(data2, geometry=gpd.points_from_xy(data2.Longitude, data2.Latitude))


fig = plt.figure()
ax_map = fig.add_axes([0, 0, 1, 1])
world.plot(ax=ax_map, color='grey', alpha=0.6)
gdf.plot(ax=ax_map, color='red', alpha=0.6, markersize=10, marker='o')
gdf2.plot(ax=ax_map, color='green', alpha=0.6, markersize=10, marker='*')

patch0 = mpatches.Patch(color='red', label='AMT22(JC079)')
patch1 = mpatches.Patch(color='green', label='AMT23(JR300)')
font = font_manager.FontProperties(family='Times new roman',  # 'Times new roman',
                                   weight='bold',
                                   style='normal')
legend = ax_map.legend(handles=[patch0, patch1], loc=1, prop=font)
legend.get_frame().set_alpha(None)
plt.savefig('/Users/u2176312/Downloads/AMTCruises.png', dpi=300)
plt.show()
