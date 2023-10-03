import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import pandas as pd
import scipy.interpolate as interp

# Indices
# 0:  Collision altitude [m]
# 1:  #Fragments [-]
# 2:  T_0 [days]
# 3:  SC Altitude Offset [m]
# 4:  Target Fraction [-]
# 5:  FoV [deg]
# 6:  Range [m]
# 7:  Incidence Angle [deg]
# 8:  Ablation Time [s]
# 9:  Scan Time [s]
# 10: Cooldown Time [s]
# 11: Fluence [J/m^2]
# 12: Removal Altitude [m]
# 13: Time Required [days]
# 14: Fraction removed [-]
# data = np.genfromtxt("runs.csv", delimiter=",", names=True, comments=None)
# # data = np.genfromtxt("runs.csv", delimiter=",", skip_header=1)
# # header = np.genfromtxt("runs.csv", dtype=str, delimiter=",", max_rows=1)
# print(data)
# data[data[:,14]<0.5] = np.nan # Any run that ran into t_max should have its fraction treated as NaN

df = pd.read_csv("runs.csv")
# df.loc[df['Fraction removed [-]'] < 0.5,'Time Required [days]'] = np.nan

fig, ax = plt.subplots(1, 2, figsize=[10,5])

x_axis = "Ablation Time [s]"
y_axis = "Range [m]"

# grid_x, grid_y = np.meshgrid(np.linspace(df[x_axis].min(), df[x_axis].max(), 1000), np.linspace(df[y_axis].min(), df[y_axis].max(), 1000))
# grid_z = interp.griddata(np.vstack([df[x_axis], df[y_axis]]).T, df["Time Required [days]"], (grid_x, grid_y), method='nearest')
# plt.imshow(grid_z, origin='lower', cmap="viridis_r", extent=[df[x_axis].min(), df[x_axis].max(), df[y_axis].min(), df[y_axis].max()], aspect="auto")

scatter = ax[0].scatter(df[x_axis], df[y_axis], s=10, c=df["Time Required [days]"], cmap="viridis_r")
contour = ax[1].tricontourf(df[x_axis], df[y_axis], df["Time Required [days]"], cmap="viridis_r")
ax[0].set_aspect("auto")
ax[1].set_aspect("auto")
ax[0].set_ylabel(y_axis)
ax[0].set_xlabel(x_axis)
ax[1].set_xlabel(x_axis)

plt.colorbar(contour, label="Removal time [days]", cax=inset_axes(ax[1],
                    width="5%",  
                    height="100%",
                    loc='right',
                    borderpad=-2
))
plt.show()