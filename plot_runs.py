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

df = pd.read_csv("runs.csv")
fig, ax = plt.subplots(1, 2, figsize=[10,5])

# Select dataset to plot, set x_axis and y_axis and set values for all other parameters that are constant in the plot (only X and Y axis are plotted obviously)
x_axis = "Ablation Time [s]"
y_axis = "Range [m]"
df = df[
    np.isclose(df["Collision altitude [m]"], 789000) &
    np.isclose(df["#Fragments [-]"], 100000) &
    np.isclose(df["T_0 [days]"], 1) &
    np.isclose(df["SC Altitude Offset [m]"], 30e3) &
    np.isclose(df["FoV [deg]"], 38.44) &
    np.isclose(df["Target Fraction [-]"], 0.5) &
    # np.isclose(df["Range [m]"], 300000) &
    np.isclose(df["Incidence Angle [deg]"], 20) &
    # np.isclose(df["Ablation Time [s]"], 70) &
    np.isclose(df["Scan Time [s]"], 5) &
    np.isclose(df["Cooldown Time [s]"], 70) &
    np.isclose(df["Fluence [J/m^2]"], 8500) &
    np.isclose(df["Removal Altitude [m]"], 340000)
]

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