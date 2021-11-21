"""

Generates a 3D plot of the photon counting distribution for resonance fluorescence 

"""

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib 
from matplotlib import cm 
import colours
matplotlib.rcParams.update({'font.size': 14})
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Import data 

directory = 'photon_count_3D_results/'

file_list = range(2,102,2)
data_list = [] # THIS HAS TO BE CHANGED TO SUPPORT MESH HERE 
photon_list = range(0,100)

for i in range(len(file_list)):

    data = np.loadtxt(directory + "photon_counting_" + str(file_list[i]) + ".txt")
    
    data_list.append(data)

X,Y = np.meshgrid(photon_list, file_list)

fig = plt.figure(1)
fig.set_size_inches(18.5, 10.5)
fig.set_alpha(0)
fig.set_facecolor("none")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_color(colours.spanish_gray)
ax.spines['bottom'].set_color(colours.spanish_gray)
# ax.xaxis.label.set_color(colours.spanish_gray)
# ax.yaxis.label.set_color(colours.spanish_gray)
# ax.zaxis.label.set_color(colours.spanish_gray)

ax.plot_surface(X,Y,np.array(data_list), cmap=cm.turbo) 
ax.set_xlabel("\nphoton number")
ax.set_ylabel("\n$\gamma t$")
ax.view_init(elev=20, azim=70)
ax.tick_params(axis='x', colors=colours.spanish_gray)
ax.tick_params(axis='y', colors=colours.spanish_gray)
ax.tick_params(axis='z', colors=colours.spanish_gray, pad=10)
# ax.set_zlabel("\n\n\nPhoton counting distribution", fontsize=11)

# # Code for wireframe plot 
# ax2 = fig.add_subplot(111, projection='3d')
# fig2 = plt.figure(2)
# fig2.set_size_inches(18.5, 10.5)
# fig2.set_alpha(0)
# fig2.set_facecolor("none")

# ax2.spines['top'].set_visible(False)
# ax2.spines['right'].set_visible(False)
# ax2.spines['left'].set_color(colours.spanish_gray)
# ax2.spines['bottom'].set_color(colours.spanish_gray)
# ax2.xaxis.label.set_color(colours.spanish_gray)
# ax2.yaxis.label.set_color(colours.spanish_gray)
# ax2.zaxis.label.set_color(colours.spanish_gray)

# ax2.plot_wireframe(X,Y,np.array(data_list), rstride=2, cstride=2) 
# ax2.set_xlabel("\nphoton number")
# ax2.set_ylabel("\n$\gamma t$")
# ax2.view_init(elev=32, azim=67)
# ax2.tick_params(axis='x', colors=colours.spanish_gray)
# ax2.tick_params(axis='y', colors=colours.spanish_gray)
# ax2.tick_params(axis='z', colors=colours.spanish_gray, pad=10)

plt.savefig("results/photon_counting_3D.pdf", facecolor=fig.get_facecolor(), transparent=True, dpi=600, bbox_inches='tight', pad_inches=0)