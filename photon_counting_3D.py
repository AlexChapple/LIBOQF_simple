"""

Generates a 3D plot of the photon counting distribution for resonance fluorescence 

"""

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib 
matplotlib.rcParams.update({'font.size': 24})
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Import data 

directory = 'photon_count_3D_results/'

file_list = range(5,105,5)
data_list = np.zeros((np.size(file_list),100)) # THIS HAS TO BE CHANGED TO SUPPORT MESH HERE 
photon_list = range(0,101)

for i in range(len(file_list)):
    
    data_list[i] = np.loadtxt(directory + "photon_counting_" + str(file_list[i]) + ".txt")

X,Y = np.meshgrid(photon_list, file_list)

ax.plot_surface(X,Y,data_list) 

plt.show()