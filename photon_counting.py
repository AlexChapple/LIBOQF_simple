"""

Plots results for Little boxes simple 

"""

import numpy as np 
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import var 
import colours
import matplotlib
import cmath as cm 
matplotlib.rcParams.update({'font.size': 24})

directory = "results/"
photon_data = np.loadtxt("photon_counting.txt")
photon_bin_cut_off = 20 
num_of_simulations = 100000


# Plot photon counting 
fig = plt.figure(1)
fig.set_size_inches(18.5, 10.5)
fig.set_alpha(0)
fig.set_facecolor("none")
ax = plt.axes()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_color(colours.spanish_gray)
ax.spines['bottom'].set_color(colours.spanish_gray)
ax.tick_params(axis='x', colors=colours.spanish_gray)
ax.tick_params(axis='y', colors=colours.spanish_gray)
ax.xaxis.label.set_color(colours.spanish_gray)
ax.yaxis.label.set_color(colours.spanish_gray)

x_list = range(np.size(photon_data))
photon_data_norm = [i / num_of_simulations for i in photon_data]


plt.bar(x_list[0:photon_bin_cut_off], photon_data_norm[0:photon_bin_cut_off])
plt.xlabel("Photon number")
plt.ylabel("Probability")
plt.savefig(directory + "photon_counting.pdf", facecolor=fig.get_facecolor(), transparent=True, dpi=600)

# Do statistics here 
nbar = 0 
for i in range(len(x_list)):

    nbar += x_list[i] * photon_data_norm[i]

variance = 0 
for i in range(len(x_list)):

    variance += ((x_list[i] - nbar)**2) * photon_data_norm[i]

Mandel_Q = (variance - nbar)/ nbar 

print("$\\overbar{n}$ = ", nbar)
print("Q = ", Mandel_Q)

# Add a poisson distribution curve onto the plot 
poisson_list = []
for k in x_list:

    p = (nbar ** k) * np.exp(-nbar) / np.math.factorial(k)
    poisson_list.append(p)

plt.plot(x_list[0:photon_bin_cut_off], poisson_list[0:photon_bin_cut_off])
plt.savefig(directory + "photon_counting_poisson.pdf", facecolor=fig.get_facecolor(), transparent=True, dpi=600)
