"""
Plot the g2 plot 
"""

import numpy as np 
import matplotlib.pyplot as plt 
import cmath as cm 
import colours
import matplotlib
matplotlib.rcParams.update({'font.size': 24})

Omega_list = [2.3, 1.0, 0.7]

time_list = np.loadtxt("results2/g2_23.txt")[:,0]

g2_data = []
for Omega in Omega_list:
    
    name = "results2/g2_" + str(Omega).replace(".", "") + ".txt"
    g2_temp = np.loadtxt(name)[:,1]

    name2 = "results2/avg_e_" + str(Omega).replace(".", "") + ".txt"
    norm_temp = np.loadtxt(name2)[:,1]

    norm = sum(norm_temp) / np.size(norm_temp)

    g2_temp /= (norm**2)

    g2_data.append(g2_temp)

# Set up analytical results 
def g2(t, Omega):

    Y = np.sqrt(2) * Omega
    delta = (1/4) * cm.sqrt(1 - 8*Y**2)

    a = np.exp(-3*t/4) * (np.cosh(delta * t) + ((3/(4*delta)) * np.sinh(delta*t)))

    return 1 - a

analytical_data = []
for Omega in Omega_list:

    g2_analytical_temp = [g2(t, Omega) for t in time_list]

    analytical_data.append(g2_analytical_temp)



# Plot numerical results 

fig1 = plt.figure(1)
fig1.set_size_inches(18.5, 10.5)
fig1.set_alpha(0)
fig1.set_facecolor("none")
ax1 = plt.axes()
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_color(colours.spanish_gray)
ax1.spines['bottom'].set_color(colours.spanish_gray)
ax1.tick_params(axis='x', colors=colours.spanish_gray)
ax1.tick_params(axis='y', colors=colours.spanish_gray)
ax1.xaxis.label.set_color(colours.spanish_gray)
ax1.yaxis.label.set_color(colours.spanish_gray)

colour_list = [colours.greek_dark_red, colours.orange_peel, colours.new_green]


for Omega_index in range(len(Omega_list)):

    plt.plot(time_list, g2_data[Omega_index], label="$\Omega$ = " + str(Omega_list[Omega_index]).replace(".", ""), lw=5, c=colour_list[Omega_index])

    if Omega_index == 0:
        plt.plot(time_list, analytical_data[Omega_index], lw=5, label="Analytical", c=colours.greek_dark_blue)     
    else:
        plt.plot(time_list, analytical_data[Omega_index], lw=5, c=colours.greek_dark_blue)     
    
plt.show()