"""
Plots g2
"""

import numpy as np 
import matplotlib.pyplot as plt 
import colours
import matplotlib
import cmath as cm 
matplotlib.rcParams.update({'font.size': 24})

directory = "results/"
end_time = 20 
bin_width = (1/50)

# Import data 
Omega_list = [0.7, 1.0, 2.3]

emission_data = []
reduced_time_list = np.linspace(0,end_time,int(np.ceil(end_time/bin_width)))

for Omega in Omega_list:

    emission_data.append(np.loadtxt(directory + "emission_tracking_" + str(Omega).replace(".", "") + ".txt"))


g2_list = []

for Omega_index in range(len(Omega_list)):

    temp_g2_list = np.zeros(int(np.ceil(end_time/bin_width)))

    simulation_initial_time = None 
    emission_counter = 0 
    new_sim = True 

    for emission in emission_data[Omega_index]:

        if emission == end_time + 50:

            new_sim = True 
            
        else:

            if new_sim == True:

                simulation_initial_time = emission
                new_sim = False 
                emission_counter += 1 

            else:

                waiting_time = emission - simulation_initial_time

                index = int(np.floor(waiting_time/bin_width))

                temp_g2_list[index] += 1
                emission_counter += 1

    temp_g2_list /= (emission_counter / (end_time / bin_width))

    print(np.size(temp_g2_list))

    g2_list.append(temp_g2_list)


### ----- Analytical result -----------------------------------

def g2(t, Omega):

    Y = np.sqrt(2) * Omega
    delta = (1/4) * cm.sqrt(1 - 8*Y**2)

    a = np.exp(-3*t/4) * (np.cosh(delta * t) + ((3/(4*delta)) * np.sinh(delta*t)))

    return 1 - a

g2_list_analytical = []

for Omega in Omega_list:

    temp_list = [g2(t, Omega) for t in reduced_time_list]
    g2_list_analytical.append(temp_list)
    print(np.size(temp_list))


### ----- Plotting ----------------------------------------------

fig = plt.figure()
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

colour_list = [colours.greek_dark_red, colours.orange_peel, colours.new_green]

for index in range(len(Omega_list)):

    if index == 2:

        plt.plot(reduced_time_list, g2_list[index], lw=3, c=colour_list[index], label="$\Omega = $" + str(Omega_list[index]))

        if index == 2:
            plt.plot(reduced_time_list, g2_list_analytical[index], lw=3, c=colours.greek_dark_blue, ls="dotted", label="analytical")
        else:
            plt.plot(reduced_time_list, g2_list_analytical[index], lw=3, c=colours.greek_dark_blue, ls="dotted")


plt.xlabel("Time (seconds)")
plt.ylabel("$g^{(2)}(t)$")
plt.grid()
plt.legend()
# plt.legend(loc="lower right")
plt.savefig(directory + "g2.pdf", facecolor=fig.get_facecolor(), transparent=True, dpi=600)
