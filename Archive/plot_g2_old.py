"""

    Plots only g2 

"""

import numpy as np 
import matplotlib.pyplot as plt 
import colours
import matplotlib
import cmath as cm 
matplotlib.rcParams.update({'font.size': 24})


directory = "results/"

# Import data 
emission_data = np.loadtxt(directory + "emission_tracking_10.txt")
num_of_simulations = 1 
end_time = 10 
bin_width = (1/50)

g2_list = np.zeros(int(np.ceil(end_time/bin_width)))
reduced_time_list = np.linspace(0,end_time,int(np.ceil(end_time/bin_width)))

# Plot waiting time distribution 

sim_init_time = None 
simulation_counter = 0 
emission_counter = 0 
new_sim = True 

for emission in emission_data:

    if emission == end_time + 50:

        new_sim = True 
        simulation_counter += 1 

    else:

        if new_sim == True: 

            sim_init_time = emission 
            new_sim = False 
            emission_counter += 1 

        else: 

            waiting_time = emission - sim_init_time

            index = int(np.floor(waiting_time / bin_width))

            g2_list[index] += 1 
            emission_counter += 1 

g2_list = g2_list / (emission_counter / (end_time / bin_width))

print(simulation_counter)
print(emission_counter)


### ----- Plotting ------------------------------

# Plots the total waiting time distribution 



### ----- Analytical result --------------------------------

Omega_list = [0.7]

def g2(t, Omega):

    Y = np.sqrt(2) * Omega
    delta = (1/4) * cm.sqrt(1 - 8*Y**2)

    a = np.exp(-3*t/4) * (np.cosh(delta * t) + ((3/(4*delta)) * np.sinh(delta*t)))

    return 1 - a

g2_list_analytical = []

for Omega in Omega_list:

    temp_list = [g2(t, Omega) for t in reduced_time_list]
    g2_list_analytical.append(temp_list)


### ----- Plotting ------------------------------------

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

plt.plot(reduced_time_list, g2_list, c=colours.greek_blue)

for index in range(len(Omega_list)):

    if index == 1:
        plt.plot(reduced_time_list, g2_list_analytical[index], lw=5, c=colours.greek_dark_blue, ls="dotted", label="analytical")
    else:
        plt.plot(reduced_time_list, g2_list_analytical[index], lw=5, c=colours.greek_dark_blue, ls="dotted")

plt.xlabel("Time $\gamma t$ (seconds)")
plt.ylabel("$g^{(2)}(t)$")
plt.grid()
plt.legend(loc="lower right")
# plt.xlim([0,13])
# plt.savefig(directory + "g2.pdf", facecolor=fig.get_facecolor(), transparent=True, dpi=600)
plt.show()
