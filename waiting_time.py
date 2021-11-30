"""

Generates a waiting time distribution 

"""

import numpy as np 
import matplotlib.pyplot as plt 
import colours
import matplotlib
import cmath as cm 
matplotlib.rcParams.update({'font.size': 24})

directory = "results/"

# Import data 
emission_data = np.loadtxt("emission_tracking.txt")
num_of_simulations = 3000 
end_time = 20 
bin_width = (1/500)

waiting_time_list = np.zeros(int(np.ceil(end_time/bin_width)))
reduced_time_list = np.linspace(0,end_time,int(np.ceil(end_time/bin_width)))

# Plot waiting time distribution 

last_time = None 
simulation_counter = 0 
new_sim = True 

for emission in emission_data:

    if emission == end_time + 50:

        new_sim = True 
        simulation_counter += 1 

    else:

        if new_sim == True: 

            last_time = emission 
            new_sim = False 

        else: 

            waiting_time = emission - last_time

            index = int(np.floor(waiting_time / bin_width))

            waiting_time_list[index] += 1 

            last_time = emission



### ----- Plotting ------------------------------

# Plots the total waiting time distribution 

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

waiting_time_norm = [i / num_of_simulations for i in waiting_time_list]

# Analytical result here 

Omega = 2.3
def waiting_analytical(t):

    Y = np.sqrt(2) * Omega
    delta_prime  = (1/2) * cm.sqrt(1 - 2*Y**2)

    a = np.exp(-t/2) * (Y**2 / (2*Y**2 - 1)) * (1 - np.cosh(delta_prime * t))

    return a 

analytical_list = [waiting_analytical(t) for t in reduced_time_list]

plt.bar(reduced_time_list, waiting_time_norm, width=0.0025, color=colours.greek_blue)
plt.plot(reduced_time_list, analytical_list, c="red")

plt.xlabel("Waiting time (seconds)")
plt.ylabel("Frequency (normalised)")
plt.xlim([0,15])
plt.savefig(directory + "waiting_time.pdf", facecolor=fig1.get_facecolor(), transparent=True, dpi=600)
