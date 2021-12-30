"""
Revised version of the waiting time distribution calculator 
"""

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
import cmath as cm 
matplotlib.rcParams.update({'font.size': 24})
import colours

# emission_data = np.loadtxt("testing/emission_tracking.txt")
emission_data = np.loadtxt("results2/emission_tracking_23.txt")
# emission_data = np.loadtxt("testing/emission_test_sample.txt")

directory = "results2/"

num_of_simulations = 100000
total_emissions = 1174201
end_time = 100
bin_increment = 20
bin_width = 0.2/bin_increment

waiting_time_list = np.zeros(int(np.ceil(end_time/bin_width)))
reduced_time_list = np.linspace(0,end_time,int(np.ceil(end_time/bin_width)))

last_time = None 
simulation_counter = 0 
waiting_time_counter = 0 
new_sim = True 
emission_counter = 0

for emission in emission_data:

    if emission == end_time + 50:
        new_sim = True
        simulation_counter += 1
    
    else:

        if new_sim:

            last_time = emission
            new_sim = False
            emission_counter += 1

        else:

            waiting_time = emission - last_time
        
            index = int(np.floor(waiting_time/bin_width))

            waiting_time_list[index] += 1

            last_time = emission

            waiting_time_counter += 1
            emission_counter += 1

print(waiting_time_counter)
waiting_time_list /= (waiting_time_counter / (np.size(reduced_time_list) / end_time))

# Plot waiting time distribution here 

Omega = 2.3
def waiting_analytical(t):

    Y = np.sqrt(2) * Omega
    delta_prime  = (1/2) * cm.sqrt(1 - 2*Y**2)

    a = np.exp(-t/2) * (Y**2 / (2*Y**2 - 1)) * (1 - np.cosh(delta_prime * t))

    return a 

analytical_list = [waiting_analytical(t) for t in reduced_time_list]
# analytical_list /= max(analytical_list)

# coherent light waiting time 
nbar = emission_counter / num_of_simulations
tau_bar = end_time / nbar 
print(tau_bar)

# tau_bar = 1
coh_waiting_list = [1/tau_bar * np.exp(-t/tau_bar) for t in reduced_time_list]
# total = sum(coh_waiting_list)
# coh_waiting_list = [i/total for i in coh_waiting_list]
# print(tau_bar, sum(coh_waiting_list))
print(sum(coh_waiting_list))

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
plt.xlim([0,15])

plt.bar(reduced_time_list, waiting_time_list, color=colours.greek_blue, label="numerical", width=0.015)
plt.plot(reduced_time_list, analytical_list, c=colours.greek_dark_red, label="analytical", lw=4)
plt.plot(reduced_time_list, coh_waiting_list, lw=4, c=colours.new_green, ls="dashed", label="coherent light")

plt.xlabel("$\gamma \\tau$")
plt.ylabel("Waiting time distribution $w(\\tau)$")
plt.legend()
plt.grid()
plt.savefig(directory + "waiting_time.pdf", dpi = 600)