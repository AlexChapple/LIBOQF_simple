"""
Revised version of the waiting time distribution calculator 
"""

import numpy as np 
import matplotlib.pyplot as plt 
import cmath as cm 

# emission_data = np.loadtxt("testing/emission_tracking.txt")
emission_data = np.loadtxt("results2/emission_tracking_23.txt")
# emission_data = np.loadtxt("testing/emission_test_sample.txt")

num_of_simulations = 50000
total_emissions = 1174201
end_time = 60
bin_increment = 10
bin_width = 0.2/bin_increment

waiting_time_list = np.zeros(int(np.ceil(end_time/bin_width)))
reduced_time_list = np.linspace(0,end_time,int(np.ceil(end_time/bin_width)))

last_time = None 
simulation_counter = 0 
waiting_time_counter = 0 
new_sim = True 

for emission in emission_data:

    if emission == end_time + 50:
        new_sim = True
        simulation_counter += 1
    
    else:

        if new_sim:

            last_time = emission
            new_sim = False

        else:

            waiting_time = emission - last_time
        
            index = int(np.floor(waiting_time/bin_width))

            waiting_time_list[index] += 1

            last_time = emission

            waiting_time_counter += 1

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

plt.plot(reduced_time_list, waiting_time_list)
plt.plot(reduced_time_list, analytical_list)

print(max(waiting_time_list) / max(analytical_list))
print(np.size(waiting_time_list))
plt.xlim([0,15])

plt.show()