"""
Plot the g2 plot
"""

import numpy as np 
import matplotlib.pyplot as plt 
import cmath as cm 

data = np.loadtxt("results2/g2_23.txt")

time_list = data[:,0]
g2_list = data[:,1] 

end_time = 60
num_of_simulations = 50000 
emission = 679658

def g2(t, Omega):

    Y = np.sqrt(2) * Omega
    delta = (1/4) * cm.sqrt(1 - 8*Y**2)

    a = np.exp(-3*t/4) * (np.cosh(delta * t) + ((3/(4*delta)) * np.sinh(delta*t)))

    return 1 - a

g2_analytical = [g2(t, 2.3) for t in time_list]

ratio = np.real(max(g2_list) / max(g2_analytical))
print(ratio)

g2_list /= ratio 





plt.figure(1)
plt.plot(time_list, g2_list)
plt.plot(time_list, g2_analytical)


# Shows if the features are the same to make sure 
g2_list /= max(g2_list)
g2_analytical /= max(g2_analytical)

plt.figure(2)
plt.plot(time_list, g2_list)
plt.plot(time_list, g2_analytical)

plt.show()