import numpy as np 
import matplotlib.pyplot as plt 
import cmath as cm 

g2_data = np.loadtxt("results_g2/g2.txt")

time_list = np.linspace(0,10,200)

def g2(t, Omega):

    Y = np.sqrt(2) * Omega
    delta = (1/4) * cm.sqrt(1 - 8*Y**2)

    a = np.exp(-3*t/4) * (np.cosh(delta * t) + ((3/(4*delta)) * np.sinh(delta*t)))

    return 1 - a

g2_analytical = [g2(t, 1.0) for t in time_list]

plt.plot(time_list, g2_data)
plt.plot(time_list, g2_analytical)
plt.show()