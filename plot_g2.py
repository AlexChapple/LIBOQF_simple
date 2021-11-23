import numpy as np 
import matplotlib.pyplot as plt 
import cmath as cm

from numpy.lib.function_base import iterable 

g2_data = np.loadtxt("results_g2/g2.txt")

emissions = 88675
iterations = 20000
end_time = 10 
time_list_length = 5000 
g2_length = 500



# divider = 0.2075
divider = 1/5

g2_data /= divider

time_list = np.linspace(0,10,5000)

def g2(t, Omega):

    Y = np.sqrt(2) * Omega
    delta = (1/4) * cm.sqrt(1 - 8*Y**2)

    a = np.exp(-3*t/4) * (np.cosh(delta * t) + ((3/(4*delta)) * np.sinh(delta*t)))

    return 1 - a

g2_analytical = [g2(t, 2.3) for t in time_list]

plt.plot(time_list, g2_data, label="numerical")
plt.plot(time_list, g2_analytical, label="analytical")
plt.legend()
plt.show()