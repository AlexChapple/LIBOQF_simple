"""
Tests to see if the revised simulation is working fine 
"""

import numpy as np 
import matplotlib.pyplot as plt 
import cmath as cm 

# z_data = np.loadtxt("testing/sigma_z_test.txt")
# L_data = np.loadtxt("testing/sigma_L_test.txt")
# R_data = np.loadtxt("testing/sigma_R_test.txt")

z_data = np.loadtxt("results2/sigma_z_23.txt")
L_data = np.loadtxt("results2/sigma_L_23.txt")
R_data = np.loadtxt("results2/sigma_R_23.txt")

time_list = z_data[:,0]
z_list = z_data[:,1]
L_list = L_data[:,1]
R_list = R_data[:,1]

# Analytical result 
def sigma_z(t, Omega ):

    Y = np.sqrt(2) * Omega

    delta = (1/4) * cm.sqrt(1 - 8*Y**2)

    a = -1 * (1/(1 + Y**2))
    b = Y**2 * np.exp(-3*t/4) 
    c = np.cosh(delta*t) + (3/(4*delta))*np.sinh(delta*t)

    return np.real(a * (1 + b*c))

z_list_analytical = [sigma_z(t, 2.3) for t in time_list]

plt.figure(1)
plt.plot(time_list, z_list)
plt.plot(time_list, z_list_analytical, ls="dashed")


def sigma_LR(t, mode, Omega):

    p = 1

    if mode == "R":
        p = -1

    Y = np.sqrt(2) * Omega
    delta = (1/4) * cm.sqrt(1 - 8*Y**2)

    a = p * (1.0j / np.sqrt(2)) * (Y / (1 + Y**2))
    b = np.exp(-3*t/4) * (np.cosh(delta*t) + (3/(4*delta))*np.sinh(delta*t))
    c = p * 1.0j * np.sqrt(2) * Y * np.exp(-3*t/4) * (1/(4*delta)) * np.sinh(delta*t)

    return np.imag(a * (1 - b) + c)

L_list_analytical = [sigma_LR(t, "L", 2.3) for t in time_list]

plt.figure(2)
plt.plot(time_list, L_list)
plt.plot(time_list, L_list_analytical, ls="dashed")

R_list_analytical = [sigma_LR(t, "R", 2.3) for t in time_list]

plt.figure(3)
plt.plot(time_list, R_list)
plt.plot(time_list, R_list_analytical, ls="dashed")

print(sum(z_list)/np.size(z_list))

plt.show()