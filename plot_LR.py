"""

Plots only sigma L,R 

"""

import numpy as np 
import matplotlib.pyplot as plt 
import colours
import matplotlib
import cmath as cm 
matplotlib.rcParams.update({'font.size': 24})

directory = "results/"
Omega_list = [2.3, 1.0, 0.7]


### ----- Import data --------------------------------------
L_list = []

for Omega in Omega_list:

    name = directory + "sigma_L_" + (str(Omega).replace(".", "")) + ".txt"

    L_list.append(np.loadtxt(name))

time_list = L_list[0][:,0]


### ----- Analytical ----------------------------------------

def sigma_LR(t, mode, Omega):

    p = 1

    if mode == "R":
        p = 1 #-1

    Y = np.sqrt(2) * Omega
    delta = (1/4) * cm.sqrt(1 - 8*Y**2)

    a = p * (1.0j / np.sqrt(2)) * (Y / (1 + Y**2))
    b = np.exp(-3*t/4) * (np.cosh(delta*t) + (3/(4*delta))*np.sinh(delta*t))
    c = p * 1.0j * np.sqrt(2) * Y * np.exp(-3*t/4) * (1/(4*delta)) * np.sinh(delta*t)

    return np.imag(a * (1 - b) + c)

L_list_analytical = []

for Omega in Omega_list:

    temp_list = [sigma_LR(t, "L", Omega) for t in time_list]
    L_list_analytical.append(temp_list)



### ----- Plotting -----------------------------------------------

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

    plt.plot(time_list, L_list[index][:,1], lw=5, c=colour_list[index], label=("$\Omega$ = " + str(Omega_list[index])))

    if index == 2:
        plt.plot(time_list, L_list_analytical[index], lw=5, c=colours.greek_dark_blue, ls="dotted", label="analytical")
    else:
        plt.plot(time_list, L_list_analytical[index], lw=5, c=colours.greek_dark_blue, ls="dotted")

plt.xlabel("Time $\gamma t$ (seconds)")
plt.ylabel("$\langle \sigma_+ (t) \\rangle$ & $\langle \sigma_- (t) \\rangle$")
plt.grid()
plt.legend(loc="lower right")
plt.xlim([0,13])
plt.savefig(directory + "sigma_L.pdf", facecolor=fig.get_facecolor(), transparent=True, dpi=600)

