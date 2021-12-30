"""

Plots only sigma L,R 

"""

import numpy as np 
import matplotlib.pyplot as plt 
import colours
import matplotlib
import cmath as cm

matplotlib.rcParams.update({'font.size': 24})

directory = "results2/"
Omega_list = [2.3, 1.0, 0.7]


### ----- Import data --------------------------------------
L_list = []
R_list = []

for Omega in Omega_list:

    name = directory + "sigma_L_" + (str(Omega).replace(".", "")) + ".txt"
    L_list.append(np.loadtxt(name)[0:13000])

    name2 = directory + "sigma_R_" + (str(Omega).replace(".", "")) + ".txt"
    R_list.append(np.loadtxt(name2)[0:13000])

time_list = L_list[0][:,0]


### ----- Analytical ----------------------------------------

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

L_list_analytical = []
R_list_analytical = []

for Omega in Omega_list:

    temp_list = [sigma_LR(t, "L", Omega) for t in time_list]
    L_list_analytical.append(temp_list)

    temp_list2 = [sigma_LR(t, "R", Omega) for t in time_list]
    R_list_analytical.append(temp_list2)



### ----- Plotting -----------------------------------------------

fig, ax = plt.subplots(2,1)
fig.set_size_inches(18.5, 10.5)
fig.set_alpha(0)
fig.set_facecolor("none")
fig.tight_layout(pad=3.0)
# ax = plt.axes()
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].spines['left'].set_color(colours.spanish_gray)
ax[0].spines['bottom'].set_color(colours.spanish_gray)
ax[0].tick_params(axis='x', colors=colours.spanish_gray)
ax[0].tick_params(axis='y', colors=colours.spanish_gray)
ax[0].xaxis.label.set_color(colours.spanish_gray)
ax[0].yaxis.label.set_color(colours.spanish_gray)
# ax[0].set_xlabel("$\gamma t$ (seconds)")
ax[0].set_ylabel("$\langle \sigma_- (t) \\rangle$")
ax[0].text(9.5, 0.4, "(a)", weight="bold")

colour_list = [colours.greek_dark_red, colours.orange_peel, colours.new_green]

for index in range(len(Omega_list)):

    ax[0].plot(time_list, L_list[index][:,1], lw=5, c=colour_list[index], label=("$\Omega$ = " + str(Omega_list[index])))

    if index == 2:
        ax[0].plot(time_list, L_list_analytical[index], lw=5, c=colours.greek_dark_blue, ls="dotted", label="analytical")
    else:
        ax[0].plot(time_list, L_list_analytical[index], lw=5, c=colours.greek_dark_blue, ls="dotted")


# plt.savefig(directory + "sigma_L.pdf", facecolor=fig.get_facecolor(), transparent=True, dpi=600)

ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['left'].set_color(colours.spanish_gray)
ax[1].spines['bottom'].set_color(colours.spanish_gray)
ax[1].tick_params(axis='x', colors=colours.spanish_gray)
ax[1].tick_params(axis='y', colors=colours.spanish_gray)
ax[1].xaxis.label.set_color(colours.spanish_gray)
ax[1].yaxis.label.set_color(colours.spanish_gray)
ax[1].set_xlabel("$\gamma t$ (seconds)")
ax[1].set_ylabel("$\langle \sigma_+ (t) \\rangle$")
ax[1].text(9.5, -0.05, "(b)", weight="bold")


for index in range(len(Omega_list)):

    ax[1].plot(time_list, R_list[index][:,1], lw=5, c=colour_list[index])
    ax[1].plot(time_list, R_list_analytical[index], lw=5, c=colours.greek_dark_blue, ls="dotted")

ax[0].legend()
# ax[1].legend()
# plt.legend(loc="lower right")
# plt.legend()
ax[0].set_xlim([0,10])
ax[1].set_xlim([0,10])
plt.savefig(directory + "sigma_LR.pdf", dpi=600)

# --- Plotting updated version 

fig = plt.figure()
fig.set_size_inches(18.5, 10.5)
fig.set_alpha(0)
fig.set_facecolor("none")
fig.tight_layout(pad=3.0)
ax = plt.axes()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_color(colours.spanish_gray)
ax.spines['bottom'].set_color(colours.spanish_gray)
ax.tick_params(axis='x', colors=colours.spanish_gray)
ax.tick_params(axis='y', colors=colours.spanish_gray)
ax.xaxis.label.set_color(colours.spanish_gray)
ax.yaxis.label.set_color(colours.spanish_gray)
ax.set_xlabel("$\gamma t$")
ax.set_ylabel("|$\langle \sigma_- (t) \\rangle$|")

colour_list = [colours.greek_dark_red, colours.orange_peel, colours.new_green]

for index in range(len(Omega_list)):

    ax.plot(time_list, L_list[index][:,1], lw=5, c=colour_list[index], label=("$\Omega$ = " + str(Omega_list[index])))

    if index == 2:
        ax.plot(time_list, L_list_analytical[index], lw=5, c=colours.greek_dark_blue, ls="dotted", label="analytical")
    else:
        ax.plot(time_list, L_list_analytical[index], lw=5, c=colours.greek_dark_blue, ls="dotted")

ax.set_xlim([0,10])
plt.legend()
plt.grid()
plt.savefig(directory + "sigma_LR_new.pdf", dpi=600)