"""

Plots results for Little boxes simple 

"""

import numpy as np 
import matplotlib.pyplot as plt 
import colours
import matplotlib
import cmath as cm 
matplotlib.rcParams.update({'font.size': 24})

directory = "results/"

### ----- Sort data ------------------------------------------------------------
z_list = np.loadtxt("sigma_z.txt")
L_list = np.loadtxt("sigma_L.txt")
R_list = np.loadtxt("sigma_R.txt")
g2_list = np.loadtxt("g2.txt")

time_list = z_list[:,0]
z_list = z_list[:,1]
L_list = L_list[:,1]
R_list = R_list[:,1]
g2_list = g2_list[:,1]

### ----- Analytical results ---------------------------------------------------

Omega = 0.025 # NOTE: the omega here is Omega / gamma 

def sigma_z(t):

    Y = np.sqrt(2) * Omega

    delta = (1/4) * cm.sqrt(1 - 8*Y**2)

    a = -1 * (1/(1 + Y**2))
    b = Y**2 * np.exp(-3*t/4) 
    c = np.cosh(delta*t) + (3/(4*delta))*np.sinh(delta*t)

    return np.real(a * (1 + b*c))

def sigma_LR(t, mode):

    p = 1

    if mode == "R":
        p = 1 #-1

    Y = np.sqrt(2) * Omega
    delta = (1/4) * cm.sqrt(1 - 8*Y**2)

    a = p * (1.0j / np.sqrt(2)) * (Y / (1 + Y**2))
    b = np.exp(-3*t/4) * (np.cosh(delta*t) + (3/(4*delta))*np.sinh(delta*t))
    c = p * 1.0j * np.sqrt(2) * Y * np.exp(-3*t/4) * (1/(4*delta)) * np.sinh(delta*t)

    return np.imag(a * (1 - b) + c)

def g2(t):

    Y = np.sqrt(2) * Omega
    delta = (1/4) * cm.sqrt(1 - 8*Y**2)

    a = np.exp(-3*t/4) * (np.cosh(delta * t) + ((3/(4*delta)) * np.sinh(delta*t)))

    return 1 - a


z_list_A = [sigma_z(t) for t in time_list]
L_list_A = [sigma_LR(t,"L") for t in time_list]
R_list_A = [sigma_LR(t,"R") for t in time_list]
g2_list_A = [g2(t) for t in time_list]

### ----- Plotting -------------------------------------------------------------

# For sigma_z 
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

plt.plot(time_list, z_list, linewidth=5, c=colours.greek_red, label="Numerical")
plt.plot(time_list, z_list_A, linewidth=5, c=colours.greek_blue, linestyle="dotted", label="Analytical")
plt.xlabel("Time (seconds)")
plt.ylabel("$\langle \sigma_z (t) \\rangle$")
plt.grid()
plt.legend()
plt.savefig(directory + "sigma_z.pdf", facecolor=fig1.get_facecolor(), transparent=True, dpi=600)

# For sigma_- 
fig2 = plt.figure(2)
fig2.set_size_inches(18.5, 10.5)
fig2.set_alpha(0)
fig2.set_facecolor("none")
ax2 = plt.axes()
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['left'].set_color(colours.spanish_gray)
ax2.spines['bottom'].set_color(colours.spanish_gray)
ax2.tick_params(axis='x', colors=colours.spanish_gray)
ax2.tick_params(axis='y', colors=colours.spanish_gray)
ax2.xaxis.label.set_color(colours.spanish_gray)
ax2.yaxis.label.set_color(colours.spanish_gray)

plt.plot(time_list, L_list, linewidth=5, c=colours.greek_red, label="Numerical")
plt.plot(time_list, L_list_A, linewidth=5, c=colours.greek_blue, linestyle="dotted", label="Analytical")
plt.xlabel("Time (seconds)")
plt.ylabel("$\langle \sigma_- (t) \\rangle$")
plt.grid()
plt.legend()
plt.savefig(directory + "sigma_L.pdf", facecolor=fig1.get_facecolor(), transparent=True, dpi=600)

# For sigma_+ 
fig3 = plt.figure(3)
fig3.set_size_inches(18.5, 10.5)
fig3.set_alpha(0)
fig3.set_facecolor("none")
ax3 = plt.axes()
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['left'].set_color(colours.spanish_gray)
ax3.spines['bottom'].set_color(colours.spanish_gray)
ax3.tick_params(axis='x', colors=colours.spanish_gray)
ax3.tick_params(axis='y', colors=colours.spanish_gray)
ax3.xaxis.label.set_color(colours.spanish_gray)
ax3.yaxis.label.set_color(colours.spanish_gray)

plt.plot(time_list, R_list, linewidth=5, c=colours.greek_red, label="Numerical")
plt.plot(time_list, R_list_A, linewidth=5, c=colours.greek_blue, linestyle="dotted", label="Analytical")
plt.xlabel("Time (seconds)")
plt.ylabel("$\langle \sigma_z (t) \\rangle$")
plt.grid()
plt.legend()
plt.savefig(directory + "sigma_R.pdf", facecolor=fig1.get_facecolor(), transparent=True, dpi=600)

# Plot g2 

fig4 = plt.figure(4)
fig4.set_size_inches(18.5, 10.5)
fig4.set_alpha(0)
fig4.set_facecolor("none")
ax4 = plt.axes()
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.spines['left'].set_color(colours.spanish_gray)
ax4.spines['bottom'].set_color(colours.spanish_gray)
ax4.tick_params(axis='x', colors=colours.spanish_gray)
ax4.tick_params(axis='y', colors=colours.spanish_gray)
ax4.xaxis.label.set_color(colours.spanish_gray)
ax4.yaxis.label.set_color(colours.spanish_gray)

plt.plot(time_list, g2_list, linewidth=5, c=colours.greek_red, label="g2")
plt.plot(time_list, g2_list_A, linewidth=5, c=colours.greek_blue, label="analytical")
plt.xlabel("Time (seconds)")
plt.ylabel("g2")
plt.grid()
plt.legend()
plt.savefig(directory + "g2.pdf", facecolor=fig1.get_facecolor(), transparent=True, dpi=600)