"""
Plot the g2 plot
"""

import numpy as np 
import matplotlib.pyplot as plt 
import cmath as cm 
import colours
import matplotlib
matplotlib.rcParams.update({'font.size': 24})

data = np.loadtxt("results2/g2_23.txt")
data2 = np.loadtxt("results2/denom.txt")
data3 = np.loadtxt("results2/z_norm.txt")
data4 = np.loadtxt("results2/avg_e.txt")

time_list = data[:,0]
G2_list = data[:,1] 
divider = data2[:,1]
z_list = data3[:,1]
avg_e = data4


"""
ratio notes:

end_time: doesn't matter 
simulation: increasing simulation by 2x increases ration by 2x (simulation count has to matter)
            This does also 2x the number of total emissions 
timesteps: doesn't matter 

intensity seems to also matter 
2.3: ratio 5223, emission 339022
1.0: ratio 2773, emission 242045
0.7: ratio 1527, emission 176087

"""

end_time = 30
num_of_simulations = 10000 
emission = 136029
norm = sum(avg_e) / np.size(avg_e)
print(norm)

z_list = [(0.5*(1 + z_list[i]))**2 for i in range(len(z_list))]
# norm = sum(z_list) / np.size(z_list)

def g2(t, Omega):

    Y = np.sqrt(2) * Omega
    delta = (1/4) * cm.sqrt(1 - 8*Y**2)

    a = np.exp(-3*t/4) * (np.cosh(delta * t) + ((3/(4*delta)) * np.sinh(delta*t)))

    return 1 - a

g2_analytical = [g2(t, 2.3) for t in time_list]

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

G2_list = [G2_list[i] / norm**2 for i in range(len(G2_list))]
print(max(G2_list)/max(g2_analytical))

plt.plot(time_list, G2_list, label="numerical", lw=5, c=colours.orange_peel)
plt.plot(time_list, g2_analytical, ls="dotted", label="analytical", lw=5, c=colours.greek_dark_blue)
plt.xlabel("Time (seconds)")
plt.ylabel("$g^{(2)}(t)$")
plt.xlim([0,12])
plt.legend()
plt.grid()
plt.savefig("results2/g2.pdf", dpi=600)