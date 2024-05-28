import matplotlib.pyplot as plt
import numpy as np
import rebound

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
plt.rc('axes', labelsize=11)  #fontsize of the x and y labels
plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}\usepackage{siunitx}\usepackage{newtx}"
plt.rcParams['errorbar.capsize'] = 4  # Einstellung Dicke der Fehlerbalkenenden
plt.rc('axes', labelsize=11, titlesize=16)
plt.rc('legend', fontsize=11)
plt.rc('xtick', labelsize=11)  #fontsize of the x tick labels
plt.rc('ytick', labelsize=11)  #fontsize of the y tick labels
plt.rcParams["figure.autolayout"] = True
plt.rcParams["axes.axisbelow"] = True
plt.rcParams["xtick.minor.visible"] = True
plt.rcParams["ytick.minor.visible"] = True
plt.rcParams["xtick.direction"] = 'in'
plt.rcParams["ytick.direction"] = 'in'
plt.rcParams["xtick.top"] = True
plt.rcParams["ytick.right"] = True

# close all plot when opened
plt.close('all')


#%% setup simulation data

# Function for Hill-Radius
def r_hill(a, M_small_body, M_central_body):
    return a * (M_small_body / 3 * (M_central_body)) ** (1 / 3)


# Function for Parameter of neptune
def tisserand(T, a):
    return 0.5 * np.sqrt(-1 + (4 * (a ** 3)) + (2 * a * T) - ((a ** 2) * (T ** 2))) / (a ** 1.5)


# Masses and Semi-major axis of Sun and Neptune
M_sun = 1
M_neptune = 5.151383772628957e-05
a_neptune = 30.237878368382898

T_neptune = 1 / a_neptune + 2 * np.sqrt(a_neptune * (1 - 0.0113 ** 2))

sa1 = rebound.Simulationarchive("archives/mercurius-dtmin1-1e7-low.bin")
sa2 = rebound.Simulationarchive("archives/mercurius-dtmin1-1e7-high.bin")

max_1 = len(sa1)
max_2 = len(sa2)

print(sa1[-1])
print(sa2[-1])

#%% load initial results
sim1_initial = sa1[-max_1]
sim2_initial = sa2[-max_2]

# last particle index

i_last = sa1[1].N

e_initial1 = []
e_initial2 = []
#
a_initial1 = []
a_initial2 = []

for i in range(2, i_last):
    e_initial1.append(sim1_initial.particles[i].e)
    a_initial1.append(sim1_initial.particles[i].a)

for k in range(2, i_last):
    e_initial2.append(sim2_initial.particles[k].e)
    a_initial2.append(sim2_initial.particles[k].a)

e_initial1 = np.array(e_initial1)
a_initial1 = np.array(a_initial1)
q_initial1 = a_initial1 * (1 - e_initial1)

e_initial2 = np.array(e_initial2)
a_initial2 = np.array(a_initial2)
q_initial2 = a_initial2 * (1 - e_initial2)

# add arrays for both simulations
e_initial = np.append(e_initial1, e_initial2)
a_initial = np.append(a_initial1, a_initial2)
q_initial = a_initial * (1 - e_initial)

#%% load final results
sim_final1 = sa1[-1]
sim_final2 = sa2[-1]

i_last1 = sa1[-1].N
i_last2 = sa2[-1].N

e_final1 = []
a_final1 = []
e_final2 = []
a_final2 = []

for i in range(2, i_last1):
    e_final1.append(sim_final1.particles[i].e)
    a_final1.append(sim_final1.particles[i].a)

for k in range(2, i_last2):
    e_final2.append(sim_final2.particles[k].e)
    a_final2.append(sim_final2.particles[k].a)

e_final1 = np.array(e_final1)
a_final1 = np.array(a_final1)
q_final1 = a_final1 * (1 - e_final1)

e_final2 = np.array(e_final2)
a_final2 = np.array(a_final2)
q_final2 = a_final2 * (1 - e_final2)

# add arrays for both simulations
e_final = np.append(e_final1, e_final2)
a_final = np.append(a_final1, a_final2)
q_final = a_final * (1 - e_final)

#%% plotting
fig, axs = plt.subplots(1, 2, constrained_layout=True, sharey=True)
#
semis = np.arange(0, max(a_final), 1)
#
axs[0].scatter(a_initial, e_initial, s=5, label=r"$t_{\text{sim}}=0$", color="tab:blue")
axs[0].scatter(a_final, e_final, s=5, alpha=0.5, label=r"$t_{\text{sim}}=5\cdot 10^6\,yr$", color="tab:orange")
axs[0].plot(semis, tisserand(T_neptune, semis), label=r"$T_\text{Neptune}$", color="tab:red")
#
axs[0].set_xlabel(r"semi-major axis $a$ [au]")
axs[0].set_ylabel(r"eccentricity $e$")
axs[0].set_xscale("log")
axs[0].legend()
#axs[0].set_ylim([0, 1])
#
axs[1].scatter(q_initial, e_initial, s=5, color="tab:blue")
axs[1].scatter(q_final, e_final, s=5, alpha=0.5, color="tab:orange")
#
axs[1].set_xlabel(r"pericentre distance $q$ [au]")
#axs[1].set_ylim([0, 1])

#plt.savefig("plots/mercurius-1e6-02.png")
plt.show()

'''#%% 2d Histogram
fig, axs = plt.subplots(1, 2, constrained_layout=True, sharey=True)
axs[0].set_xlim(min(q_final), max(q_final))
axs[0].hist2d(q_initial, e_initial, bins=(100, 100), cmap="inferno")
axs[1].hist2d(q_final, e_final, bins=(100, 100), cmap="inferno")
#
axs[0].set_xlabel(r"pericentre distance $q$ [au]")
axs[1].set_xlabel(r"pericentre distance $q$ [au]")
axs[0].set_ylabel(r"eccentricity $e$")
plt.show()'''

#%% using density scatter plot
# from https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density/53865762#53865762https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density/53865762#53865762
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy.interpolate import interpn


#@numba.jit
def density_scatter(x, y, ax=None, sort=True, bins=20):
    """
    Scatter plot colored by 2d histogram
    """
    if ax is None:
        fig, ax = plt.subplots(constrained_layout=True)
    data, x_e, y_e = np.histogram2d(x, y, bins=bins, density=True)
    z = interpn((0.5 * (x_e[1:] + x_e[:-1]), 0.5 * (y_e[1:] + y_e[:-1])), data, np.vstack([x, y]).T, method="splinef2d",
                bounds_error=False)

    #To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0

    # Sort the points by density, so that the densest points are plotted last
    if sort:
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    ax.scatter(x, y, c=z)
    ax.set_xlabel(r"pericentre distance $q$ [au]")
    ax.set_ylabel(r"eccentricity $e$")

    norm = Normalize(vmin=np.min(z), vmax=np.max(z))
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm), ax=ax)
    cbar.ax.set_ylabel('Density')

    return ax


#density_scatter(q_final, e_final, bins=[100, 100])
#plt.show()

#density_scatter(q_initial, e_initial, bins=[100, 100])
#plt.show()

#%%
from scipy.optimize import curve_fit

def fit_func(x, a, b, c, d):
    return a * np.exp(-b * (x/c) ** d)

'''def fit_func(x, a, b):
    return 1 - x * a * b'''


xdata = []
ydata = []

for ind in range(len(q_final)):
    if q_final[ind] > 30:
        xdata.append(q_final[ind])
        ydata.append(abs(e_final[ind] - e_initial[ind]))

xdata = np.array(xdata)
ydata = np.array(ydata)

p0 = [999999, 10, 25, 2]

popt, pcov = curve_fit(fit_func, xdata, ydata, p0=p0)

fig, ax = plt.subplots(constrained_layout=True)
ax.scatter(xdata, ydata, s=5, color="tab:blue")
#ax.plot(np.sort(xdata), fit_func(np.sort(xdata), *popt), color="tab:orange")
ax.set_xlabel(r"pericentre distance $q$ [au]")
ax.set_ylabel(r"eccentricity deviation $\Delta e$")
#ax.set_xscale("log")
#ax.set_yscale("log")
#ax.set_xlim(30, max(q_final))
plt.show()
