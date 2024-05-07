import matplotlib.pyplot as plt
import numpy as np
import rebound

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
plt.rc('axes', labelsize=11) #fontsize of the x and y labels
plt.rcParams['text.latex.preamble']= r"\usepackage{amsmath}\usepackage{siunitx}\usepackage{newtx}"
plt.rcParams['errorbar.capsize']=4 # Einstellung Dicke der Fehlerbalkenenden
plt.rc('axes', labelsize=11, titlesize=16)
plt.rc('legend', fontsize=11)
plt.rc('xtick', labelsize=11) #fontsize of the x tick labels
plt.rc('ytick', labelsize=11) #fontsize of the y tick labels
plt.rcParams["figure.autolayout"] = True
plt.rcParams["axes.axisbelow"] = True
plt.rcParams["xtick.minor.visible"] = True
plt.rcParams["ytick.minor.visible"] = True
plt.rcParams["xtick.direction"] = 'in'
plt.rcParams["ytick.direction"] = 'in'
plt.rcParams["xtick.top"] = True
plt.rcParams["ytick.right"] = True



#%% setup simulation data

# Function for Hill-Radius
def r_hill(a, M_small_body, M_central_body):
    return a * (M_small_body / 3 * (M_central_body)) ** (1/3)

# Masses and Semi-major axis of Sun and Neptune
M_sun = 1
M_neptune = 5.151383772628957e-05
a_neptune = 30.237878368382898

sa = rebound.Simulationarchive("archives/archive-t_1e6-mercurius-10_rhill_02.bin")

max = len(sa)

#%% load initial results
sim_initial = sa[-max]

e_initial = []
a_initial = []
for i in range(2, 576):
    e_initial.append(sim_initial.particles[i].e)
    a_initial.append(sim_initial.particles[i].a)

e_initial = np.array(e_initial)
a_initial = np.array(a_initial)
q_initial = a_initial * (1-e_initial)

#%% load final results
sim_final = sa[-1]
print("sim-final: ", sim_final)
e_final = []
a_final = []
for i in range(2, 576):
    e_final.append(sim_final.particles[i].e)
    a_final.append(sim_final.particles[i].a)

e_final = np.array(e_final)
a_final = np.array(a_final)
q_final = a_final * (1-e_final)

#%% plotting
fig, axs = plt.subplots(1, 2, constrained_layout=True, sharey=True)

axs[0].scatter(a_initial, e_initial, s=5, label=r"$t_{\text{sim}}=0$")
axs[0].scatter(a_final, e_final, s=5, alpha=0.5, label=r"$t_{\text{sim}}=10^6\,yr$")
axs[0].set_xlabel(r"semi-major axis $a$ [au]")
axs[0].set_ylabel(r"eccentricity $e$")
axs[0].set_xscale("log")
axs[0].legend()
axs[0].set_ylim([0, 1])

axs[1].scatter(q_initial, e_initial, s=5)
axs[1].scatter(q_final, e_final, s=5, alpha=0.5)
axs[1].set_xlabel(r"pericentre distance $q$ [au]")
axs[1].set_ylim([0, 1])

plt.savefig("plots/mercurius-1e6-02.png")
plt.show()