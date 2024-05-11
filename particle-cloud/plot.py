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

# close all plot when opened
plt.close('all')

#%% setup simulation data

# Function for Hill-Radius
def r_hill(a, M_small_body, M_central_body):
    return a * (M_small_body / 3 * (M_central_body)) ** (1/3)

# Masses and Semi-major axis of Sun and Neptune
M_sun = 1
M_neptune = 5.151383772628957e-05
a_neptune = 30.237878368382898

sa1 = rebound.Simulationarchive("archives/custom-10e6-e-high-04.bin")
sa2 = rebound.Simulationarchive("archives/custom-10e6-e-low-04.bin")

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

for k in range(2,i_last):
    e_initial2.append(sim2_initial.particles[k].e)
    a_initial2.append(sim2_initial.particles[k].a)

e_initial1 = np.array(e_initial1)
a_initial1 = np.array(a_initial1)
q_initial1 = a_initial1 * (1-e_initial1)

e_initial2 = np.array(e_initial2)
a_initial2 = np.array(a_initial2)
q_initial2 = a_initial2 * (1-e_initial2)

#%% load final results
sim_final1 = sa1[-1]
sim_final2 = sa2[-1]


e_final1 = []
a_final1 = []
e_final2 = []
a_final2 = []

for i in range(2, i_last):
    e_final1.append(sim_final1.particles[i].e)
    a_final1.append(sim_final1.particles[i].a)

for k in range(2, i_last):
    e_final2.append(sim_final2.particles[k].e)
    a_final2.append(sim_final2.particles[k].a)


e_final1 = np.array(e_final1)
a_final1 = np.array(a_final1)
q_final1 = a_final1 * (1-e_final1)

e_final2 = np.array(e_final2)
a_final2 = np.array(a_final2)
q_final2 = a_final2 * (1-e_final2)

#%% plotting
fig, axs = plt.subplots(1, 2, constrained_layout=True, sharey=True)
#
axs[0].scatter(a_initial1, e_initial1, s=5, label=r"$t_{\text{sim}}=0$", color = "tab:blue")
axs[0].scatter(a_initial2, e_initial2, s=5,color = "tab:blue")
axs[0].scatter(a_final1, e_final1, s=5, alpha=0.5, label=r"$t_{\text{sim}}=10^6\,yr$", color = "tab:orange")
axs[0].scatter(a_final2, e_final2, s=5, alpha=0.5, color = "tab:orange")
#
axs[0].set_xlabel(r"semi-major axis $a$ [au]")
axs[0].set_ylabel(r"eccentricity $e$")
axs[0].set_xscale("log")
axs[0].legend()
axs[0].set_ylim([0, 1])
#
axs[1].scatter(q_initial1, e_initial1, s=5, color = "tab:blue")
axs[1].scatter(q_initial2, e_initial2, s=5, color = "tab:blue")
axs[1].scatter(q_final1, e_final1, s=5, alpha=0.5, color = "tab:orange")
axs[1].scatter(q_final2, e_final2, s=5, alpha=0.5, color = "tab:orange")
#
axs[1].set_xlabel(r"pericentre distance $q$ [au]")
axs[1].set_ylim([0, 1])

#plt.savefig("plots/mercurius-1e6-02.png")
plt.show()