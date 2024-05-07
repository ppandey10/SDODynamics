import matplotlib.pyplot as plt
import numpy as np
import rebound

#%% setup simulation data


# Function for Hill-Radius
def r_hill(a, M_small_body, M_central_body):
    return a * (M_small_body / 3 * (M_central_body)) ** (1/3)

# Masses and Semi-major axis of Sun and Neptune
M_sun = 1
M_neptune = 5.151383772628957e-05
a_neptune = 30.237878368382898

sa = rebound.Simulationarchive("archives/archive-t_1e6-mercurius-10_rhill_04.bin")

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

e_final = []
a_final = []
for i in range(2, 576):
    e_final.append(sim_final.particles[i].e)
    a_final.append(sim_final.particles[i].a)

e_final = np.array(e_final)
a_final = np.array(a_final)
q_final = a_final * (1-e_final)

#%% plotting
fig, axs = plt.subplots(1,2,constrained_layout=True,sharey=True)
axs[0].scatter(a_initial, e_initial, s=10, label="initial distribution")
axs[0].scatter(a_final, e_final, s=10, label="final distribution")
axs[0].legend()
axs[0].set_ylim([0, 1])

axs[1].scatter(q_initial, e_initial, s=10, label="initial distribution")
axs[1].scatter(q_final, e_final, s=10, label="final distribution")
axs[1].set_ylim([0, 1])
plt.show()