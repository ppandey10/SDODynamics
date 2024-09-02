import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import rebound

plt.style.use("~/custom.mplstyle")
# Plot the filtered data
#sns.set_style("whitegrid", {'axes.grid' : False, 'axes.edgecolor': 'black'})


# Load the simulation data
# sim = rebound.Simulationarchive("archives/mercurius_tsim_1e7_dtias15_5e-4_dt_5_rhill_2-no_inc.bin")
sim = rebound.Simulationarchive("archives/inc-final-1e8.bin")
# sim = rebound.Simulationarchive("archives/no_inc-final_1e7.bin")


initial_sim = sim[0]  # First snapshot
final_sim = sim[len(sim)-1]   # Last snapshot (10 Myr later)

print(initial_sim)
print(final_sim)

# Extract the proper elements
particles_initial = initial_sim.particles
particles_final = final_sim.particles

# Initialize lists to store proper elements
proper_elements_initial = []
proper_elements_final = []

# Iterate through particles to extract (a, e) at initial and final times
for i in range(2, len(particles_initial)):  # first particle is the Sun
    p_initial = particles_initial[i]
    p_final = particles_final[i]
    
    a_initial, e_initial, i_initial = p_initial.a, p_initial.e, p_initial.inc
    a_final, e_final, i_final = p_final.a, p_final.e, p_initial.inc
    
    proper_elements_initial.append((a_initial, e_initial, i_initial))
    proper_elements_final.append((a_final, e_final, i_final))

proper_elements_initial = np.array(proper_elements_initial)
proper_elements_final = np.array(proper_elements_final)

# Remove particles with e >= 1 and count how many were removed
mask = proper_elements_final[:, 1] <= 1
removed_particles_count = len(proper_elements_final) - np.sum(mask)
proper_elements_final = proper_elements_final[mask]
proper_elements_initial = proper_elements_initial[mask]

print(f"Number of particles removed: {removed_particles_count}")

filtered_p_values_first = proper_elements_initial[:,0] * (1-proper_elements_initial[:,1])
filtered_p_values_last = proper_elements_final[:,0] * (1-proper_elements_final[:,1])
filtered_eccentricity_values_last = proper_elements_final[:,1]

# add lines for tisserand paramter and resonances

# Masses and Semi-major axis of Sun and Neptune
M_sun = 1
M_neptune = 5.151383772628957e-05
a_neptune = 30.237878368382898

T_neptune = 1 / a_neptune + 2 * np.sqrt(a_neptune * (1 - 0.0113 ** 2))

# Function for Hill-Radius
def r_hill(a, M_small_body, M_central_body):
    return a * (M_small_body / 3 * (M_central_body)) ** (1 / 3)


# Function for Parameter of neptune
def tisserand(T, a):
    return 0.5 * np.sqrt(-1 + (4 * (a ** 3)) + (2 * a * T) - ((a ** 2) * (T ** 2))) / (a ** 1.5)

def res_semi(j, l, a1):
    return ( ((j+l)/j)**2 * a1**3) ** (1/3)

def res_peri(j, l, a1, q):
    return 1 - (q / (( ((j+l)/j)**2 * a1**3) ** (1/3)))

semis = np.arange(0, max(proper_elements_final[:,0]), 1)



fig, axs = plt.subplots(1,3, sharey=True, figsize=(12,3.4))
axs[0].scatter(filtered_p_values_first, proper_elements_initial[:,1], s=2)
axs[0].scatter(filtered_p_values_last, proper_elements_final[:,1], s=2)
axs[0].plot(np.linspace(5,40,100), res_peri(2, 1, a_neptune, np.linspace(5,40,100)), color='tab:grey', alpha=0.7, lw=1)
axs[0].plot(np.linspace(5,45,100), res_peri(1, 1, a_neptune, np.linspace(5,45,100)), color='tab:grey', alpha=0.7, lw=1)
axs[0].text(11, 0.6, "3:2", rotation=-35, color='tab:grey', alpha=0.7)
axs[0].text(11, 0.75, "2:1", rotation=-35, color='tab:grey', alpha=0.7)
axs[0].set_xlabel(r"pericentre distance $q$ [AU]")
axs[0].set_ylabel(r"eccentricity $e$")


axs[1].scatter(proper_elements_initial[:,0], proper_elements_initial[:,1], s=2)
axs[1].scatter(proper_elements_final[:,0], proper_elements_final[:,1], s=2)
axs[1].set_xscale("log")
axs[1].set_xlabel(r"semi-major axis $a$ [AU]")
# axs[1].plot(semis, tisserand(T_neptune, semis))
axs[1].axvline(res_semi(2, 1, a_neptune), 0, 1, color='tab:grey', alpha=0.7, lw=1)
axs[1].axvline(res_semi(1, 1, a_neptune), 0, 1, color='tab:grey', alpha=0.7, lw=1)
axs[1].text(30, 0.9, "3:2", rotation='vertical', color='tab:grey', alpha=0.7)
axs[1].text(49, 0.9, "2:1", rotation='vertical', color='tab:grey', alpha=0.7)
# Inset
inset_ax = axs[1].inset_axes([0.5, 0.1, 0.45, 0.6], transform=axs[1].transAxes)
inset_ax.scatter(proper_elements_initial[:, 0], proper_elements_initial[:, 1], s=2)
inset_ax.scatter(proper_elements_final[:, 0], proper_elements_final[:, 1], s=2)
inset_ax.axvline(res_semi(2, 1, a_neptune), 0, 1, color='tab:grey', alpha=0.7, lw=1)
inset_ax.axvline(res_semi(1, 1, a_neptune), 0, 1, color='tab:grey', alpha=0.7, lw=1)
inset_ax.set_xlim([38, 50])
inset_ax.set_ylim([0, 0.3])
# inset_ax.set_xticklabels([])
# inset_ax.set_yticklabels([])
inset_ax.set_xticks([40, 45, 50])
inset_ax.set_yticks([0, 0.1, 0.2, 0.3])


axs[2].scatter(proper_elements_initial[:,2] * 180/np.pi , proper_elements_initial[:,1], s=2)
axs[2].scatter(proper_elements_final[:,2]* 180/np.pi , proper_elements_final[:,1], s=2)
axs[2].set_xlabel(r"inclination $i$ [°]")

plt.savefig("plots/scatterplot-all.pdf")

plt.show()
