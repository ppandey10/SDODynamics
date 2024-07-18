import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import rebound

plt.style.use("~/custom.mplstyle")
# Plot the filtered data
sns.set_style("whitegrid", {'axes.grid' : False, 'axes.edgecolor': 'black'})


# Load the simulation data
# sim = rebound.Simulationarchive("archives/mercurius_tsim_1e7_dtias15_5e-4_dt_5_rhill_2-no_inc.bin")
sim = rebound.Simulationarchive("archives/inc-final-1e7.bin")
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



fig, axs = plt.subplots(1,3, sharey=True, figsize=(16,9))
plt.subplots_adjust()
axs[0].scatter(filtered_p_values_first, proper_elements_initial[:,1])
axs[1].scatter(proper_elements_initial[:,0], proper_elements_initial[:,1])
axs[2].scatter(proper_elements_initial[:,2], proper_elements_initial[:,1])

# Create the joint plot
# Add KDE plot
joint = sns.jointplot(
    x=filtered_p_values_last, y=filtered_eccentricity_values_last, 
    kind="scatter", color='red', marginal_kws=dict(bins=50), height=6, s=8,
    marker='x', alpha=0.9
)
kdeplot = sns.kdeplot(
    x=filtered_p_values_last, 
    y=filtered_eccentricity_values_last,
    # bw_adjust=1.5, 
    thresh=0, 
    levels=10,
    cmap='coolwarm', fill=True, alpha=0.5, 
    cbar=True, cbar_kws={'shrink': 0.9, 'aspect':40, 'label': 'KDE density'},
    clip=((min(filtered_p_values_last), max(filtered_p_values_last)), (-0.1, 1.2)),  # Adjusted limits for the KDE plot
)

plt.subplots_adjust(left=0.1, right=0.8, top=0.99, bottom=0.1)
pos_joint_ax = joint.ax_joint.get_position()
pos_marg_x_ax = joint.ax_marg_x.get_position()
joint.ax_joint.set_position([pos_joint_ax.x0, pos_joint_ax.y0, pos_marg_x_ax.width, pos_joint_ax.height])
joint.fig.axes[-1].set_position([.83, pos_joint_ax.y0, .07, pos_joint_ax.height])

# plt.axvline(x=4.5, color='red', linestyle='--')  # Vertical line at x=25
# plt.axvline(x=5.1, color='red', linestyle='--')  # Vertical line at x=31
plt.xlabel(r'Perihelion Distance $q$ [AU]')
plt.ylabel(r'Eccentricity $e$')
#plt.suptitle(f'Time Step: {time_step_index}', y=0.95)
# plt.xlim(16.5, 40)
# plt.ylim(0,0.93)
# plt.tight_layout()
# plt.savefig('try2.png')


## scatterplot with final and last timesteps


plt.show()
