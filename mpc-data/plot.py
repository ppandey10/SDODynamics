import pandas as pd
import matplotlib.pyplot as plt
import rebound
import numpy as np
plt.style.use("~/custom.mplstyle")

# Define the path to the text file
file_path = 'mpc-data-scattering-centaurs.txt'

# Read the text file, skipping the first line and using fixed-width format
column_names = [
    "Designation", "Prov", "Des", "q", "Q", "H", "Epoch", "M", "Peri", 
    "Node", "Incl", "e", "a", "Opps", "Ref", "Designation_name", "Discovery_date_site_discoverer"
]

# Define the fixed widths for each column based on the structure of the file
column_widths = [
    27, 5, 7, 7, 10, 7, 10, 7, 6, 6, 6, 6, 9, 9, 19, 33, 48
]

# Read the file using pandas.read_fwf
df = pd.read_fwf(file_path, widths=column_widths, names=column_names, skiprows=1)

# import data from simulation

sim = rebound.Simulationarchive("../particle-cloud/archives/mercurius_tsim_1e7_dtias15_5e-4_dt_5_rhill_2_final.bin")
initial_sim = sim[0]  # First snapshot
final_sim = sim[len(sim)-2]   # Last snapshot (10 Myr later)

print(initial_sim)
print(final_sim)

# Extract the proper elements
particles_initial = initial_sim.particles
particles_final = final_sim.particles

# Initialize lists to store proper elements
proper_elements_initial = []
proper_elements_final = []

# Iterate through particles to extract (a, e) at initial and final times
for i in range(1, len(particles_initial)):  # first particle is the Sun
    p_initial = particles_initial[i]
    p_final = particles_final[i]
    
    a_initial, e_initial = p_initial.a, p_initial.e
    a_final, e_final = p_final.a, p_final.e
    
    proper_elements_initial.append((a_initial, e_initial))
    proper_elements_final.append((a_final, e_final))

proper_elements_initial = np.array(proper_elements_initial)
proper_elements_final = np.array(proper_elements_final)

# Remove particles with e >= 1 and count how many were removed
mask = proper_elements_final[:, 1] <= 1
removed_particles_count = len(proper_elements_final) - np.sum(mask)
proper_elements_final = proper_elements_final[mask]
proper_elements_initial = proper_elements_initial[mask]

a_sim = proper_elements_final[:,0]
e_sim = proper_elements_final[:,1]

fig, ax = plt.subplots()
ax.scatter(df["a"], df["e"], s=3, color="grey", alpha=0.3, label=r"SDOs from MPC")
ax.scatter(a_sim, e_sim, s=3, label=r"$10\,\text{Myr}$ simulation")
ax.set_xscale("log")
ax.set_xlabel(r"semi-major axis $a$ [au]")
ax.set_ylabel(r"eccentricity $e$")
ax.legend(markerscale=3)
plt.savefig("../particle-cloud/plots/mpc-comp.pdf")
plt.show()

# TODO: remove centaurs, calculate perihelion for 3:2, 2:1 resonances
