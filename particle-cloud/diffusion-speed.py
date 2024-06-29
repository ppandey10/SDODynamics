import rebound
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm, LogNorm

plt.style.use("~/custom.mplstyle")

# Load the simulation data
sim = rebound.Simulationarchive("archives/mercurius_tsim_1e7_dtias15_5e-4_dt_5_rhill_2_final.bin")
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

# for i in range(len(proper_elements_initial)):
#     print("(", proper_elements_initial[i,0], proper_elements_initial[i,1], ")", 
#           "(", proper_elements_final[i,0], proper_elements_final[i,1], ")") 

print(f"Number of particles removed: {removed_particles_count}")

# Define the phase space grid
bin_num = 20
a_bins = np.linspace(proper_elements_initial[:, 0].min(), proper_elements_initial[:, 0].max(), num=bin_num)
e_bins = np.linspace(proper_elements_initial[:, 1].min(), proper_elements_initial[:, 1].max(), num=bin_num)

# Create a dictionary to store the differences for each box
diffs = {(i, j): [] for i in range(len(a_bins)-1) for j in range(len(e_bins)-1)}

# Function to get the box index for a given (a, e)
def get_box_index(a, e):
    i = np.digitize(a, a_bins) - 1
    j = np.digitize(e, e_bins) - 1
    return i, j

# Calculate differences and store them in the appropriate bins
for (a_initial, e_initial), (a_final, e_final) in zip(proper_elements_initial, proper_elements_final):
    delta_a = a_final - a_initial
    delta_e = e_final - e_initial
    i, j = get_box_index(a_initial, e_initial)
    
    # Store the differences in the corresponding box
    if (i, j) in diffs:
        diffs[(i, j)].append((delta_a, delta_e))

# print("diffs", diffs)

# Compute the average differences for each box
diffusion_speeds = {}
for key, values in diffs.items():
    if values:
        avg_delta_a = np.mean([v[0] for v in values])
        avg_delta_e = np.mean([v[1] for v in values])
        diffusion_speeds[key] = (avg_delta_a, avg_delta_e)

# Convert the result to a DataFrame 
result = []
for (i, j), (avg_delta_a, avg_delta_e) in diffusion_speeds.items():
    result.append({'a_bin': (a_bins[i], a_bins[i+1]), 'e_bin': (e_bins[j], e_bins[j+1]), 
                   'avg_delta_a': avg_delta_a, 'avg_delta_e': avg_delta_e})

result_df = pd.DataFrame(result)

# Save the result to a CSV file
result_df.to_csv('diffusion_speeds.csv', index=False)

# Extract the mid-points of the bins for plotting
a_centers = (a_bins[:-1] + a_bins[1:]) / 2
e_centers = (e_bins[:-1] + e_bins[1:]) / 2

# Create 2D arrays for the mid-points and the average delta_a
a_grid, e_grid = np.meshgrid(a_centers, e_centers)
avg_delta_a_grid = np.full_like(a_grid, np.nan, dtype=np.double)
avg_delta_e_grid = np.full_like(e_grid, np.nan, dtype=np.double)

for (i, j), (avg_delta_a, avg_delta_e) in diffusion_speeds.items():
    avg_delta_a_grid[j, i] = avg_delta_a
    avg_delta_e_grid[j, i] = avg_delta_e

# calculate Tisserand-curve
a_N = 30.07
T = -1.5
a_values = np.linspace(a_bins.min(), a_bins.max(), 500)
e_values = np.sqrt(abs(1 - (a_N/(4*a_values) * (T - a_N/a_values)**2) ))

# calculate curves for fixed pericentre distances
def peri_curve(q):
    e = 1-q/a_values
    return e

### Here starts the plotting ##########################################################################################

# Scatterplot with distribution of particles
fig, ax = plt.subplots()
ax.scatter(proper_elements_initial[:,0], proper_elements_initial[:,1], s=3, label=r"$t_\text{sim}=0\,\text{Myr}$")
ax.scatter(proper_elements_final[:,0], proper_elements_final[:,1], s=3, label=r"$t_\text{sim}=10\,\text{Myr}$")
# ax.plot(a_values, e_values, label=r"$q=40\,\text{au}$")
ax.plot(a_values, peri_curve(40), color="k",label=r"$q=40\,\text{au}$")
ax.set_ylabel(r"eccentricity e")
ax.set_xlabel(r"semi-major axis a [au]")
ax.set_ylim(0,1.05)
ax.set_xscale("log")
ax.legend()
# plt.savefig("plots/10-Myr-comp.pdf")

### Plot the 2D histogram

# normalization
sym_norm = SymLogNorm(linthresh=0.03, linscale=0.03, vmin=-1000, vmax=1000, base=10)

## for semi-major-axis

fig, ax = plt.subplots()

c = ax.pcolormesh(a_grid, e_grid, avg_delta_a_grid, shading='auto', cmap='coolwarm',norm=sym_norm)
cbar = fig.colorbar(c, ax=ax, label=r"$\Delta a$ [au/$10\,\text{Myr}$]" )
cbar.set_ticks(ticks=[-1000, -100, -10, -1, -0.1, -0, 0.1, 1, 10, 100, 1000])

# ax.plot(a_values, peri_curve(30), color ="k", label=r"$q=30\,\text{au}$", linewidth=1)
# ax.plot(a_values, peri_curve(35), color ="k", label=r"$q=35\,\text{au}$", linewidth=1)
# ax.plot(a_values, peri_curve(40), color ="k", label=r"$q=40\,\text{au}$", linewidth=1)
# ax.plot(a_values, peri_curve(45), color ="k", label=r"$q=45\,\text{au}$", linewidth=1)

ax.set_xlabel(r"semi-major axis $a$ [au]")
ax.set_ylabel(r"eccentricity $e$")

ax.set_ylim(0,1.05)
ax.legend()

plt.savefig("plots/Delta_a-10Myr.pdf")

## for eccentricities

fig, ax = plt.subplots()

# Find the maximum absolute value in avg_delta_e_grid
max_abs = max(abs(avg_delta_e_grid.min()), abs(avg_delta_e_grid.max()))

c = ax.pcolormesh(a_grid, e_grid, avg_delta_e_grid, shading='auto', cmap='coolwarm', vmin=-max_abs, vmax=max_abs)
cbar = fig.colorbar(c, ax=ax, label=r"$\Delta e$ [au/$10\,\text{Myr}$]" )

# ax.plot(a_values, peri_curve(30), color ="k", label=r"$q=30\,\text{au}$", linewidth=1)
# ax.plot(a_values, peri_curve(35), color ="k", label=r"$q=35\,\text{au}$", linewidth=1)
# ax.plot(a_values, peri_curve(40), color ="k", label=r"$q=40\,\text{au}$", linewidth=1)
# ax.plot(a_values, peri_curve(45), color ="k", label=r"$q=45\,\text{au}$", linewidth=1)

ax.set_xlabel(r"semi-major axis $a$ [au]")
ax.set_ylabel(r"eccentricity $e$")

ax.set_ylim(0,1.05)
ax.legend() 
plt.savefig("plots/Delta_e-10Myr.pdf")


## absolute change of orbital elements 

# for semi-major axis

fig, ax = plt.subplots()

norm = LogNorm(vmin=0.01, vmax=1000) 

c = ax.pcolormesh(a_grid, e_grid,abs(avg_delta_a_grid), shading='auto', cmap='Reds',norm=norm)
cbar = fig.colorbar(c, ax=ax, label=r"$|\Delta a|$ [au/$10\,\text{Myr}$]" )

# ax.plot(a_values, peri_curve(30), color ="k", label=r"$q=30\,\text{au}$", linewidth=1)
# ax.plot(a_values, peri_curve(35), color ="k", label=r"$q=35\,\text{au}$", linewidth=1)
# ax.plot(a_values, peri_curve(40), color ="k", label=r"$q=40\,\text{au}$", linewidth=1)
# ax.plot(a_values, peri_curve(45), color ="k", label=r"$q=45\,\text{au}$", linewidth=1)

ax.set_xlabel(r"semi-major axis $a$ [au]")
ax.set_ylabel(r"eccentricity $e$")

ax.set_ylim(0,1.05)
ax.legend() 
plt.savefig("plots/abs-Delta_a-10Myr.pdf")


# for eccentrity

fig, ax = plt.subplots()

# Find the maximum absolute value in avg_delta_e_grid
max_abs = max(abs(avg_delta_e_grid.min()), abs(avg_delta_e_grid.max()))

c = ax.pcolormesh(a_grid, e_grid, abs(avg_delta_e_grid), shading='auto', cmap='Reds') #, vmin=-max_abs, vmax=max_abs)
cbar = fig.colorbar(c, ax=ax, label=r"$|\Delta e|$ [au/$10\,\text{Myr}$]" )

# ax.plot(a_values, peri_curve(30), color ="k", label=r"$q=30\,\text{au}$", linewidth=1)
# ax.plot(a_values, peri_curve(35), color ="k", label=r"$q=35\,\text{au}$", linewidth=1)
# ax.plot(a_values, peri_curve(40), color ="k", label=r"$q=40\,\text{au}$", linewidth=1)
# ax.plot(a_values, peri_curve(45), color ="k", label=r"$q=45\,\text{au}$", linewidth=1)

ax.set_xlabel(r"semi-major axis $a$ [au]")
ax.set_ylabel(r"eccentricity $e$")

ax.set_ylim(0,1.05)

plt.savefig("plots/abs-Delta_e-10Myr.pdf")

plt.show()
