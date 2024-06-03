# If you divide the e(q)-space into n squares how many objects are in each square for every timestep?
#%% libraries
import numpy as np
import matplotlib.pyplot as plt
import rebound
from load_orb_elem import LoadOrbitalElements
plt.style.use("custom.mplstyle")

#%% get data from archive files
orbital_elements = LoadOrbitalElements("archives/mercurius-dtmin1-1e7-low.bin", "archives/mercurius-dtmin1-1e7-high.bin")

sa1 = orbital_elements["sa1"]
sa2 = orbital_elements["sa2"]

e_initial = orbital_elements["e_initial"]
a_initial = orbital_elements["a_initial"]
q_initial = orbital_elements["q_initial"]

e_final = orbital_elements["e_final"]
a_final = orbital_elements["a_final"]
q_final = orbital_elements["q_final"]


### Version with fixed, hardcoded edges
xedges = np.linspace(20, 45, 4)
yedges = np.linspace(0, 1, 4)


print("xedges: ",xedges)
print("yedges: ",yedges)

H_t = []

for s in range(len(sa1)):
    if s % 10 == 0:
        sim1 = sa1[s]
        sim2 = sa2[s]  # iterate through each snapshot in sa
        ps1 = sim1.particles  # simplify  referencing
        ps2 = sim2.particles
        # create list of orbital elements of all objects using list comprehension
        x_data = [ps1[i].orbit(primary=sim1.particles[0]).a * (1 - ps1[i].orbit(primary=sim1.particles[0]).e) for i in range(1, len(ps1))] + [ps2[i].orbit(primary=sim2.particles[0]).a * (1 - ps2[i].orbit(primary=sim2.particles[0]).e) for i in range(1, len(ps2))]  # perihelion distance
        y_data = [ps1[i].orbit(primary=sim1.particles[0]).e for i in range(1, len(ps1))] + [ps2[i].orbit(primary=sim2.particles[0]).e for i in range(1, len(ps2))]  # eccentricity
        H, xedges, yedges = np.histogram2d(x_data, y_data, bins=(xedges, yedges))
        # correct to original coordinates
        H = H.transpose()
        H_t.append(H)

# example colormesh of one 2d histogram
fig, ax = plt.subplots(constrained_layout=True)
X, Y = np.meshgrid(xedges, yedges)
# plot bars for 2d histogram
mesh = ax.pcolormesh(X, Y, H_t[0])
cbar = fig.colorbar(mesh, ax=ax)
cbar.set_label('number of particles')#, rotation=270)
# scatterplots of orbital elements of actual objects
ax.scatter(q_initial, e_initial,s=3)
ax.scatter(q_final, e_final,s=3)
ax.set_ylim(0,1)
ax.set_xlim(20,45)
ax.set_xlabel(r"$q$ [au]")
ax.set_ylabel(r"$e$")


# look at time evolution of number of particles in one square
def num_evo(i,j):
    # i,j: indices of matrix of 2d histogram
    n_t = []
    for k in range(1,len(H_t)):
        # use difference between initial matrix
        n_t.append(H_t[k][i,j]-H_t[1][i,j])
    return n_t

# plot evolution of number of particles in each square
t_steps = np.arange(1, len(H_t), 1)
fig, axs = plt.subplots(3,3, constrained_layout=True, sharex=True, sharey=False)
#plt.title("Evolution of number of particles for 1Myr")
for i in range(3):
    for j in range(3):
        axs[i,j].plot(t_steps, num_evo(i, j))

        axs[i, j].set_title(f"square ({i}, {j})")
        if i > 1:
            axs[i, j].set_xlabel("number of timesteps")
        if j == 0:
            axs[i, j].set_ylabel(r"$\Delta n$")
        axs[i,j].set_xscale('log')
#plt.savefig("plots/nep-Dn-1Myr-fewer-snapshots.pdf")
plt.show()
