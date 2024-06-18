# If you divide the e(q)-space into n squares how many objects are in each square for every timestep?
#%% libraries
import numpy as np
from numba import jit
import matplotlib.pyplot as plt
from load_orb_elem import LoadOrbitalElements
plt.style.use("custom.mplstyle")

#%% get data from archive files
#orbital_elements = LoadOrbitalElements("archives/testa.bin", "archives/testb.bin")
#orbital_elements = LoadOrbitalElements("archives/jup-mercurius-dtmin05-1e6-low.bin", "archives/jup-mercurius-dtmin05-1e6-high.bin")
orbital_elements = LoadOrbitalElements("archives/mercurius-testa-10e6.bin", "archives/mercurius-testb-10e6.bin")

sa1 = orbital_elements["sa1"]
sa2 = orbital_elements["sa2"]

print(sa1)
print(sa2)

e_initial = orbital_elements["e_initial"]
a_initial = orbital_elements["a_initial"]
q_initial = orbital_elements["q_initial"]

e_final = orbital_elements["e_final"]
a_final = orbital_elements["a_final"]
q_final = orbital_elements["q_final"]


#%% Version with fixed, hardcoded edges
nbins = 4
qmin = 10
qmax = 50
emin = 0
emax = 1
xedges = np.linspace(qmin, qmax, nbins+1)
yedges = np.linspace(emin, emax, nbins+1)


print("xedges: ",xedges)
print("yedges: ",yedges)


#@jit(nopython=True)
def LoopOverSim(sa1, sa2, xedges, yedges):
    # empty arrays for storage of time-evolutions
    H_t = []
    q_t_list = []
    e_t_list = []
    result = "Dictionary not defined" 
    
    # loop
    for s in range(len(sa1)):
        if s % 10 == 0:
            sim1 = sa1[s]
            sim2 = sa2[s]  # iterate through each snapshot in sa
            ps1 = sim1.particles  # simplify  referencing
            ps2 = sim2.particles
            # create list of orbital elements of all objects using list comprehension
            q_data = [ps1[i].orbit(primary=sim1.particles[0]).a * (1 - ps1[i].orbit(primary=sim1.particles[0]).e) for i in range(1, len(ps1))] + [ps2[i].orbit(primary=sim2.particles[0]).a * (1 - ps2[i].orbit(primary=sim2.particles[0]).e) for i in range(1, len(ps2))]  # perihelion distance
            e_data = [ps1[i].orbit(primary=sim1.particles[0]).e for i in range(1, len(ps1))] + [ps2[i].orbit(primary=sim2.particles[0]).e for i in range(1, len(ps2))]  # eccentricity
            H, xedges, yedges = np.histogram2d(q_data, e_data, bins=(xedges, yedges))
            # correct to original coordinates
            H = H.transpose()
            H_t.append(H)
            # save e- and q-value for corresponding time-steps
            q_t_list.append(q_data)
            e_t_list.append(e_data)
            q_t = np.array(q_t_list)
            e_t = np.array(e_t_list)

    result = {
        "H_t": H_t,
        "q_t": q_t,
        "e_t": e_t,
        }
    return result

data_all = LoopOverSim(sa1, sa2, xedges, yedges)
H_t = data_all["H_t"]
q_t = data_all["q_t"]
e_t = data_all["e_t"]


#print("q_t: ", q_t[0])
#print(H_t)

# plot population for one "line" of q distribution
q_t_single = []
e_t_single = []
equal_q_indices = []

# Identify indices where q_t[0][n] is approximately 30.3313
for n in range(len(q_t[1])):
    if round(q_t[0][n], 4) == 30.3313:
        equal_q_indices.append(n)

# Extract values for these indices
for k in range(len(q_t)):
    q_tmp = []
    e_tmp = []
    for ind in equal_q_indices:
        q_tmp.append(q_t[k][ind])
        e_tmp.append(e_t[k][ind])
    q_t_single.append(q_tmp) 
    e_t_single.append(e_tmp)
q_t_single = np.array(q_t_single)
e_t_single = np.array(e_t_single)

print(q_t_single[0])

# Plot lines
fig, ax = plt.subplots(constrained_layout=True)
for x in range(len(q_t_single)):
    if x  == 10 or x == 11:
    #if x  == 10:
        ax.scatter((1 - e_t_single[x][:]), q_t_single[x][:], s=10, label=f"ts {x}")
        #ax.plot((1 - e_t_single[x][:])**2, q_t_single[x][:] )
        ax.legend()
ax.set_xlabel(r"$(1-e)$")
ax.set_ylabel(r"$q$")

plt.show()

# print("qt: " , len(q_t[0]))
# print("q_t_single: ", q_t_single[0])
# print("length q_t_single: ", len(q_t_single[0]))


# 2d histogram at last timestep
fig, ax = plt.subplots(constrained_layout=True)
X, Y = np.meshgrid(xedges, yedges)
# plot bars for 2d histogram
mesh = ax.pcolormesh(X, Y, H_t[len(H_t)-1])
cbar = fig.colorbar(mesh, ax=ax)
cbar.set_label('number of particles')#, rotation=270)
# scatterplots of orbital elements of actual objects
ax.scatter(q_initial, e_initial,s=3)
ax.scatter(q_final, e_final,s=3)
ax.set_ylim(emin,emax)
ax.set_xlim(qmin,qmax)
ax.set_xlabel(r"$q$ [au]")
ax.set_ylabel(r"$e$")
#plt.savefig("plots/nep-merc-1Myr-4bins-hist.pdf")


# look at time evolution of number of particles in one square
def num_evo(i,j):
    # i,j: indices of matrix of 2d histogram
    n_t = []
    for k in range(1,len(H_t)):
        # use difference between initial matrix
        n_t.append(H_t[k][i,j])#-H_t[1][i,j])
    return n_t

# plot evolution of number of particles in each square
t_steps = np.arange(1, len(H_t), 1)
fig, axs = plt.subplots(nbins,nbins, constrained_layout=True, sharex=True, sharey=False, figsize=(16,8))
#plt.title("Evolution of number of particles for 1Myr")
for i in range(nbins):
    for j in range(nbins):
        axs[i,j].plot(t_steps, num_evo((nbins-1)-i, j)) # invert y-axis of plotting direction 
        #axs[i,j].scatter(t_steps, num_evo((nbins-1)-i, j))
        #axs[i, j].set_title(f"square ({(nbins-1)-i}, {j})")
        if i > 2:
            axs[i, j].set_xlabel("number of timesteps")
        if j == 0:
            axs[i, j].set_ylabel(r"$n$")
        #axs[i,j].set_xscale('log')
#plt.savefig("plots/nep-merc-1Myr-4bins-abs.pdf")
#plt.show()
plt.close()
