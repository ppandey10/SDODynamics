# If you divide the e(q)-space into n squares how many objects are in each square for every timestep?
#%% libraries
import numpy as np
import matplotlib.pyplot as plt
import rebound
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
plt.rc('axes', labelsize=11)  #fontsize of the x and y labels
plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}\usepackage{siunitx}\usepackage{newtx}"
plt.rcParams['errorbar.capsize'] = 4  # Einstellung Dicke der Fehlerbalkenenden
plt.rc('axes', labelsize=11, titlesize=11)
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


#%% Iterate through timesteps and save matrix H
sa1 = rebound.Simulationarchive("archives/mercurius-dtmin1-1e7-low.bin")
sa2 = rebound.Simulationarchive("archives/mercurius-dtmin1-1e7-high.bin")

# boundaries for edges from initial distribution
sim1_initial = sa1[-len(sa1)]
sim2_initial = sa2[-len(sa2)]

# last particle index
i_last1 = sa1[1].N
i_last2 = sa2[1].N

# create arrays for initial a, e and q
e_initial1 = e_initial2 = a_initial1 = a_initial2 = np.empty([])

for i in range(2, i_last1):
    e_initial1 = np.append(e_initial1, sim1_initial.particles[i].e)
    a_initial1 = np.append(a_initial1, sim1_initial.particles[i].a)

for k in range(2, i_last2):
    e_initial2 = np.append(e_initial2, sim2_initial.particles[k].e)
    a_initial2 = np.append(a_initial2, sim2_initial.particles[k].a)

q_initial1 = a_initial1 * (1 - e_initial1)
q_initial2 = a_initial2 * (1 - e_initial2)

q_initial = np.append(q_initial1, q_initial2)
a_initial = np.append(a_initial1, a_initial2)
e_initial = np.append(e_initial1, e_initial2)

# set edges of histogram
#q_initial = q_initial[1:]
x_min = 26.328623 #min(q_initial)
x_max = max(q_initial)


xedges = np.linspace(x_min, x_max, 4)  # 4 edges to create 3 bins

# Calculate y-axis edges
y_min = 0.000207#min(e_initial)
y_max = max(e_initial)

yedges = np.linspace(y_min, y_max, 4)  # 4 edges to create 3 bins
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
ax.pcolormesh(X, Y, H_t[0])
ax.set_xlabel(r"$q$ [au]")
ax.set_ylabel(r"$e$")
plt.show()

# look at time evolution of number of particles in one square
def num_evo(i,j):
    # i,j: indices of matrix of 2d histogram
    n_t = []
    for k in range(1,len(H_t)):
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
