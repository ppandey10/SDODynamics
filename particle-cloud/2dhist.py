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

#%% import data
sa1 = rebound.Simulationarchive("archives/5e6-dt50-rhill2-elow.bin")
sa2 = rebound.Simulationarchive("archives/5e6-dt100-rhill2-ehigh.bin")

max_1 = len(sa1)
max_2 = len(sa2)

print(sa1[-1])
print(sa2[-1])

# load initial results
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

for k in range(2, i_last):
    e_initial2.append(sim2_initial.particles[k].e)
    a_initial2.append(sim2_initial.particles[k].a)

e_initial1 = np.array(e_initial1)
a_initial1 = np.array(a_initial1)
q_initial1 = a_initial1 * (1 - e_initial1)

e_initial2 = np.array(e_initial2)
a_initial2 = np.array(a_initial2)
q_initial2 = a_initial2 * (1 - e_initial2)

# add arrays for both simulations
e_initial = np.append(e_initial1, e_initial2)
a_initial = np.append(a_initial1, a_initial2)
q_initial = a_initial * (1 - e_initial)

# load final results
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
q_final1 = a_final1 * (1 - e_final1)

e_final2 = np.array(e_final2)
a_final2 = np.array(a_final2)
q_final2 = a_final2 * (1 - e_final2)

# add arrays for both simulations
e_final = np.append(e_final1, e_final2)
a_final = np.append(a_final1, a_final2)
q_final = a_final * (1 - e_final)

#%% histogram

e = e_final
q = q_final

xedges = np.arange(int(min(q)) - 0.25, int(max(q)) + 1.75, 1)
yedges = np.arange(-0.01, 0.94, 0.05)
H, xedges, yedges = np.histogram2d(q, e, bins=(xedges, yedges))
H = H.transpose()

fig, ax = plt.subplots(constrained_layout=True)
X, Y = np.meshgrid(xedges, yedges)
ax.pcolormesh(X, Y, H)
ax.scatter(q, e, color='black')
plt.show()

#%% substract both matrices
H1, xedges1, yedges1 = np.histogram2d(q_initial, e_initial, bins=(xedges, yedges))
H2, xedges2, yedges2 = np.histogram2d(q_final, e_final, bins=(xedges, yedges))

H_diff = H1 - H2
H_diff.transpose()

fig, ax = plt.subplots(constrained_layout=True)
X, Y = np.meshgrid(xedges, yedges)
ax.pcolormesh(X, Y, H)
plt.show()

#%% Iterate through timesteps and save matrix H
sa1 = rebound.Simulationarchive("archives/nep-mercurius-dtmin1-1e6-high.bin")
sa2 = rebound.Simulationarchive("archives/nep-mercurius-dtmin1-1e6-low.bin")

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
y_min = min(e_initial)
y_max = max(e_initial)

yedges = np.linspace(y_min, y_max, 4)  # 4 edges to create 3 bins
print("xedges: ",xedges)
print("yedges: ",yedges)

H_t = []

for s in range(len(sa1)):
    if s % 50 == 0:
        sim1 = sa1[s]
        sim2 = sa2[s]  # iterate through each snapshot in sa
        ps1 = sim1.particles  # intermediate object to simplify the referencing
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
        #axs[i,j].set_xscale('log')
plt.savefig("plots/nep-Dn-1Myr-fewer-snapshots.pdf")
plt.show()



#%% Histogram Animation

from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots(figsize=(5, 5), dpi=100, constrained_layout=True)

# Initialize the plot object for the 2D histogram (empty at start)
xedges = np.linspace(5, 60, 101)  # 100 bins in the range [5, 60]
yedges = np.linspace(0, 1.25, 101)  # 100 bins in the range [0, 1.25]
X, Y = np.meshgrid(xedges, yedges)
mesh = ax.pcolormesh(X, Y, np.zeros((100, 100)), cmap='viridis', shading='auto')

ax.set_xlabel("pericentre distance $q$ [au]")
ax.set_ylabel("eccentricity $e$")

# Define the initialization function
def init():
    mesh.set_array(np.zeros((100, 100)).ravel())  # Clear the previous plot
    return mesh,

# Define the update function
def update(s):
    if s % 10000 == 0:
        sim1 = sa1[s]
        sim2 = sa2[s]  # iterate through each snapshot in sa
        ps1 = sim1.particles  # intermediate object to simplify the referencing
        ps2 = sim2.particles
        x_data = [ps1[i].orbit(primary=sim1.particles[0]).a * (1 - ps1[i].orbit(primary=sim1.particles[0]).e) for i in range(1, len(ps1))] + [ps2[i].orbit(primary=sim2.particles[0]).a * (1 - ps2[i].orbit(primary=sim2.particles[0]).e) for i in range(1, len(ps2))]  # perihelion distance
        y_data = [ps1[i].orbit(primary=sim1.particles[0]).e for i in range(1, len(ps1))] + [ps2[i].orbit(primary=sim2.particles[0]).e for i in range(1, len(ps2))]  # eccentricity

        # Create a 2D histogram
        hist, xedges, yedges = np.histogram2d(x_data, y_data, bins=[100, 100], range=[[5, 60], [0, 1.25]])

        # Update the mesh with new data
        mesh.set_array(hist.T.ravel())
        mesh.set_clim([hist.min(), hist.max()])  # Adjust color limits

        # Update the title and axes
        ax.set_title(r"$t_{\text{sim}}=$" + f"{int(sa1[s].t)} yr")

    return mesh,

# Create the animation
ani = FuncAnimation(fig, update, frames=len(sa1), init_func=init, blit=True)

# Save the animation as a GIF
ani.save('gifs/2dhist.gif', writer='pillow', fps=24)

plt.close()

#%% Animation
from matplotlib.animation import FuncAnimation

# Plot for e(q) change
# Initialize a figure and axis
fig, ax = plt.subplots(figsize=(5, 5), dpi=100, constrained_layout=True)

# Initialize the plot objects
line, = ax.plot([], [], '.', ms=3, color="tab:blue")
#ax.vlines(a_neptune + 10 * r_hill(a_neptune, M_neptune, M_sun), 0, 1.25, colors='r', linestyles='dashed')
#ax.vlines(a_neptune - 10 * r_hill(a_neptune, M_neptune, M_sun), 0, 1.25, colors='r', linestyles='dashed')
ax.set_xlabel("pericentre distance $q$ [au]")
ax.set_ylabel("eccentricity $e$")


# Define the initialization function
def init():
    line.set_data([], [])  # Clear the previous plot
    return line,


# Define the update function
def update(s):
    for s in range(len(sa.t)):
        if s % 10000 == 0:
            sim1 = sa1[s]
            sim2 = sa2[s]  # iterate through each snapshot in sa
            ps1 = sim1.particles  # intermediate object to simplify the referencing
            ps2 = sim2.particles
            x_data = [ps1[i].orbit(primary=sim1.particles[0]).a * (1 - ps1[i].orbit(primary=sim1.particles[0]).e) for i in range(1, len(ps1))] + [ps2[i].orbit(primary=sim2.particles[0]).a * (1 - ps2[i].orbit(primary=sim2.particles[0]).e) for i in range(1, len(ps2))]  # perihelion distancey
            y_data = [ps1[i].orbit(primary=sim1.particles[0]).e for i in range(1, len(ps1))] + [ps2[i].orbit(primary=sim2.particles[0]).e for i in range(1, len(ps2))]  # eccentricity
            line.set_data(x_data, y_data)
            ax.set_title(r"$t_{\text{sim}}=$ + "f"{int(sa1[s].t)} yr")  # Add a title with the frame number
            # Set axis limits
            ax.set_xlim(5, 60)
            ax.set_ylim(0, 1.25)
            return line,


# Create the animation
ani = FuncAnimation(fig, update, frames=len(sa1), init_func=init, blit=True)

# Save the animation as a GIF
ani.save('gifs/2dhist.gif', writer='pillow', fps=24)

plt.close()
