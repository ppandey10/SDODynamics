#%% Libraries
import rebound
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

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

#%% setup simulation data

# Function for Hill-Radius
def r_hill(a, M_small_body, M_central_body):
    return a * (M_small_body / 3 * (M_central_body)) ** (1/3)

# Masses and Semi-major axis of Sun and Neptune
M_sun = 1
M_neptune = 5.151383772628957e-05
a_neptune = 30.237878368382898

sa = rebound.Simulationarchive("archives/archive-t_1e6-mercurius-10_rhill_01.bin")

#%% Plot for e(q) change
# Initialize a figure and axis
fig, ax = plt.subplots(figsize=(5, 5), dpi=300, constrained_layout=True)

# Initialize the plot objects
line, = ax.plot([], [],'.', ms=3, color="tab:blue")
ax.vlines(a_neptune + 10 * r_hill(a_neptune, M_neptune, M_sun), 0, 1.25, colors='r', linestyles='dashed')
ax.vlines(a_neptune - 10 * r_hill(a_neptune, M_neptune, M_sun), 0, 1.25, colors='r', linestyles='dashed')
ax.set_xlabel("pericentre distance $q$ [au]")
ax.set_ylabel("eccentricity $e$")

# Define the initialization function
def init():
    line.set_data([], [])  # Clear the previous plot
    return line,

# Define the update function
def update(s):
    sim = sa[s+10]  # iterate through each snapshot in sa
    ps = sim.particles  # intermediate object to simplify the referencing
    x_data = [ps[i].orbit(primary=sim.particles[0]).a * (1-ps[i].orbit(primary=sim.particles[0]).e) for i in range(1, len(ps))]  # perihelion distance
    y_data = [ps[i].orbit(primary=sim.particles[0]).e for i in range(1, len(ps))]  # eccentricity
    line.set_data(x_data, y_data)
    ax.set_title(r"$t_{\text{sim}}=$ + "f"{int(sa[s+10].t)} yr")  # Add a title with the frame number
    # Set axis limits
    ax.set_xlim(5, 60)
    ax.set_ylim(0, 1.25)
    return line,

# Create the animation
ani = FuncAnimation(fig, update, frames=len(sa), init_func=init, blit=True)

# Save the animation as a GIF
ani.save('gifs/q-t_1e6-mercurius-10_rhill_01_term.gif', writer='pillow', fps=24)

plt.close()

print('First gif saved')

'''fig, ax = plt.subplots(figsize=(5, 5), dpi=350, constrained_layout=True)

# Initialize the plot objects
line, = ax.plot([], [],'.', ms=3, color="tab:blue")
ax.set_xlabel("semi major axis $a$ [au]")
ax.set_ylabel("eccentricity $e$")

# Define the initialization function
def init():
    line.set_data([], [])  # Clear the previous plot
    return line,

# Define the update function
def update(s):
    sim = sa[s]  # iterate through each snapshot in sa
    ps = sim.particles  # intermediate object to simplify the referencing
    x_data = [ps[i].orbit(primary=sim.particles[0]).a for i in range(1, len(ps))]  # perihelion distance
    y_data = [ps[i].orbit(primary=sim.particles[0]).e for i in range(1, len(ps))]  # eccentricity
    line.set_data(x_data, y_data)
    ax.set_title(f"time: {int(sa[s].t)} yr")  # Add a title with the frame number
    # Set axis limits
    ax.set_xlim(20, max(x_data))
    ax.set_ylim(0, 1)
    return line,

# Create the animation
ani = FuncAnimation(fig, update, frames=len(sa), init_func=init, blit=True)

# Save the animation as a GIF
ani.save('gifs/a-t_1e6-mercurius-10_rhill_01.gif', writer='pillow', fps=24)

plt.close()

print('second gif saved')'''