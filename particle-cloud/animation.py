#%% Libraries
import rebound
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

plt.style.use("custom.mplstyle")

#%% setup simulation data

# Function for Hill-Radius
def r_hill(a, M_small_body, M_central_body):
    return a * (M_small_body / 3 * (M_central_body)) ** (1/3)

# Masses and Semi-major axis of Sun and Neptune
M_sun = 1
M_neptune = 5.151383772628957e-05
a_neptune = 30.237878368382898

sa = rebound.Simulationarchive("archives/mercurius_tsim_1e7_dtias15_5e-4_dt_5_rhill_2_final.bin")

#%% Plot for e(q) change
# Initialize a figure and axis
fig, ax = plt.subplots()

# Initialize the plot objects
line, = ax.plot([], [],'.', ms=3, color="tab:blue")
# ax.vlines(a_neptune + 11 * r_hill(a_neptune, M_neptune, M_sun), 0, 1.25, colors='r', linestyles='dashed')
# ax.vlines(a_neptune - 5 * r_hill(a_neptune, M_neptune, M_sun), 0, 1.25, colors='r', linestyles='dashed')
ax.set_xlabel("pericentre distance $q$ [au]")
ax.set_ylabel("eccentricity $e$")

# Define the initialization function
def init():
    line.set_data([], [])  # Clear the previous plot
    return line,

# Define the update function
def update(s):
    if s % 1 == 0:
        sim = sa[s]  # iterate through each snapshot in sa
        ps = sim.particles  # intermediate object to simplify the referencing
        x_data = [ps[i].orbit(primary=sim.particles[0]).a * (1-ps[i].orbit(primary=sim.particles[0]).e) for i in range(1, len(ps))]  # perihelion distance
        y_data = [ps[i].orbit(primary=sim.particles[0]).e for i in range(1, len(ps))]  # eccentricity
        line.set_data(x_data, y_data)
        ax.set_title(r"$t_{\text{sim}}=$ + "f"{int(sa[s].t)} yr")  # Add a title with the frame number
        # Set axis limits
        ax.set_xlim(5, 60)
        ax.set_ylim(0, 1.25)
    return line,

# Create the animation
ani = FuncAnimation(fig, update, frames=len(sa), init_func=init, blit=True)

ani.save(
   "gifs/1e7evolution-q.mp4",
   writer="ffmpeg",
   fps=24,
   bitrate=1500,
)# Save the animation as a GIF

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
