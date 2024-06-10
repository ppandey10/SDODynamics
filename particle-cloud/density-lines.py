# If you divide the e(q)-space into n squares how many objects are in each square for every timestep?
# %% libraries
import numpy as np
import matplotlib.pyplot as plt
from load_orb_elem import LoadOrbitalElements
from scipy.stats import gaussian_kde

plt.style.use("custom.mplstyle")

# %% get data from archive files
# orbital_elements = LoadOrbitalElements("archives/testa.bin", "archives/testb.bin")
# orbital_elements = LoadOrbitalElements("archives/jup-mercurius-dtmin05-1e6-low.bin", "archives/jup-mercurius-dtmin05-1e6-high.bin")
orbital_elements = LoadOrbitalElements(
    "archives/mercurius-testa-10e6.bin", "archives/mercurius-testb-10e6.bin"
)

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


# %% Version with fixed, hardcoded edges
nbins = 4
qmin = 10
qmax = 50
emin = 0
emax = 1
xedges = np.linspace(qmin, qmax, nbins + 1)
yedges = np.linspace(emin, emax, nbins + 1)


print("xedges: ", xedges)
print("yedges: ", yedges)


## Settings for used KDE
def my_kde_bandwidth(obj, fac=1):
    """We use Scott's Rule, multiplied by a constant factor."""
    return np.power(obj.n, -1.0 / (obj.d + 4)) * fac


def LoopOverSim(sa1, sa2, xedges, yedges):
    # empty arrays for storage of time-evolutions
    result = "Dictionary not defined"
    # initialise output parameters to avoid unbound error
    kde_values = []
    y_values = []
    y_t = []
    # loop
    for s in range(len(sa1)):
        if s % 2 == 0:
            sim1 = sa1[s]
            sim2 = sa2[s]  # iterate through each snapshot in sa
            ps1 = sim1.particles  # simplify  referencing
            ps2 = sim2.particles
            # create list of orbital elements of all objects using list comprehension
            q_data = [
                ps1[i].orbit(primary=sim1.particles[0]).a
                * (1 - ps1[i].orbit(primary=sim1.particles[0]).e)
                for i in range(1, len(ps1))
            ] + [
                ps2[i].orbit(primary=sim2.particles[0]).a
                * (1 - ps2[i].orbit(primary=sim2.particles[0]).e)
                for i in range(1, len(ps2))
            ]  # perihelion distance
            e_data = [
                ps1[i].orbit(primary=sim1.particles[0]).e for i in range(1, len(ps1))
            ] + [
                ps2[i].orbit(primary=sim2.particles[0]).e for i in range(1, len(ps2))
            ]  # eccentricity

            # calculate density in defined bin
            q_min = 29
            q_max = 31
            q_data = np.array(q_data)
            e_data = np.array(e_data)

            # Filter points within the interval
            mask = (q_data >= q_min) & (q_data <= q_max)
            x_filtered = q_data[mask]
            y_filtered = e_data[mask]

            if len(y_filtered) > 0:
                kde = gaussian_kde(y_filtered, bw_method=my_kde_bandwidth)
                y = np.linspace(min(y_filtered), max(y_filtered), 1000)
                kde_values.append(kde(y))
                y_values.append(y)
                y_t.append(y_filtered)

    return kde_values, y_values, y_t


kde_values, y_values, y_t = LoopOverSim(sa1, sa2, xedges, yedges)


fig, ax = plt.subplots()
for i in range(len(kde_values)):
    if i % 1000 == 0:
        ax.plot(
            (1 - y_values[i]) ** 2,
            kde_values[i] / kde_values[0],
            label=f"timestep: {i}",
        )
        ax.plot(y_t[i], np.zeros(y_t[i].shape), "b+", ms=10)
# ax.set_xscale("log")
# ax.set_yscale("log")
# ax.set_ylim(-0.05, 2)
ax.legend()

# Animation for every timestep

from matplotlib.animation import FuncAnimation

# Animation setup
fig, ax = plt.subplots()
(line,) = ax.plot([], [], label="timestep: 0")
(scatter,) = ax.plot([], [], "b+", ms=10)
ax.set_xlim(emin, emax)
ax.set_ylim(-0.05, max(max(kde) for kde in kde_values))
ax.set_xlabel(r"$e$")
ax.set_ylabel(r"$n_{\text{KDE}}$")


def init():
    line.set_data([], [])
    scatter.set_data([], [])
    return line, scatter


def update(frame):
    line.set_data(y_values[frame], kde_values[frame] / kde_values[0])
    scatter.set_data(y_t[frame], np.zeros(y_t[frame].shape))
    line.set_label(f"timestep: {frame*10}")
    ax.legend()
    return line, scatter


ani = FuncAnimation(fig, update, frames=len(kde_values), init_func=init, blit=True)
print("start saving animation")
ani.save("gifs/density-evolutionq29-31.mp4", writer="ffmpeg", fps=24, bitrate=1500)
print("finished saving animation")
plt.show()
