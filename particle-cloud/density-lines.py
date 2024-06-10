# %% libraries
import numpy as np
import matplotlib.pyplot as plt
from load_orb_elem import LoadOrbitalElements
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks, peak_widths

plt.style.use("custom.mplstyle")

# %% get data from archive files
# orbital_elements = LoadOrbitalElements("archives/testa.bin", "archives/testb.bin")
# orbital_elements = LoadOrbitalElements("archives/jup-mercurius-dtmin05-1e6-low.bin", "archives/jup-mercurius-dtmin05-1e6-high.bin")
orbital_elements = LoadOrbitalElements(
    "archives/mercurius-testa-10e6.bin", "archives/mercurius-testb-10e6.bin"
)

sa1 = orbital_elements["sa1"]
sa2 = orbital_elements["sa2"]


## Settings for used KDE (for higher values of "fac" oscillations increase)
def my_kde_bandwidth(obj, fac=1):
    """We use Scott's Rule, multiplied by a constant factor."""
    return np.power(obj.n, -1.0 / (obj.d + 4)) * fac


def LoopOverSim(sa1, sa2):
    # empty arrays for storage of time-evolutions
    # initialise output parameters to avoid unbound error
    kde_values = []
    y_values = []
    y_t = []
    q_min = 0
    q_max = 0
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
            q_max = 30
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

    return kde_values, y_values, y_t, q_min, q_max


kde_values, y_values, y_t, q_min, q_max = LoopOverSim(sa1, sa2)

# Finding peaks and FWHM
fwhm_t = []
peak_t = []
peak_x = []
peak_y = []
width_heights = []

for i in range(len(kde_values)):
    if i % 1 == 0:
        # Find peaks
        peaks, _ = find_peaks(kde_values[i] / kde_values[0])
        if len(peaks) > 0:
            # Take the highest peak
            peak = peaks[np.argmax(kde_values[i][peaks] / kde_values[0][peaks])]
            peak_x = y_values[i][peak]
            peak_y = kde_values[i][peak] / kde_values[0][peak]

            # Find FWHM
            widths, width_heights, left_ips, right_ips = peak_widths(
                kde_values[i] / kde_values[0], [peak], rel_height=0.5
            )
            left_ip = y_values[i][int(left_ips[0])]
            right_ip = y_values[i][int(right_ips[0])]

            # save peaks and fwhm
            left_ip = np.array(left_ip)
            right_ip = np.array(right_ip)

            peak_t.append(peak_x)
            fwhm_t.append(abs(left_ip - right_ip))

# Plot time evolution of FWMH and peaks

fig, axs = plt.subplots(2, 1, sharex=True)
axs[0].scatter(np.arange(0, len(fwhm_t), 1), fwhm_t, s=2)
axs[0].set_ylabel(r"$\Delta e_\text{FWHM}$")
axs[1].scatter(np.arange(0, len(fwhm_t), 1), peak_t, s=2)
axs[1].set_ylabel(r"$e_{\text{peak}}$")
# axs[1].hlines(0.25, 0, len(fwhm_t), colors="tab:orange", linestyles="dashed")
axs[1].set_xlabel(r"timesteps")
plt.savefig(f"plots/peaks-fwhm-evolution-q{q_min}_{q_max}.pdf")

# Plotting setup
fig, ax = plt.subplots()
colors = plt.cm.hot(np.linspace(0, 1, len(kde_values)))

for i in range(len(kde_values)):
    if i % 100 == 0:
        # plot KDE densities normalised with distribiton at t=0
        ax.plot(
            y_values[i],
            kde_values[i] / kde_values[0],
            color=colors[i],
            label=f"timestep: {i*10}",
        )
        ax.set_ylabel(r"$n_{\text{KDE}}/n^0_\text{KDE}$")
        ax.set_xlabel(r"e")
        # ax.plot(y_t[i], np.zeros(y_t[i].shape), "b+", ms=10)

        # Find peaks
        peaks, _ = find_peaks(kde_values[i] / kde_values[0])
        if len(peaks) > 0:
            # Take the highest peak
            peak = peaks[np.argmax(kde_values[i][peaks] / kde_values[0][peaks])]
            peak_x = y_values[i][peak]
            peak_y = kde_values[i][peak] / kde_values[0][peak]
            ax.plot(peak_x, peak_y, "ro")

            # Find FWHM
            widths, width_heights, left_ips, right_ips = peak_widths(
                kde_values[i] / kde_values[0], [peak], rel_height=0.5
            )
            left_ip = y_values[i][int(left_ips[0])]
            right_ip = y_values[i][int(right_ips[0])]
            ax.plot([left_ip, right_ip], [width_heights[0]] * 2, "g--")
            ax.text(
                (left_ip + right_ip) / 2,
                width_heights[0] + 0.005,
                f"FWHM: {right_ip - left_ip:.2f}",
                color="g",
            )

ax.legend()
plt.savefig(f"plots/example-peaks-fwhm-q{q_min}_{q_max}.pdf")


# calculate peaks and FWHM
peaks, _ = find_peaks(kde_values[10])
results_half = peak_widths(kde_values[10], peaks, rel_height=0.5)
results_full = peak_widths(kde_values[10], peaks, rel_height=1)

"""
fig, ax = plt.subplots()
for i in range(len(kde_values)):
    if i == 10:
        ax.plot(
            (y_values[i]),
            kde_values[i] / kde_values[0],
            label=f"timestep: {i}",
        )
        # ax.plot(y_t[i], np.zeros(y_t[i].shape), "b+", ms=10)
        plt.hlines(*results_full[1:], color="C3")
        plt.plot(peaks, kde_values[peaks], "x")
# ax.set_xscale("log")
# ax.set_yscale("log")
ax.set_ylim(-0.05, 2)
ax.legend()
#plt.savefig(f"plots/example-peaks-fwhm-q{q_min}_{q_max}.pdf")
"""

## Animation for every timestep

from matplotlib.animation import FuncAnimation

# Animation setup
fig, ax = plt.subplots()
(line,) = ax.plot([], [], label="timestep: 0")
(scatter,) = ax.plot([], [], "b+", ms=10)
ax.set_xlim(0, 1)
ax.set_ylim(-0.05, max(max(kde) for kde in kde_values / kde_values[0]))
ax.set_xlabel(r"$e$")
ax.set_ylabel(r"$n_{\text{KDE}}/n^0_\text{KDE}$")


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
ani.save(
    f"gifs/density-evolution-q-{q_min}_{q_max}.mp4",
    writer="ffmpeg",
    fps=24,
    bitrate=1500,
)
print("finished saving animation")
plt.close()
