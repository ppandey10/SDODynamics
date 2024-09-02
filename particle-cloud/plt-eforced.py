import numpy as np
import rebound
import matplotlib.pyplot as plt


plt.style.use("custom.mplstyle")
# load simulation archive
sa = rebound.Simulationarchive("archives/mercurius-e0-amin-10e6.bin")

# loop over archive


def LoopOverSim(sa):
    # empty arrays for storage of time-evolutions
    # initialise output parameters to avoid unbound error
    pomega_t = []
    e_t = []
    a_t = []
    q_t = []
    # loop
    for s in range(len(sa)):
        if s % 2 == 0:
            sim = sa[s]
            ps = sim.particles  # simplify  referencing
            # create list of orbital elements of all objects using list comprehension
            q_data = [
                ps[i].orbit(primary=sim.particles[0]).a
                * (1 - ps[i].orbit(primary=sim.particles[0]).e)
                for i in range(1, len(ps))
            ]  # perihelion distance
            e_data = [
                ps[i].orbit(primary=sim.particles[0]).e for i in range(1, len(ps))
            ]  # eccentricity

            pomega_data = [
                ps[i].orbit(primary=sim.particles[0]).pomega for i in range(1, len(ps))
            ]  # longitude of pericentre

            a_data = [
                ps[i].orbit(primary=sim.particles[0]).a for i in range(1, len(ps))
            ]  # semi-major axis

            # save data
            e_t.append(e_data)
            a_t.append(a_data)
            q_t.append(q_data)
            pomega_t.append(pomega_data)

    # transform to numpy array (faster)
    e_t = np.array(e_t)
    a_t = np.array(a_t)
    q_t = np.array(q_t)
    pomega_t = np.array(pomega_t)

    result_data = {
        "pomega_t": pomega_t,
        "e_t": e_t,
        "a_t": a_t,
        "q_t": q_t,
    }

    return result_data


print(sa)

sim_data = LoopOverSim(sa)

pomega_t = sim_data["pomega_t"]
e_t = sim_data["e_t"]

print(e_t[1])

fig, ax = plt.subplots()
ax.scatter(e_t[1] * np.cos(pomega_t[1]), e_t[1] * np.sin(pomega_t[1]), s=2)
plt.show()
