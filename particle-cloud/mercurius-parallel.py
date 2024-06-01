from numba import jit
import rebound
import numpy as np
from threading import Thread
from tqdm import tqdm

# %% Simulation
def init_simulation(scatterer):
    sim = rebound.Simulation()
    sim.units = ('yr', 'AU', 'Msun')
    sim.add('Sun')
    sim.add(scatterer)
    return sim

sim1 = init_simulation("Jupiter")
sim2 = init_simulation("Jupiter")

print("added Sun and scatterer")

# Functions for Hill-Radius and distribution
def r_hill(a, M_body, M_center):
    return a * (M_body / (3 * M_center)) ** (1 / 3)

@jit(nopython=True)
def generate_uniform_distribution(e, q):
    results = []
    for ei in e:
        for qi in q:
            results.append([qi, ei, qi / (1 - ei)])
    return np.array(results)

# Integration function
def custom_integration(e, a, dt, r_hill_crit, t_sim, filenames):
    integrator = "mercurius"
    print(f"Using {integrator} integrator")

    sim1.integrator = sim2.integrator = integrator
    progress_bar = [tqdm(total=t_sim, desc=f"Thread {i + 1} Progress", position=i) for i in range(2)]

    def integrate(sim, e_filter, progress, filename):
        for ei, ai in zip(e, a):
            if e_filter(ei):
                rand = np.random.random() * 2 * np.pi
                sim.add(a=ai, e=ei, primary=sim.particles[0], Omega=0, omega=rand, f=rand)
        print(f"N particles ({filename}): {sim.N}")

        sim.dt = dt if e_filter(0.4) else 1.5 * dt
        sim.ri_ias15.min_dt = 0.1
        sim.ri_mercurius.r_crit_hill = r_hill_crit
        sim.move_to_com()
        E0 = sim.energy()

        removed_particles = 0
        for ts in range(0, int(t_sim), 1):
            sim.integrate(ts, exact_finish_time=0)
            if ts % 1000 == 0:  # Adjust the frequency of removals
                sim.save_to_file(filename)
                # Remove particles with e > 1.02
                # TODO: additional condition for distance of particle after close encounter
                particles_to_remove = [p for p in sim.particles[1:] if p.e > 1.02]
                for p in particles_to_remove:
                    sim.remove(p)
                    removed_particles += 1
            progress.update(1)

        dE = abs((sim.energy() - E0))
        print(f"dE ({filename}) {dE}")
        print(f"Removed Particles ({filename}): {removed_particles}")

    threads = [Thread(target=integrate, args=(sim, f, p, fn)) for sim, f, p, fn in zip(
        [sim1, sim2], [lambda e: e <= 0.4, lambda e: e > 0.4], progress_bar, filenames)]

    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()

# Integrate
a_scatterer = sim1.particles[1].a
M_scatterer = sim1.particles[1].m
M_sun = sim1.particles[0].m

e_generate = np.arange(0, 0.95, 0.05)
q_generate = np.arange(a_scatterer - (5 * r_hill(a_scatterer, M_scatterer, M_sun)),
                       a_scatterer + (11 * r_hill(a_scatterer, M_scatterer, M_sun)), 0.5)

orbital_elements = generate_uniform_distribution(e_generate, q_generate)

print(f"added {len(orbital_elements)} test particles")

custom_integration(orbital_elements[:, 1], orbital_elements[:, 2], 0.5, 3, 1e6,
               ["archives/testa.bin", "archives/testb.bin"])
