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

# initialize Simulations
sim1 = init_simulation("Neptune")
sim2 = init_simulation("Neptune")

print("added Sun and scatterer")

# Functions for Hill-Radius and distribution
def r_hill(a, M_body, M_center):
    return a * (M_body / (3 * M_center)) ** (1 / 3)


def generate_uniform_distribution(e, q):
    results = []
    for ei in e:
        for qi in q:
            results.append([qi, ei, qi / (1 - ei)])
    return np.array(results)

# Integration function
def custom_integration(e, a, dt, r_hill_crit, t_sim, filenames):
    
    """        
    Params:
    e:  (1d-array) eccentricities
    a:  (1d-array) with semi major axis
    dt: (double) timestep (keep units defined earlier in mind)
    r_hill_crit: (double) Hill-radius at with mercurius switches from whfast to ias15
    t_sim:  (int) total runtime of simulation
    filenames: (["file1.bin","file2.bin"]) filenames/paths for binary files

    Returns: 
    file1.bin
    file2.bin
    """
    
    # specify integrator
    integrator = "mercurius"
    #integrator = "whfast"
    print(f"Using {integrator} integrator")

    sim1.integrator = sim2.integrator = integrator

    # initialize progress bar
    progress_bar = [tqdm(total=t_sim, desc=f"Thread {i + 1} Progress", position=i) for i in range(2)]

    def integrate(sim, e_filter, progress, filename):
        # iterates through eccentricies and semi-major axis 
        # "zip" combines two arrays to one new (actually returns iterator)
        for ei, ai in zip(e, a):
            if e_filter(ei):
                rand = np.random.random() * 2 * np.pi
                sim.add(a=ai, e=ei, primary=sim.particles[0], Omega=0, omega=rand, f=rand)
        print(f"N particles ({filename}): {sim.N}")

        sim.dt = dt if e_filter(0.4) else 2 * dt
        print("dt: ", sim.dt)
        sim.ri_ias15.min_dt = 1e-6 # ensure that close encounters do not stall the integration 
        sim.ri_mercurius.r_crit_hill = r_hill_crit
        sim.move_to_com()
        E0 = sim.energy()

        # Integration and keeping remove ejected particles
        removed_particles = 0
        for ts in range(0, int(t_sim), 1): #int(sim.dt))
            sim.integrate(ts, exact_finish_time=0)
            if ts % 1000 == 0:  # Adjust the frequency of removals
                sim.save_to_file(filename)
                # Remove particles with e > 1.02 and distance of particle over 35 au two avoid artefacts
                # at very close encounters (eccentricity might be >1 for a short period of time)
                particles_to_remove = [p for p in sim.particles[1:] if p.e > 1.02] #& (p.a*(1-p.e**2))/(1-p.e*np.cos(p.f)) > 35]
                for p in particles_to_remove:
                    sim.remove(p)
                    removed_particles += 1
            progress.update(1)

        # Evaluate Energy difference of simulation and print ejected particles
        dE = abs((sim.energy() - E0))
        print(f"dE ({filename}) {dE}")
        print(f"Removed Particles ({filename}): {removed_particles}")

    # Passing function "integrate" in two separate threads. lambda functions for e_filter
    threads = [Thread(target=integrate, args=(sim, f, p, fn)) for sim, f, p, fn in zip(
        [sim1, sim2], [lambda e: e <= 0.4, lambda e: e > 0.4], progress_bar, filenames)]

    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()

#%% Actual integration
# Parameters needed to calculate Hill radius
a_scatterer = sim1.particles[1].a
M_scatterer = sim1.particles[1].m
M_sun = sim1.particles[0].m

# Generate uniform e(q) distribution of test-particles 
e_generate = np.arange(0, 0.95, 0.05)
q_generate = np.arange(a_scatterer - (5 * r_hill(a_scatterer, M_scatterer, M_sun)),
                       a_scatterer + (11 * r_hill(a_scatterer, M_scatterer, M_sun)), 0.5)

orbital_elements = generate_uniform_distribution(e_generate, q_generate)

print(f"added {len(orbital_elements)} test particles")

# Integrate test-population
custom_integration( orbital_elements[:, 1], # e
                    orbital_elements[:, 2], # a
                    1,                      # dt
                    2,                      # r_hill_crit
                    1e6,                    # tsim
                    ["archives/mercurius-testa-10e6.bin", "archives/mercurius-testb-10e6.bin"])
