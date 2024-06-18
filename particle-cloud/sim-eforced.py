import rebound
import numpy as np
from tqdm import tqdm


def r_hill(a, M_body, M_center):
    return a * (M_body / (3 * M_center)) ** (1 / 3)


# Function for distribution
def generate_uniform_distribution(e_, q_):
    # takes in two arrays "e" and "q" and outputs numpy array with "q", "e" and "a"
    o_e = []
    for i in range(0, len(e_)):
        for k in range(0, len(q_)):
            a = q_[k] / (1 - e_[i])
            o_e.append([q_[k], e_[i], a])
    o_e = np.array(o_e)
    return o_e


sim = rebound.Simulation()
sim.units = ("yr", "AU", "Msun")
sim.add("Sun")
sim.add("Neptune")

a_scatterer = sim.particles[1].a
M_scatterer = sim.particles[1].m
M_sun = sim.particles[0].m

# create arrays for e and q with Hill-radius as reference
e_generate = np.arange(0, 0.95, 0.05)
q_generate = np.arange(
    a_scatterer - (5 * r_hill(a_scatterer, M_scatterer, M_sun)),
    a_scatterer + (20 * r_hill(a_scatterer, M_scatterer, M_sun)),
    1,
)

orbital_elements = generate_uniform_distribution(e_generate, q_generate)

peri = orbital_elements[:, 0]
ecc = orbital_elements[:, 1]
semi = orbital_elements[:, 2]

n = 150

a_obj = a_scatterer - (5 * r_hill(a_scatterer, M_scatterer, M_sun))
e_obj = 0

for i in range(0, len(q_generate)):
    rand1 = np.random.random() * 2 * np.pi
    rand2 = np.random.random() * 2 * np.pi
    rand3 = np.random.random() * 2 * np.pi
    sim.add(a=semi[i], e=ecc[i], Omega=rand1, omega=rand2, f=rand3)

sim.dt = 5
sim.ri_ias15.min_dt = 5e-4  # ensure that close encounters do not stall the integration
sim.ri_mercurius.r_crit_hill = 2
sim.move_to_com()
E0 = sim.energy()

# Integration with progress bar
t_final = 1e8
t_step = 1  # Progress bar update interval

for t in tqdm(np.arange(0, t_final, t_step), desc="Integration Progress"):
    sim.integrate(t, exact_finish_time=0)

    if t % 1000 == 0:  # Adjust the frequency of removals
        sim.save_to_file("archives/mercurius_tsim_1e9_dtias15_5e-4_dt_5_rhill_2.bin")

dE = abs((sim.energy() - E0))
print(f"dE = {dE}")
