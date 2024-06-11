import rebound
import numpy as np
from tqdm import tqdm


def r_hill(a, M_body, M_center):
    return a * (M_body / (3 * M_center)) ** (1 / 3)


sim = rebound.Simulation()
sim.units = ("yr", "AU", "Msun")
sim.add("Sun")
sim.add("Neptune")

n = 150

a_scatterer = sim.particles[1].a
M_scatterer = sim.particles[1].m
M_sun = sim.particles[0].m

a_obj = a_scatterer - (5 * r_hill(a_scatterer, M_scatterer, M_sun))
e_obj = 0

for i in range(n):
    rand1 = np.random.random() * 2 * np.pi
    rand2 = np.random.random() * 2 * np.pi
    rand3 = np.random.random() * 2 * np.pi
    sim.add(
        a=a_obj, e=e_obj, primary=sim.particles[0], Omega=rand1, omega=rand2, f=rand3
    )

sim.dt = 2
sim.ri_ias15.min_dt = 1e-6  # ensure that close encounters do not stall the integration
sim.ri_mercurius.r_crit_hill = 2
sim.move_to_com()
E0 = sim.energy()
# sim.save_to_file("archives/mercurius-e0-amin-10e6.bin", interval=1000, delete_file=True)

# Integration with progress bar
t_final = 1e6
t_step = 1  # Progress bar update interval

for t in tqdm(np.arange(0, t_final, t_step), desc="Integration Progress"):
    sim.integrate(t, exact_finish_time=0)

    if t % 1000 == 0:  # Adjust the frequency of removals
        sim.save_to_file("archives/mercurius-e0-amin-10e6.bin")

dE = abs((sim.energy() - E0))
print(f"dE = {dE}")
