import rebound
import numpy as np

# Functions for Hill-Radius and distribution
def r_hill(a, M_body, M_center):
    return a * (M_body / (3 * M_center)) ** (1 / 3)


def generate_uniform_distribution(e, q):
    results = []
    for ei in e:
        for qi in q:
            results.append([qi, ei, qi / (1 - ei)])
    return np.array(results)

sim = rebound.Simulation()
sim.units = ("yr", "AU", "Msun")
sim.add("Sun")
sim.add("Neptune")

# Parameters needed to calculate Hill radius
a_scatterer = sim.particles[1].a
M_scatterer = sim.particles[1].m
M_sun = sim.particles[0].m

# Generate uniform e(q) distribution of test-particles
e_generate = np.linspace(0, 1, 20, endpoint=False)
q_generate = np.linspace(30, 42, 20)

# for uniform perihelion distribution
orbital_elements = generate_uniform_distribution(e_generate, q_generate)

for i in range(400):
    rand1 = np.random.random() * 2 * np.pi
    rand2 = np.random.random() * 2 * np.pi
    rand3 = np.random.random() * 2 * np.pi
    sim.add(a=orbital_elements[i,2], e=orbital_elements[i,1], Omega=rand1, omega=rand2, f=rand3)
    
sim.integrator = "mercurius"
sim.dt = 5
sim.ri_ias15.min_dt = 5e-4  # ensure that close encounters do not stall the integration
sim.ri_mercurius.r_crit_hill = 2
sim.move_to_com()
E0 = sim.energy()

t_final = 1e7

sim.save_to_file("no_inc-final_1e7.bin", interval=1000)

sim.integrate(t_final)

sim.status()

dE = abs((sim.energy() - E0))
print(f"dE = {dE}")


'''
### PLOTTING #####################

import matplotlib.pyplot as plt

sim = rebound.Simulationarchive("archives/tmp-test-final.bin")

initial_sim = sim[0]  # First snapshot
final_sim = sim[len(sim)-2]   # Last snapshot (10 Myr later)

print(initial_sim)
print(final_sim)

# Extract the proper elements
particles_initial = initial_sim.particles
particles_final = final_sim.particles

# Initialize lists to store proper elements
proper_elements_initial = []
proper_elements_final = []

# Iterate through particles to extract (a, e) at initial and final times
for i in range(2, len(particles_initial)):  # first particle is the Sun
    p_initial = particles_initial[i]
    p_final = particles_final[i]
    
    a_initial, e_initial = p_initial.a, p_initial.e
    a_final, e_final = p_final.a, p_final.e
    
    proper_elements_initial.append((a_initial, e_initial))
    proper_elements_final.append((a_final, e_final))

proper_elements_initial = np.array(proper_elements_initial)
proper_elements_final = np.array(proper_elements_final)

# Remove particles with e >= 1 and count how many were removed
mask = proper_elements_final[:, 1] <= 1
removed_particles_count = len(proper_elements_final) - np.sum(mask)
proper_elements_final = proper_elements_final[mask]
proper_elements_initial = proper_elements_initial[mask]


print(f"Number of particles removed: {removed_particles_count}")

print(len(proper_elements_initial[:,0]))

plt.scatter(proper_elements_initial[:,1], proper_elements_initial[:,0] * (1 - proper_elements_initial[:,1]))
plt.show() 
'''
