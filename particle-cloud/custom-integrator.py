#%% Libraries
import rebound
import numpy as np
from threading import Thread
from tqdm import tqdm

#%% Simulation
# initialise two simulations which can run parallel later on
sim1 = rebound.Simulation()
sim2 = rebound.Simulation()

## correction of units
sim1.units = ('yr', 'AU', 'Msun')
sim2.units = ('yr', 'AU', 'Msun')

## Add Sun and Neptune
sim1.add('Sun')
sim1.add('Neptune')
sim2.add('Sun')
sim2.add('Neptune')
print("added Sun and Neptune")


## Functions for test objects and Integration

# Function for Hill-Radius
def r_hill(a_small_body, M_small_body, M_central_body):
    return a_small_body * (M_small_body / (3 * M_central_body)) ** (1 / 3)


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


# custom function for Integration

def custom_integration(e, a, integrator, dt_min_ias15, r_hill_crit, t_sim, filename1, filename2):
    # specify integrator
    integrators = ["ias15", "mercurius"]
    integrators = integrators[integrator]
    print("Using ", integrators, "integrator")
    sim1.integrator = integrators
    sim2.integrator = integrators

    # Define two functions for integration to pass them in threads later on.
    def first_integration(progress_bar):

        integrators = "mercurius"

        # add particles
        for ind in range(0, len(e)):
            if e[ind] <= 0.4:
                rand = np.random.random() * 2 * np.pi
                sim1.add(a=a[ind], e=e[ind], Omega=0, omega=rand, f=rand)
        print("N particles (sim1): ", sim1.N)

        # different options depending on choice of integrator (can be removed maybe if mercurius is chosen anyway)
        if integrators == "mercurius":
            #sim1.ri_ias15.min_dt = dt_min_ias15
            sim1.dt = dt_min_ias15  # * sim1.particles[1].P / (2 * np.pi)
            sim1.ri_mercurius.r_crit_hill = r_hill_crit

        if integrators == "ias15":
            sim1.dt = dt_min_ias15

        sim1.move_to_com()
        sim1.save_to_file(filename1, interval=1e3, delete_file=True)
        sim1.integrate(t_sim)

        #sim1.move_to_com()
        E0 = sim1.energy()
        time_steps = np.arange(0, t_sim, 1)

'''        # integrate of each timestep
        for ts in time_steps:
            sim1.integrate(ts, exact_finish_time=0)
            #check after 1000 timesteps if objects have been ejected
            if ts % 500 == 0:
                sim1.save_to_file(filename1)
                removed_particles = 0
                for k in range(2, sim1.N):
                    if sim1.particles[k].e > 1.1:
                        sim1.remove(k)
                        removed_particles += 1
            progress_bar[0].update(1)

        dE = abs((sim1.energy() - E0) / E0)
        print("dE (sim1) ", dE)
        print("Removed Particles (sim 1): ", removed_particles)'''

    def second_integration(progress_bar):
        integrators = "mercurius"

        # add particles
        for ind in range(0, len(e)):
            if e[ind] > 0.4:
                rand = np.random.random() * 2 * np.pi
                sim2.add(a=a[ind], e=e[ind], Omega=0, omega=rand, f=rand)
        print("N particles (sim2): ", sim2.N)

        if integrators == "mercurius":
            #sim2.ri_ias15.min_dt =  dt_min_ias15
            sim2.dt = 2 * dt_min_ias15  #  * sim2.particles[1].P / (2 * np.pi) # timesteps are 5-times larger then for low eccentricities
            sim2.ri_mercurius.r_crit_hill = r_hill_crit

        if integrators == "ias15":
            sim2.dt = 5 * dt_min_ias15

        sim2.move_to_com()
        sim2.save_to_file(filename2, interval=1e3, delete_file=True)
        sim2.integrate(t_sim)

        #sim2.move_to_com()
        E0 = sim2.energy()

'''        for ts in time_steps:
            sim2.integrate(ts, exact_finish_time=0)
            if ts % 500 == 0:
                sim2.save_to_file(filename2)
                removed_particles = 0
                for k in range(2, sim2.N):
                    if sim2.particles[k].e > 1.1:
                        sim2.remove(k)
                        removed_particles += 1
                    # TODO: particles are not removed (maybe because only initial parameters accessed?x)

            progress_bar[1].update(1)

        dE = abs((sim2.energy() - E0) / E0)
        print("dE (sim2) ", dE)
        print("Removed Particles (sim 2): ", removed_particles)'''

    # start both integrations in two different threads
    progress_bar = [tqdm(total=t_sim, desc="Thread 1 Progress", position=0),
                    tqdm(total=t_sim, desc="Thread 2 Progress", position=0)]

    threadA = Thread(target=first_integration, args=(progress_bar,))
    threadB = Thread(target=second_integration, args=(progress_bar,))

    threadA.start()
    threadB.start()

    threadA.join()
    threadB.join()

    # Uncomment this and comment thread.start() + thread.join() to let threads run after each other
    # threadA.run()
    # threadB.run()


# Integrate
# necessary for r_hill
a_neptune = sim1.particles[1].a
M_neptune = sim1.particles[1].m
M_sun = sim1.particles[0].m

# create arrays for e and q with Hill-radius as reference
e_generate = np.arange(0, 0.95, 0.05)
q_generate = np.arange(a_neptune - (5 * r_hill(a_neptune, M_neptune, M_sun)),
                       a_neptune + (20 * r_hill(a_neptune, M_neptune, M_sun)), 1)

orbital_elements = generate_uniform_distribution(e_generate, q_generate)

peri = orbital_elements[:, 0]
ecc = orbital_elements[:, 1]
semi = orbital_elements[:, 2]

print(f"added {len(orbital_elements)} test particles")



custom_integration(ecc,
                   semi,
                   1,
                   50,
                   2,
                   1e6,
                   "archives/test1.bin",
                   "archives/test2.bin")
