#%% Libraries
import rebound
import numpy as np
import pandas as pd
import time
#%%
# import dataset of objects

df = pd.read_csv("test-populations/orbital_elements.txt", header=None, sep=",")

a = df[2]
e = df[1]
q = df[0]

'''#plot initial distribution
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
ax.plot(q, e, 'o')
plt.show()'''


#%% Simulation
start_time = time.perf_counter()

sim = rebound.Simulation()

# right units
sim.G = 6.674e-11
sim.units = ('yr', 'AU', 'Msun')

# Add Sun and Neptune
sim.add('Sun')
sim.add('Neptune')

# add object from csv
for i in range(len(a)):
    rand = np.random.random()*2*np.pi
    sim.add(a=a[i], e=e[i], Omega=0, omega=rand, f=rand)

sim.move_to_com()
E0 = sim.energy()

sim.save_to_file("archives/archive-t_1e6-mercurius-10_rhill_03.bin", interval=5e2, delete_file=True)

sim.integrator = "mercurius"

# settings for mercurius integrator
sim.dt = sim.particles[1].P * 0.1   # time-step of used WH-fast integrator
sim.ri_ias15.min_dt = 0.6           # minimal timestep of used ias15-integrator
sim.ri_mercurius.r_crit_hill = 4    # hill radius when integrator is switched

# Integrate
sim.integrate(1e6)

end_time = time.perf_counter()
print(f"Execution Time in second : {end_time - start_time}")
dE = abs((sim.energy() - E0)/E0)
print("dE=", dE) # energy can be used to estimate if simulation has been accurate (should be very low, about 1e-15 or so)
