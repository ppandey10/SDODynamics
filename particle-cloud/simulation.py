#%% Libraries
import rebound
import numpy as np
import polars as pl
import time
#%%
# import dataset of objects

df = pl.read_csv("test-populations/orbital_elements-5_rhill.csv")

a = df["semi"]
e = df["ecc"]
q = df["peri"]


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

sim.save_to_file("archives/archive-t_1e6-ias15-dt_1-5_rhill.bin",interval=1e2,delete_file=True)

sim.integrator = "ias15"
sim.ri_ias15.min_dt = 1
sim.integrate(1e6)
end_time = time.perf_counter()
print(f"Execution Time in second : {end_time - start_time}")
dE = abs((sim.energy() - E0)/E0)
print("dE=", dE) # energy can be used to estimate if simulation has been accurate (should be very low, about 1e-15 or so)

