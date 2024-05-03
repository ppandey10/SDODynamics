import matplotlib.pyplot as plt
import rebound
import numpy as np

# constant(s)
P_NEPTUNE = 60189 / 365 # in yr


# import archive
sa = rebound.Simulationarchive("sun-neptune-scatterer-dt0f01-10e2P.bin")
# index of object that should be analysed
object_index = 2

##  Storing data of orbital elements at every timestep in arrays
# initialise arrays 
eccentricities = np.zeros(len(sa))
times = np.zeros(len(sa))
inclinations = np.zeros(len(sa))
omegas = np.zeros(len(sa))
Omegas = np.zeros(len(sa))
true = np.zeros(len(sa))
semi = np.zeros(len(sa))
dt = np.zeros(len(sa))
# fill arrays with data from archive
for i, sim in enumerate(sa):
    eccentricities[i] = sim.particles[object_index].orbit(primary=sim.particles[0]).e
    times[i] = sim.particles[object_index].orbit(primary=sim.particles[0]).P
    inclinations[i] = sim.particles[object_index].orbit(primary=sim.particles[0]).inc
    omegas[i] = sim.particles[object_index].orbit(primary=sim.particles[0]).omega
    Omegas[i] = sim.particles[object_index].orbit(primary=sim.particles[0]).Omega
    true[i] = sim.particles[object_index].orbit(primary=sim.particles[0]).f
    semi[i] = sim.particles[object_index].orbit(primary=sim.particles[0]).a
    dt[i] = sa[i].t

    #print(sim.particles[object_index])





#print(len(x_vals))
print(len(eccentricities))

for index in range(len(sa)):
    #print(sa[index])
    print(dt[index])



'''## plotting
# define data that should be plotted
# convert index of each time step in time [yr] 
x_vals = dt #* P_NEPTUNE 
y_vals_1 = eccentricities
y_vals_2 = inclinations
y_vals_3 = omegas
y_vals_4 = Omegas
y_vals_5 = semi

#plot
fig, axs= plt.subplots(nrows=2, ncols=2,figsize=(16,9), layout='constrained')
axs[0,0].scatter(x_vals, y_vals_1, s=2)
axs[0,0].set_xscale('log')
axs[0,0].set_ylabel('eccentricity')

axs[1,0].scatter(x_vals, y_vals_2, s=2)
axs[1,0].set_xscale('log')
axs[1,0].set_ylabel('inclination [deg]')
axs[1,0].set_xlabel('time [yr]')

axs[0,1].scatter(x_vals, y_vals_3, s=2)
axs[0,1].set_xscale('log')
axs[0,1].set_ylabel('argument of pericenter [deg]')

axs[1,1].scatter(x_vals, y_vals_5, s=2)
axs[1,1].set_xscale('log')
axs[1,1].set_ylabel('semi major axis [au/a.u.]')
axs[1,1].set_xlabel('time [yr]')

#plt.savefig('dt0f01P.pdf')

plt.show()
'''
