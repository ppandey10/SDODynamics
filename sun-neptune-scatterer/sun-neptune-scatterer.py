import rebound 
import numpy as np
 

## initialise simulation
sim = rebound.Simulation()

## scaling factor for final time
def scale_time(P,t):
    return (t * 2*np.pi) / P

## add objects
sim.add('Sun')
sim.add('Neptune')
# object that should be scattered
sim.add(a=20, e=0.7, inc=np.pi/3)
# move to center of mass
sim.move_to_com()
# define timesteps as fractions of the orbital period of Neptune (1/2pi) is conversion factor to get years)
#sim.dt = sim.particles[1].P / (2*np.pi) * 0.01 # fixed time steps don't work with ias15 but with whfast

# chose integrator
sim.integrator = "ias15" #"whfast" 
#sim.ri_ias15.safe_mode = 1

# save data for each time step in file
#sim.save_to_file("sun-neptune-scatterer-dt0f01-10e2P.bin",step=1,delete_file=True)

sim.save_to_file("sun-neptune-scatterer-test.bin",step=1,delete_file=True)

## integrate
# orbital period of reference body in yr
REF_ORBITAL_PERIOD = sim.particles[1].P / (2*np.pi)

# integrate
#sim.integrate(1e4 * sim.particles[1].P / (2*np.pi) * 0.01, exact_finish_time=0) # scaling the final time with orbital period

sim.integrate(10 * REF_ORBITAL_PERIOD, exact_finish_time=0) # scaling the final time with orbital period

print(sim.particles[1].P / (2*np.pi) )




## Just for testing
sa = rebound.Simulationarchive("sun-neptune-scatterer-test.bin")

for i in range(len(sa)):
    print(sa[i].t)