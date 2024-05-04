import numpy as np

# Function for Hill-Radius
def r_hill(a, M_small_body, M_central_body):
    return a * (M_small_body / 3 * (M_central_body)) ** (1/3)

# Masses and Semi-major axis of Sun and Neptune
M_sun = 1
M_neptune = 5.151383772628957e-05
a_neptune = 30.237878368382898


# create arrays for e and q
e = np.arange(0, 1, .01)
q = np.arange(a_neptune - r_hill(a_neptune,M_neptune,M_sun), a_neptune + r_hill(a_neptune,M_neptune,M_sun), .01)

orbital_elements = []

for i in range(0, len(e)):
    for k in range(0, len(q)):
        a = q[k] / (1-e[i])
        orbital_elements.append([q[k], e[i], a])

orbital_elements = np.array(orbital_elements)
print(orbital_elements)
np.savetxt('test-populations/orbital_elements.txt', orbital_elements, delimiter=',')