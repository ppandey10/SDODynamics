import numpy as np
import rebound

def LoadOrbitalElements (Filepath1, Filepath2):

	'''
	In: 
	Filepath1: Path to File from first simulation 
	Filepath2: Path to File from second simulation 

	Out: Matrix, 	first row -> e_initial, a_initial, q_initial
					second row -> e_final, a_final, q_final

	'''

	sa1 = rebound.Simulationarchive(Filepath1)
	sa2 = rebound.Simulationarchive(Filepath2)
	
	#%% load initial orbital elements
	sim1_initial = sa1[-len(sa1)]
	sim2_initial = sa2[-len(sa2)]
	
	# last particle index
	i_last1 = sa1[1].N
	i_last2 = sa2[1].N
	
	# create arrays for initial a, e and q
	e_initial1 = e_initial2 = a_initial1 = a_initial2 = np.empty([])
	
	for i in range(2, i_last1):
	    e_initial1 = np.append(e_initial1, sim1_initial.particles[i].e)
	    a_initial1 = np.append(a_initial1, sim1_initial.particles[i].a)
	
	for k in range(2, i_last2):
	    e_initial2 = np.append(e_initial2, sim2_initial.particles[k].e)
	    a_initial2 = np.append(a_initial2, sim2_initial.particles[k].a)
	
	q_initial1 = a_initial1 * (1 - e_initial1)
	q_initial2 = a_initial2 * (1 - e_initial2)
	
	q_initial = np.append(q_initial1, q_initial2)
	a_initial = np.append(a_initial1, a_initial2)
	e_initial = np.append(e_initial1, e_initial2)
	
	# set edges of histogram
	#q_initial = q_initial[1:]
	x_min = 26.328623 #min(q_initial)
	x_max = max(q_initial)
	
	#xedges = np.linspace(x_min, x_max, 4)  # 4 edges to create 3 bins
	
	# Calculate y-axis edges
	y_min = 0.000207#min(e_initial)
	y_max = max(e_initial)
	
	#yedges = np.linspace(y_min, y_max, 4)  # 4 edges to create 3 bins
	
	#%% load final orbital parameters
	sim_final1 = sa1[-1]
	sim_final2 = sa2[-1]
	
	i_last1 = sa1[-1].N
	i_last2 = sa2[-1].N
	
	e_final1 = []
	a_final1 = []
	e_final2 = []
	a_final2 = []
	
	for i in range(2, i_last1):
	    e_final1.append(sim_final1.particles[i].e)
	    a_final1.append(sim_final1.particles[i].a)
	
	for k in range(2, i_last2):
	    e_final2.append(sim_final2.particles[k].e)
	    a_final2.append(sim_final2.particles[k].a)
	
	e_final1 = np.array(e_final1)
	a_final1 = np.array(a_final1)
	q_final1 = a_final1 * (1 - e_final1)
	
	e_final2 = np.array(e_final2)
	a_final2 = np.array(a_final2)
	q_final2 = a_final2 * (1 - e_final2)
	
	# add arrays for both simulations
	e_final = np.append(e_final1, e_final2)
	a_final = np.append(a_final1, a_final2)
	q_final = a_final * (1 - e_final)

	result = {
		"sa1": sa1,
		"sa2": sa2,
		"e_initial": e_initial,
		"a_initial": a_initial,
		"q_initial": q_initial,
		"e_final": e_final,
		"a_final": a_final,
		"q_final": q_final,
	}
	
	return result
