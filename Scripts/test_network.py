#%%
## Import libraries
import numpy as np
import os
import matplotlib.pyplot as plt

from DS.solvers.diff_eqn_system import ReactionNetwork as rn
from DS.analyses.steady_state import diversity_carbon, normalize_carbon

print ("All libraries loaded.")

project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
simulations_dir = os.path.join(project_dir, "simulations")
results_dir = os.path.join(project_dir, "results")
figures_dir = os.path.join(project_dir, "figures")

print ("Project directories set up.")

#%%
# Initialize the entire system
# declare a time vector (time window)
t_span = [0,100000]
t = np.arange(t_span[0], t_span [1],0.01)

# Number of DOM/Carbon species:

dom_n = np.random.randint(5,10,1)[0]
bio_n = np.random.randint(4,10,1)[0]
print(dom_n, bio_n)
# Initialize the same number of parameters:

enzparams = np.random.uniform(0.2, 0.4, bio_n)
vparams = np.random.normal(0.8, 0.001, bio_n*dom_n)#np.random.uniform(0.001,15,bio_n*dom_n)
kparams = np.random.normal(600, 200, bio_n*dom_n)#np.random.randint(1000,8000,bio_n*dom_n)
yparams = np.random.normal(0.5, 0.01, bio_n)#np.random.uniform(0.1,0.99,bio_n)
mparams = np.mean(vparams.reshape(dom_n, bio_n),axis=0)/5

# initial conditions
dom_initial = np.random.normal(400, 100, dom_n)#np.random.randint(500,1000,dom_n)
biomass_initial = np.random.normal(400, 50, bio_n)#np.random.randint(100,1000,bio_n)
x0 = np.append(dom_initial,biomass_initial)

# To enforce steady state, identify the simplest carbon compound
simplest_C = np.argsort(np.mean(kparams.reshape(dom_n, bio_n), axis=1))[-1]
complex_C = np.argsort(np.mean(kparams.reshape(dom_n, bio_n), axis=1))[0]
# Compute fresh carbon input commensurate to the mineralization of the simplest carbon compound
compute_carbon_input = 0#[100,50,10,0,50]

trial = rn(maximum_capacity=5,
        #carbon_mol_bio = 10,
        carbon_num = dom_n,
        bio_num = bio_n,
        carbon_input = compute_carbon_input,
        sigmoid_coeff_stolpovsky = 0.01,
        enzyme_production_rate_constant = enzparams[0],
        efficiency_bio_uptake = enzparams[1],
        necromass_distribution="notequal")
trial.set_rate_constants(vparams, kparams, yparams,mparams)
trial.identify_components_natures()
solution = trial.solve_network(x0, t_span, t)

tim = solution.t
sim_array = solution.y
C1 = sim_array[trial.most_labile_c,:]
C2 = sim_array[trial.least_labile_c,:]
C = sim_array[:dom_n,:].T
B = sim_array[dom_n:,:].T
#plt.plot(tim, C1, linestyle = '-', label = "Simp_C")
#plt.plot(tim, C2, linestyle = '-', label = "Comp_C")
plt.figure()
plt.plot (tim,B, linestyle = '--')#, label = "Bacteria")
plt.plot (tim,C, linestyle = '-')#, label = "Bacteria")
#plt.ylim(bottom=0.0)


#plt.figure()
#for i in list(range(dom_n)):
#    plt.plot(C[i,:], linestyle = '-', label = i)
#plt.legend()

#%%# Shannon diversity
S, DOC, TOC = diversity_carbon(sim_array, dom_n, bio_n)

#%%
plt.scatter(S, TOC)
#plt.yscale("log")

#%%
plt.hist(yparams, bins = 20)