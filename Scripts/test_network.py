#%%
## Import libraries
import numpy as np
import os
import matplotlib.pyplot as plt

from DS.solvers.diff_eqn_system import ReactionNetwork as rn

project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
simulations_dir = os.path.join(project_dir, "simulations")
results_dir = os.path.join(project_dir, "results")
figures_dir = os.path.join(project_dir, "figures")

#%%
# Initialize the entire system

# declare a time vector (time window)
t = np.arange(0,1000,0.01)
# Number of DOM/Carbon species:

dom_n = np.random.randint(4,10,1)[0]
bio_n = np.random.randint(2,10,1)[0]

# Initialize the same number of parameters:

vparams = np.random.uniform(0.2,10,bio_n*dom_n)
kparams = np.random.randint(1000,8000,bio_n*dom_n)
yparams = np.random.uniform(0.1,0.99,bio_n)
mparams = np.mean(vparams.reshape(dom_n, bio_n),axis=0)/10

# initial conditions
dom_initial = np.random.randint(500,1000,dom_n)
biomass_initial = np.random.randint(100,1000,bio_n)
x0 = np.append(dom_initial,biomass_initial)

# To enforce steady state, identify the simplest carbon compound
simplest_C = np.argsort(np.mean(kparams.reshape(dom_n, bio_n), axis=1))[0]
complex_C = np.argsort(np.mean(kparams.reshape(dom_n, bio_n), axis=1))[-1]
# Compute fresh carbon input commensurate to the mineralization of the simplest carbon compound
compute_carbon_input = dom_initial[simplest_C]*sum((1-yparams)*vparams.reshape(bio_n, dom_n)[:,simplest_C]*biomass_initial/(kparams.reshape(bio_n, dom_n)[:,simplest_C]+dom_initial[simplest_C]))

trial = rn(#maximum_capacity=300,
        #carbon_mol_bio = 10,
        carbon_num = dom_n,
        bio_num = bio_n,
        carbon_input = compute_carbon_input,
        sigmoid_coeff_stolpovsky = 0.01)
trial.set_rate_constants(vparams, kparams, yparams,mparams)
sim_array = trial.solve_network(x0, t)

C1 = sim_array[:,simplest_C]
C2 = sim_array[:,complex_C]
B = sim_array[:,dom_n:]
plt.plot(t, C1, linestyle = '-', label = "Simp_C")
plt.plot(t, C2, linestyle = '-', label = "Comp_C")
plt.plot (t,B, linestyle = '--', label = "Biomass")
