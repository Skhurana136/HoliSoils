# %%
## Import libraries
import numpy as np
import matplotlib.pyplot as plt

#from solvers import diff_eqn_system as diffsys
from DS.solvers.diff_eqn_system import ReactionNetwork as rn

# Number of DOM/Carbon species:

dom_n = np.random.randint(4,10,1)[0]
bio_n = np.random.randint(2,10,1)[0]
print("Carbon species: ", dom_n)
print("Biomass species: ", bio_n)

# Initialize the same number of parameters:

vparams = np.random.uniform(0.2,10,bio_n*dom_n)
kparams = np.random.randint(1000,8000,bio_n*dom_n)
yparams = np.random.uniform(0.1,0.99,bio_n)
mparams = np.mean(vparams.reshape(dom_n, bio_n),axis=0)/10

# initial conditions
dom_initial = np.random.randint(500,1000,dom_n)
biomass_initial = np.random.randint(100,1000,bio_n)
x0 = np.append(dom_initial,biomass_initial)

# declare a time vector (time window)
t = np.arange(0,10,0.01)

trial = rn(maximum_capacity=5000,
           #carbon_mol_bio = 10,
           carbon_num = dom_n,
           bio_num = bio_n,
           carbon_input = np.max(mparams),
           sigmoid_coeff_stolpovsky = 0.1)

trial.set_rate_constants(vparams, kparams, yparams,mparams)
x = trial.solve_network(x0, t)
## PLOT RESULTS

# Time series of DOM concentration profile with biomass concnetration profile
C = x[:,:dom_n]
B = x[:,dom_n:]
plt.plot(t,C, linestyle = "-", label = "DOM")
plt.plot(t,B, linestyle = "--", label = "Biomass")
plt.xlabel ("Time")
plt.ylabel ("Concentration")
#plt.legend()
# Shannon diversity vs carbon stock
S, DOC, TOC = trial.diversity_carbon ()
plt.figure()
plt.plot(S, DOC, label ="DOC")
plt.plot(S, TOC, label ="TOC")
plt.xlabel ("Shannon")
plt.ylabel ("Carbon")
plt.legend()
#%%
plt.figure()
plt.plot(S, label = "Shannon")
plt.legend()

plt.figure()
plt.plot(TOC, label = "TOC")
plt.ylim (bottom = 0)
plt.legend()

