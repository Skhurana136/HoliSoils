## Import libraries
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py

from DS.solvers.diff_eqn_system import diversity_carbon

## LOAD RESULTS
project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
results_dir = os.path.join(project_dir, "simulations")
hw = h5py.File(os.path.join(results_dir,"simulations.h5"), mode = 'r')

# Load a particular dataset:
data = hw['0']

# Identify number of dom and microbial species:
dom_n, bio_n = data['species_number']
x = data['solution']

# Time series of DOM concentration profile with biomass concnetration profile
C = x['dom']
B = x['biomass']
t = np.shape(C)[0]
plt.plot(C, linestyle = "-", label = "DOM")
plt.plot(B, linestyle = "--", label = "Biomass")
plt.xlabel ("Time")
plt.ylabel ("Concentration")
#plt.legend()

# Shannon diversity vs carbon stock
data = np.append(C,B,axis=1)

S, DOC, TOC = diversity_carbon(data, dom_n, bio_n)
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

