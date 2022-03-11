#%%
## Import libraries
import os
import pandas as pd
import numpy as np
import h5py

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns

## LOAD RESULTS
#project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
project_dir = "C:/Users/swkh9804/Documents/Projects/HoliSoils/data"
results_dir = os.path.join(project_dir, "simulations")
output_dir = os.path.join(project_dir, "results")
figures_dir = os.path.join(project_dir, "figures")

filename = os.path.join(output_dir, "diversity_data.pkl")
diversity_data = pd.read_pickle(filename)

#%%
# Is there a relationship between number of biomass species and Shannon?
sns.scatterplot(x = 'biomass_species', y = 'Shannon', data = diversity_data)
plt.xlabel ('Biomass_number')
plt.ylabel("Shannon diversity index")
plt.savefig(os.path.join(figures_dir, "biomass_diversity_shannon.png"), dpi = 300)

#%%
# Is there a relationship between number of biomass species and Shannon?
sns.scatterplot(x = 'carbon_species', y = 'Shannon', data = diversity_data)
plt.xlabel ('Carbon_number')
plt.ylabel("Shannon diversity index")
plt.savefig(os.path.join(figures_dir, "dom_diversity_shannon.png"), dpi = 300)

#%%
# Is there a relationship between the combination of number of biomass and carbon species and Shannon?
diversity_data['carbon_biomass'] = diversity_data.carbon_species * diversity_data.biomass_species
sns.scatterplot(x = 'carbon_biomass', y = 'Shannon', data = diversity_data)
plt.xlabel ('Carbon x Biomass number')
plt.ylabel("Shannon diversity index")
plt.savefig(os.path.join(figures_dir, "domxbiomass_diversity_shannon.png"), dpi = 300)

#%%
plt.figure()
sns.scatterplot(x = 'Shannon',y='DOC', data = diversity_data, size = 'carbon_species',
hue = 'biomass_species')
plt.xlabel ('Shannon diversity index')
plt.ylabel("DOC")
plt.legend(ncol = 2)
plt.savefig(os.path.join(figures_dir, "diversity_conc_steady_state.png"), dpi = 300)
#%%
plt.figure()
sns.scatterplot(x = 'carbon_biomass',y='DOC', data = diversity_data, size = 'Shannon')
plt.xlabel ('Carbon x Biomass number')
plt.ylabel("DOC")
plt.legend(ncol = 2)
plt.savefig(os.path.join(figures_dir, "diversity_cxb_conc_steady_state.png"), dpi = 300)