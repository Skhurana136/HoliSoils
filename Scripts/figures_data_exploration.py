#%%
## Import libraries
import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns

## LOAD RESULTS
#project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
project_dir = "C:/Users/swkh9804/Documents/Projects/HoliSoils/data"
details_subfolder = 'carbon_18' + '_ip_0'
simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
results_dir = os.path.join(project_dir, "results", details_subfolder)
figures_dir = os.path.join(project_dir, "figures", details_subfolder)

filename = os.path.join(results_dir, "diversity_data.pkl")
diversity_data = pd.read_pickle(filename)

#%%
# Additional data processing
diversity_data['DOC_removal'] = (1 - diversity_data.DOC_end/diversity_data.DOC_input) * 100
diversity_data['carbon_biomass'] = diversity_data.carbon_species * diversity_data.biomass_species
diversity_data['t_50_days']  = diversity_data.T_50/100
diversity_data['t_50_b1_days']  = diversity_data.T_50_B1/100
diversity_data['t_50_ratio'] = diversity_data.T_50/diversity_data.T_50_B1
print(diversity_data.shape)
print(diversity_data.columns)

#%%
# Does initial diversity impact the removal of DOC?
sns.scatterplot(x = 'S_initial', y = 'DOC_removal',
data = diversity_data)
plt.xlabel ('Initial Shannon diversity')
plt.ylabel("DOC removal (%)")
plt.savefig(os.path.join(figures_dir, "S_i_DOC_norm.png"), dpi = 300)

#%%
sns.scatterplot(x = 'S_initial', y = 't_50_days', data = diversity_data)
plt.xlabel ('Initial Shannon diversity')
plt.ylabel("t50_days")
plt.savefig(os.path.join(figures_dir, "S_i_t_50.png"), dpi = 300)

#%%
sns.scatterplot(x = 'S_initial', y = 't_50_ratio', data = diversity_data)
plt.xlabel ('Initial Shannon diversity')
plt.ylabel("Ratio_t_50")
plt.savefig(os.path.join(figures_dir, "S_i_ratio_t_50.png"), dpi = 300)

#%%
# Is there a relationship between initial and max Shannon?
sns.scatterplot(x = 'S_initial', y = 'DOC_T50_B1', data = diversity_data)
plt.xlabel ('Shannon_initial')
plt.ylabel("DOC_baselinet50")
plt.savefig(os.path.join(figures_dir, "S_i_DOC_t50_b1.png"), dpi = 300)

#%%
# Is there a relationship between initial and max Shannon?
sns.scatterplot(x = 'S_initial', y = 'Biomass_t_50_b1', data = diversity_data)
plt.xlabel ('Shannon_initial')
plt.ylabel("Biomass_baselinet50")
plt.savefig(os.path.join(figures_dir, "S_i_Biomass_t50_b1.png"), dpi = 300)

#%%
# Is there a relationship between number of biomass species and Shannon?
sns.scatterplot(x = 'carbon_species', y = 'S_max', data = diversity_data)
plt.xlabel ('#Carbon species')
plt.ylabel("DOC_end")
plt.savefig(os.path.join(figures_dir, "dom_n_S_max.png"), dpi = 300)

#%%
# Is there a relationship between number of biomass species and Shannon?
sns.scatterplot(x = 'biomass_species', y = 'S_max', data = diversity_data)
plt.xlabel ('#Biomass species')
plt.ylabel("Shannon_max")
plt.savefig(os.path.join(figures_dir, "bio_n_S_max.png"), dpi = 300)