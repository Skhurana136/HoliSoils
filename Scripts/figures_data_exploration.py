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
details_subfolder = 'paras_adjust_v2'
simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
results_dir = os.path.join(project_dir, "results", details_subfolder)
figures_dir = os.path.join(project_dir, "figures", details_subfolder)

filename = os.path.join(results_dir, "diversity_data.pkl")
diversity_data = pd.read_pickle(filename)

#%%
# Additional data processing
diversity_data['DOC_removal'] = (1 - diversity_data.DOC_end/diversity_data.DOC_input) * 100
diversity_data['carbon_biomass'] = diversity_data.carbon_species * diversity_data.biomass_species
print(diversity_data.shape)
print(diversity_data.columns)

#%%
# Does initial diversity impact the removal of DOC?
sns.scatterplot(x = 'S_initial', y = 'DOC_removal', data = diversity_data)
plt.xlabel ('Initial Shannon diversity')
plt.ylabel("DOC removal (%)")
plt.savefig(os.path.join(figures_dir, "S_i_DOC_norm.png"), dpi = 300)


#%%
# Is there a relationship between number of biomass species and Shannon?
sns.scatterplot(x = 'S_initial', y = 'S_at_max_DOC', data = diversity_data)
plt.xlabel ('Shannon_initial')
plt.ylabel("Shannon_at_max_DOC")
plt.savefig(os.path.join(figures_dir, "S_i_at_max_DOC.png"), dpi = 300)

#%%
# Is there a relationship between number of biomass species and Shannon?
sns.scatterplot(x = 'S_initial', y = 'S_max', data = diversity_data)
plt.xlabel ('Shannon_initial')
plt.ylabel("Shannon_max")
plt.savefig(os.path.join(figures_dir, "S_i_S_max.png"), dpi = 300)

#%%
# Is there a relationship between number of biomass species and Shannon?
sns.scatterplot(x = 'S_initial', y = 'DOC_end', data = diversity_data)
plt.xlabel ('Shannon_initial')
plt.ylabel("DOC_end")
plt.savefig(os.path.join(figures_dir, "S_i_DOC_end.png"), dpi = 300)

#%%
# Is there a relationship between number of biomass species and Shannon?
sns.scatterplot(x = 'S_initial', y = 'Max_DOC', data = diversity_data)
plt.xlabel ('Shannon_initial')
plt.ylabel("DOC_max")
plt.savefig(os.path.join(figures_dir, "S_i_DOC_max.png"), dpi = 300)

#%%
# Is there a relationship between number of biomass species and Shannon?
sns.scatterplot(x = 'carbon_species', y = 'S_max', data = diversity_data)
plt.xlabel ('Carbon_number')
plt.ylabel("Max Shannon diversity index")
plt.savefig(os.path.join(figures_dir, "dom_diversity_shannon.png"), dpi = 300)

#%%
# Is there a relationship between the combination of number of biomass and carbon species and Shannon?

sns.scatterplot(x = 'carbon_biomass', y = 'S_end', data = diversity_data)
plt.xlabel ('Carbon x Biomass number')
plt.ylabel("Max Shannon diversity index")
plt.savefig(os.path.join(figures_dir, "domxbiomass_diversity_shannon.png"), dpi = 300)

#%%
# Is diversity between initial and end states changing?
sns.scatterplot(x = 'S_initial', y = 'S_end', data = diversity_data)
plt.xlabel ('Initial diversity')
plt.ylabel("Final diversity")
plt.savefig(os.path.join(figures_dir, "S_i_S_end.png"), dpi = 300)