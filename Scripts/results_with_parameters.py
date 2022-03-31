#%%
## Import libraries
import os
import pandas as pd
import numpy as np
import h5py as h

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
params_file = os.path.join(simulations_dir, "simulations.h5")
hr = h.File(params_file, mode = 'r')
# Load only those datasets that are captured in diversity_data dataset:
params_sim = {}
for sim in list(diversity_data.Sim_series.unique()):
    sim_data = hr[str(sim)]
    params = sim_data['parameters']
    ox_state = params['oxidation_state']
    exo_enz = params['exo_enzyme_rate']
    c_uptake = params['carbon_uptake']
    max_consumption = params['max_rate']
    half_k = params['half_saturation']
    yield_coef = params['yield_coefficient']
    initial_state = sim_data['initial_conditions']
    dom = initial_state['dom']
    bio = initial_state['biomass']
    dom_n = sim_data['species_number'][0]
    bio_n = sim_data['species_number'][1]
    params_sim[sim] = {'oxidation_state':ox_state, 'exo_enzyme_rate':exo_enz, 'carbon_uptake':c_uptake,
    'max_rate':max_consumption, 'half_saturation':half_k, 'yield_coefficient':yield_coef, 'dom':dom, 'bio':bio,
    'dom_n':dom_n, 'bio_n':bio_n}

#%%
# Take salient features of parameters into the diversity dataframe
for attr in list(params_sim[list(params_sim.keys())[0]].keys()):
    diversity_data['mean_'+attr] = diversity_data.apply (lambda row: np.mean(params_sim[row.Sim_series][attr]), axis=1)
#%%
# Correlation between any of the parameters in the dataframe?
corr = diversity_data.corr()
mask = np.triu(np.ones_like(corr, dtype=bool))
f,ax = plt.subplots(figsize=(11,11))
cmap = sns.diverging_palette(230,20, as_cmap=True)

sns.heatmap(corr, mask = mask, cmap = cmap, vmax = 0.9, center = 0, square=True)
plt.savefig(os.path.join(figures_dir, "correlation.png"), dpi = 300)

#%%
#parameter distributions:
for attr in list(params_sim[list(params_sim.keys())[0]].keys()):
    p = []
    for s in list(diversity_data.Sim_series.unique()):
        p.append(params_sim[s][attr][:].ravel().tolist())
    flat_list = [item for sublist in p for item in sublist]
    plt.figure()
    sns.kdeplot(data = flat_list)
    plt.xlabel(attr)
    plt.savefig(os.path.join(figures_dir, attr+".png"), dpi = 300)
#%%
diversity_data['DOC_biomass_ratio'] = diversity_data.DOC_end/diversity_data.Biomass_end
diversity_data['carbon_biomass_species_ratio'] = diversity_data.carbon_species/diversity_data.biomass_species
# Is biomass comparable to DOC in any of the simulations?
sns.scatterplot(x = 'carbon_biomass_species_ratio', y = 'DOC_biomass_ratio', hue = 'S_initial', data = diversity_data)
plt.xlabel("#carbon/#biomass")
plt.ylabel ('DOC/Biomass')
plt.savefig(os.path.join(figures_dir, "ratios_S_i.png"), dpi = 300)

#%%
# Does initial diversity and oxidation state impact the removal of DOC?
sns.scatterplot(x = 'S_initial', y = 'DOC_removal', hue = 'mean_oxidation_state', data = diversity_data)
plt.xlabel ('Initial Shannon diversity')
plt.ylabel("mean_oxidation_state")
plt.savefig(os.path.join(figures_dir, "S_i_ox_state.png"), dpi = 300)

#%%
sns.scatterplot(x = 'S_initial', y = 'DOC_removal', data = diversity_data)
plt.ylim(bottom = 0)
plt.xlabel ('Initial Shannon diversity')
plt.ylabel("mean_oxidation_state")
plt.savefig(os.path.join(figures_dir, "S_i_DOC_removal.png"), dpi = 300)
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