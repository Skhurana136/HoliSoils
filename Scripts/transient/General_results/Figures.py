#%%
## Import libraries
import os
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

## LOAD RESULTS
## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient", "gen_spec_lognorm")
simulations_dir = os.path.join(project_dir, "simulations")
results_dir = os.path.join(project_dir, "results")
figures_dir = os.path.join(project_dir, "figures")

null_prefixes = ["expo_", "gamma_", "1c_", ""]
adapt_prefixes = ["expo_", "gamma_", "1c_", "gen_"]

null_files = []
adapt_files = []
for p,a in zip(null_prefixes, adapt_prefixes):
    filename = os.path.join(results_dir, p + "null_combined_dataset.pkl")
    n_data = pd.read_pickle(filename)
    null_files.append(n_data)
    filename = os.path.join(results_dir, a + "adaptation_combined_dataset.pkl")
    a_data = pd.read_pickle(filename)
    adapt_files.append(a_data)

null_data = pd.concat(null_files)
adapt_data = pd.concat(adapt_files)

null_data = null_data.drop_duplicates()
null_data['DOC_initial_int'] = round(null_data.DOC_initial, -3)
null_data['S_initial_int'] = round(null_data.S_initial, 1)
null_data['ratio_t_50'] = null_data.T_50/null_data.T_50_B1
adapt_data = adapt_data.drop_duplicates()
adapt_data['DOC_initial_int'] = round(adapt_data.DOC_initial, -3)
adapt_data['S_initial_int'] = round(adapt_data.S_initial, 1)
adapt_data['ratio_t_50'] = adapt_data.T_50/adapt_data.T_50_B1

a_n_data = [null_data, adapt_data]
all_data = pd.concat(a_n_data)

#%%
diversity_data = adapt_data
#%%
# Distribution of number of carbon compounds and biomass species that we tested.
sns.kdeplot(x = "carbon_species", data = diversity_data, label = "#carbon")
sns.kdeplot(x = "biomass_species", data = diversity_data, label = "#biomass")
plt.legend()

#%%
sns.kdeplot(x = "DOC_removal", data = diversity_data)
#%%
sns.kdeplot(x = "ratio_t_50", data = diversity_data)
#%%
diversity_data['S_initial_int'] = round(diversity_data.S_initial,1)
print(diversity_data.S_initial_int.unique())
#%%
for s in list(diversity_data.S_initial_int.unique()):
    sns.kdeplot(x = 100 - diversity_data[diversity_data.S_initial_int==s]["DOC_removal"], data = diversity_data[diversity_data.S_initial_int==s], label=s)
plt.legend()

#%%
for s in list(diversity_data.S_initial_int.unique()):
    sns.kdeplot(x = diversity_data[diversity_data.S_initial_int==s]["ratio_t_50"], data = diversity_data[diversity_data.S_initial_int==s], label=s)
plt.legend()
#%%
sns.kdeplot(x = "T_50", data = diversity_data, hue = "carbon_species")
plt.xscale("log")
plt.xlabel("Time for 37% loss (days)")
plt.tight_layout()
plt.savefig(os.path.join(figures_dir, "adapt_t_63_distribution.png"), dpi = 300)
#%%
sns.kdeplot(x = "S_initial", data = diversity_data, label = "Initial S")
sns.kdeplot(x = "S_max", data = diversity_data, label = "Max S")
sns.kdeplot(x = "S_t_50", data = diversity_data, label = "S_at_t50")
sns.kdeplot(x = "S_t_50_b1", data = diversity_data, label = "S_at_t50_b1")
sns.kdeplot(x = "S_end", data = diversity_data, label = "S_end")
plt.xlabel("Shannon diversity index")
plt.legend()
#%%
sns.kdeplot(x = "Biomass_initial", data = diversity_data, label = "Initial Biomass")
sns.kdeplot(x = "Biomass_max", data = diversity_data, label = "Max Biomass")
sns.kdeplot(x = "Biomass_t_50", data = diversity_data, label = "Biomass_at_t50")
sns.kdeplot(x = "Biomass_t_50_b1", data = diversity_data, label = "Biomass_at_t50_b1")
sns.kdeplot(x = "Biomass_end", data = diversity_data, label = "Biomass_end")
plt.xlabel("Biomass")
plt.legend()
#%%
diversity_data['ratio_S_i_end'] = diversity_data.S_initial/diversity_data.S_end
plt.hist('ratio_S_i_end', data = diversity_data, label = "S_end")
plt.xlabel("Shannon diversity index: initial/final")
#%%
diversity_data['ratio_S_i_max'] = diversity_data.S_initial/diversity_data.S_max
plt.hist('ratio_S_i_max', data = diversity_data, label = "S_end")
plt.xlabel("Shannon diversity index: initial/max")
#%%
diversity_data['ratio_biomass_i_max'] = diversity_data.Biomass_initial/diversity_data.Biomass_max
plt.hist('ratio_biomass_i_max', data = diversity_data, label = "S_end")
plt.xlabel("Biomass: initial/max")
#%%
diversity_data['ratio_biomass_i_end'] = diversity_data.Biomass_initial/diversity_data.Biomass_end
plt.hist('ratio_biomass_i_end', data = diversity_data, label = "S_end")
plt.xlabel("Biomass: initial/final")

#%%
g = sns.FacetGrid(data = diversity_data, col = 'activity', row = 'DOC_initial_int', hue = "carbon_species",sharey="row", sharex=True)
g.map(sns.scatterplot, "S_initial", "DOC_removal", alpha = 0.7)
g.add_legend()
plt.savefig(os.path.join(figures_dir, "adapt_DOC_removal_S_i_compl.png"), dpi = 300)
g = sns.FacetGrid(data = diversity_data, col = 'activity', row = 'DOC_initial_int', hue = "carbon_species")
g.map(sns.scatterplot, "S_initial", "DOC_removal_mid", alpha = 0.7)
g.add_legend()
plt.savefig(os.path.join(figures_dir, "adapt_DOC_removal_mid_S_i_compl.png"), dpi = 300)
g = sns.FacetGrid(data = diversity_data, col = 'activity', row = 'DOC_initial_int', hue = "carbon_species")
g.map(sns.scatterplot, "S_initial", "ratio_t_50", alpha = 0.7)
g.add_legend()
plt.savefig(os.path.join(figures_dir, "adapt_t_50_ratio_S_i_compl.png"), dpi = 300)
g = sns.FacetGrid(data = diversity_data, col = 'activity', row = 'DOC_initial_int', hue = "carbon_species")
g.map(sns.scatterplot, "S_initial", "t_50_days", alpha = 0.7)
g.add_legend()
plt.savefig(os.path.join(figures_dir, "adapt_t_50_days_S_i_compl.png"), dpi = 300)

#%%
g = sns.FacetGrid(data = diversity_data, col = 'S_initial_int', row = 'DOC_initial_int', hue = "carbon_species",sharey="row", sharex=True)
g.map(sns.scatterplot, "activity", "ratio_t_50", alpha = 0.7)
g.add_legend()
plt.savefig(os.path.join(figures_dir, "adapt_ratio_t_50_act_compl.png"), dpi = 300)