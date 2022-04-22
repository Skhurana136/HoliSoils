#%%
## Import libraries
import os
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

## LOAD RESULTS
#project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
project_dir = "D:\Projects\HoliSoils\data"
simulations_dir = os.path.join(project_dir, "simulations")
results_dir = os.path.join(project_dir, "results")
figures_dir = os.path.join(project_dir, "figures")

filename = os.path.join(results_dir, "combined_dataset.pkl")
diversity_data = pd.read_pickle(filename)
diversity_data['DOC_initial_int'] = round(diversity_data.DOC_initial, -3)
diversity_data['ratio_t_50'] = diversity_data.T_50/diversity_data.T_50_B1
#%%
# Distribution of number of carbon compounds and biomass species that we tested.
sns.kdeplot(x = "carbon_species", data = diversity_data, label = "#carbon")
sns.kdeplot(x = "biomass_species", data = diversity_data, label = "#biomass")
plt.legend()

#%%
sns.kdeplot(x = "DOC_removal", data = diversity_data)

#%%
sns.kdeplot(x = "T_50", data = diversity_data, hue = "carbon_species")
plt.xscale("log")
plt.xlabel("Time for 37% loss (days)")
plt.tight_layout()
plt.savefig(os.path.join(figures_dir, "t_63_distribution.png"), dpi = 300)
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
sns.scatterplot(x = "carbon_species", y = "S_end", data = diversity_data, label = "S_end")
sns.scatterplot(x = "carbon_species", y = "S_initial", data = diversity_data, label = "S_i")
sns.scatterplot(x = "carbon_species", y = "S_t_50", data = diversity_data, label = "S_at_t50")
sns.scatterplot(x = "carbon_species", y = "S_t_50_b1", data = diversity_data, label = "S_at_t50_b1")
plt.ylabel("Shannon diversity index")
plt.xlabel ("number of carbon compounds")

#%%
sns.scatterplot(x = "biomass_species", y = "S_end", data = diversity_data, label = "S_end")
sns.scatterplot(x = "biomass_species", y = "S_initial", data = diversity_data, label = "S_i")
sns.scatterplot(x = "biomass_species", y = "S_t_50", data = diversity_data, label = "S_at_t50")
sns.scatterplot(x = "biomass_species", y = "S_t_50_b1", data = diversity_data, label = "S_at_t50_b1")
plt.ylabel("Shannon diversity index")
plt.xlabel ("number of biomass species")

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
g = sns.FacetGrid(diversity_data, col = 'activity', row = 'DOC_initial_int', hue = "carbon_species")
g.map(sns.scatterplot, "S_initial", "DOC_removal", alpha = 0.7)
g.add_legend()
plt.savefig(os.path.join(figures_dir, "DOC_removal_S_i.png"), dpi = 300)

#%%
g = sns.FacetGrid(diversity_data, col = 'activity', row = 'DOC_initial_int', hue = "carbon_species")
g.map(sns.scatterplot, "S_initial", "ratio_t_50", alpha = 0.7)
g.add_legend()
plt.savefig(os.path.join(figures_dir, "t_50_ratio_S_i.png"), dpi = 300)
#%%
g = sns.FacetGrid(diversity_data, col = 'activity', row = 'DOC_initial_int', hue = "carbon_species")
g.map(sns.scatterplot, "S_initial", "t_50_days", alpha = 0.7)
g.add_legend()
plt.savefig(os.path.join(figures_dir, "t_50_days_S_i.png"), dpi = 300)
#%%
# In some cases, the removal of DOC for c18 is low, or even negative.
# This was when even the baseline simulation fails. So removal all cases where T_50_B1 is NA.
all_base = diversity_data.dropna(subset = ['T_50_B1'])
g = sns.FacetGrid(all_base, col = 'activity', row = 'DOC_initial_int', hue = "carbon_species")
g.map(sns.scatterplot, "S_initial", "DOC_removal", alpha = 0.7)
g.add_legend()
plt.savefig(os.path.join(figures_dir, "DOC_removal_S_i_filtered.png"), dpi = 300)
g = sns.FacetGrid(all_base, col = 'activity', row = 'DOC_initial_int', hue = "carbon_species")
g.map(sns.scatterplot, "S_initial", "ratio_t_50", alpha = 0.7)
g.add_legend()
plt.savefig(os.path.join(figures_dir, "t_50_ratio_S_i_filtered.png"), dpi = 300)
g = sns.FacetGrid(all_base, col = 'activity', row = 'DOC_initial_int', hue = "carbon_species")
g.map(sns.scatterplot, "S_initial", "t_50_days", alpha = 0.7)
g.add_legend()
plt.savefig(os.path.join(figures_dir, "t_50_days_S_i_filtered.png"), dpi = 300)