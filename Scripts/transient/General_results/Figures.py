#%%
## Import libraries
import os
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
simulations_dir = os.path.join(project_dir, "simulations")
results_dir = os.path.join(project_dir, "results")
figures_dir = os.path.join(project_dir, "figures")

filename = os.path.join(results_dir, "combined_dataset.pkl")
diversity_data = pd.read_pickle(filename)
diversity_data = diversity_data.drop_duplicates()
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
# Filter out all failed simulations:
compl = diversity_data[diversity_data.Status>=0]
g = sns.FacetGrid(compl, col = 'activity', row = 'DOC_initial_int', hue = "carbon_species")
g.map(sns.scatterplot, "S_initial", "DOC_removal", alpha = 0.7)
g.add_legend()
plt.savefig(os.path.join(figures_dir, "DOC_removal_S_i_compl.png"), dpi = 300)
g = sns.FacetGrid(compl, col = 'activity', row = 'DOC_initial_int', hue = "carbon_species")
g.map(sns.scatterplot, "S_initial", "ratio_t_50", alpha = 0.7)
g.add_legend()
plt.savefig(os.path.join(figures_dir, "t_50_ratio_S_i_compl.png"), dpi = 300)
g = sns.FacetGrid(compl, col = 'activity', row = 'DOC_initial_int', hue = "carbon_species")
g.map(sns.scatterplot, "S_initial", "t_50_days", alpha = 0.7)
g.add_legend()
plt.savefig(os.path.join(figures_dir, "t_50_days_S_i_compl.png"), dpi = 300)

#%%
# List simulations that failed and explore parameters
failed = diversity_data[diversity_data.Status<0].reset_index()
print(failed.shape)
failed['bcratio'] = failed["biomass_species"]/failed["carbon_species"]
failed['docbcratio'] = failed.DOC_initial_int*100/(failed.carbon_species*failed.biomass_species*failed.activity)
compl['bcratio'] = compl["biomass_species"]/compl["carbon_species"]
compl['docbcratio'] = compl.DOC_initial_int*100/(compl.carbon_species*compl.biomass_species*compl.activity)
#%%
failed = diversity_data[diversity_data.Status<0].reset_index()
failed['Seed'] = failed['Seed'].map(str)
failed['carbon_species'] = failed['carbon_species'].map(str)
failed['biomass_species'] = failed['biomass_species'].map(str)
failed['DOC_initial_int'] = failed['DOC_initial_int'].astype(int).map(str)
failed['id'] = failed[['Seed', 'carbon_species', 'biomass_species', 'DOC_initial_int']].agg('-'.join, axis=1)
print(failed.id.unique())
#%%
#seeds_paras_exploration
failed['sid'] = failed['Sim_series']+"/bio_n_"+failed['biomass_species']+"/dom_initial_"+failed['DOC_initial_int']+"/seed_"+failed['Seed']
#%%
unique_c_n = list(failed['carbon_species'].unique())
for c_n in unique_c_n[:1]:
    subset = failed[failed.carbon_species==c_n]
    seed_cn_list = list(subset["Seed"].unique())
    key_sid = subset.sid.tolist()
    for seed in seed_cn_list:
        details_subfolder = 'carbon_' + c_n + '_' + seed + '_ip_0'
        simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
        filename = os.path.join(simulations_dir, "seeds_randoms.pkl")
        seed_details = pd.read_pickle(filename)
        for k in key_sid:
            print(k, seed_details[k]['message'])
