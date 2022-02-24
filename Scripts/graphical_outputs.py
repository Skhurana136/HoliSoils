#%%
## Import libraries
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

## LOAD RESULTS
project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
results_dir = os.path.join(project_dir, "results")
figures_dir = os.path.join(project_dir, "figures")

# Shannon diversity vs carbon stock
data = pd.read_pickle(os.path.join(results_dir, "diversity_data.pkl"))
data_norm = pd.read_pickle(os.path.join(results_dir, "diversity_normalized_data.pkl"))
#%%
sns.scatterplot(x="Shannon", y = "DOC", data = data)
#plt.yscale("log")

#%%
sns.scatterplot(x="Shannon", y = "DOC", data = data_norm, hue = "carbon_species", size = "biomass_species")
plt.ylabel("Normalised DOC")
#plt.yscale("log")

#%%
data_norm["S_DOC"] = data_norm["Shannon"]/data_norm["DOC"]
sns.scatterplot(x="S_DOC", y = "DOC", data = data_norm, hue = "carbon_species", size = "biomass_species", edgecolor = None)
plt.ylabel("Normalised DOC")
#plt.yscale("log")