#%%
import os
import pandas as pd
from scipy import stats

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
compl = diversity_data[diversity_data.Status>=0]
#%%
init_doc_list = list(compl.DOC_initial_int.unique())
activity_list = list(compl.activity.unique())

res = []
for doci in init_doc_list:
    for act in activity_list:
        sub = diversity_data[(diversity_data['DOC_initial_int']==doci) & (diversity_data['activity']==act)]
        X = sub['S_initial'].to_numpy(dtype = float)
        y = sub['DOC_removal'].to_numpy(dtype = float)
        b = sub["Biomass_initial"].to_numpy(dtype = float)[0]
        result = stats.linregress(X, y)
        res.append([doci, act, b, result.slope, result.intercept, result.rvalue, result.pvalue, result.stderr, result.intercept_stderr])
resdf = pd.DataFrame.from_records(res, columns = ["DOC_initial", "activity", "Biomass_initial", "slope", "intercept", "score", "pval", "slope_error", "intercept_error"])

#%%
# Plotting various possibilities of the x variable
resdf["x_variable"] = resdf["DOC_initial"]*100/(resdf["Biomass_initial"]*resdf["activity"])

fig, axes = plt.subplots(2,2, figsize = (7,6), sharex = True)#, sharey = "row")
axes[0,0].scatter(resdf["x_variable"], resdf["slope"], alpha = 0.6)
axes[0,0].set_ylabel ("Slope")
axes[0,1].scatter(resdf["x_variable"], resdf["intercept"], alpha = 0.6)
axes[0,1].set_ylabel ("Intercept")
axes[1,0].scatter(resdf["x_variable"], resdf["score"], alpha = 0.6)
axes[1,0].set_ylabel ("Score")
axes[1,1].scatter(resdf["x_variable"], resdf["slope_error"], alpha = 0.6)
axes[1,1].set_ylabel ("Slope_stderror")
#axes[1,1].set_yscale("log")
for a in axes[1,:]:
    a.set_xlabel ('Average C available\nfor active biomass (-)')
fig.tight_layout()
plt.savefig(os.path.join(figures_dir, "regression_performance.png"), dpi = 300)

#%%
res = []
for doci in init_doc_list:
    for act in activity_list:
        sub = diversity_data[(diversity_data['DOC_initial_int']==doci) & (diversity_data['activity']==act)]
        for c_n in [3,6,12,18]:
            X = sub[sub.carbon_species == c_n]['S_initial'].to_numpy(dtype = float)
            y = sub[sub.carbon_species == c_n]['DOC_removal'].to_numpy(dtype = float)
            b = sub["Biomass_initial"].to_numpy(dtype = float)[0]
            result = stats.linregress(X, y)
            res.append([doci, act, c_n, b, result.slope, result.intercept, result.rvalue, result.pvalue, result.stderr, result.intercept_stderr])

resdf = pd.DataFrame.from_records(res, columns = ["DOC_initial", "activity", "carbon_species", "Biomass_initial", "slope", "intercept", "score", "pval", "slope_error", "intercept_error"])

#%%
# Plotting various possibilities of the x variable
resdf["x_variable"] = resdf["DOC_initial"]*100/(resdf["carbon_species"]*resdf["Biomass_initial"]*resdf["activity"])

fig, axes = plt.subplots(2,2, figsize = (7,6), sharex = True)#, sharey = "row")
sns.scatterplot(resdf["x_variable"], resdf["slope"], alpha = 0.6, hue = resdf["carbon_species"], ax = axes[0,0])
axes[0,0].set_ylabel ("Slope")
sns.scatterplot(resdf["x_variable"], resdf["intercept"], alpha = 0.6, hue = resdf["carbon_species"], ax = axes[0,1])
axes[0,1].set_ylabel ("Intercept")
sns.scatterplot(resdf["x_variable"], resdf["score"], alpha = 0.6, hue = resdf["carbon_species"], ax = axes[1,0])
axes[1,0].set_ylabel ("Score")
sns.scatterplot(resdf["x_variable"], resdf["slope_error"], alpha = 0.6, hue = resdf["carbon_species"], ax = axes[1,1])
axes[1,1].set_ylabel ("Slope_stderror")
for a in axes[1,:]:
    a.set_xlabel ('Average DOC concentration per chemical\navailable for active biomass (-)')
for a in axes.flatten():
    a.set_xscale("log")
    a.set_xticks([0.4,1,2,5,10])
    a.set_xticklabels([0.4,1,2,5,10])
for a in axes.flatten():
    a.get_legend().remove()
axes[1,0].legend(title = "Carbon species", ncol = 4, bbox_to_anchor = (1.1,-0.3))
fig.tight_layout()
plt.savefig(os.path.join(figures_dir, "regression_performance_cnum_diff.png"), dpi = 300)