#%%
import os
import pandas as pd
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#%%
## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
simulations_dir = os.path.join(project_dir, "simulations")
results_dir = os.path.join(project_dir, "results")
figures_dir = os.path.join(project_dir, "figures")

filename = os.path.join(results_dir, "null_cont_combined_dataset.pkl")
diversity_data = pd.read_pickle(filename)
diversity_data = diversity_data.drop_duplicates()
diversity_data['DOC_initial_int'] = round(diversity_data.DOC_initial, -3)
diversity_data['ratio_t_50'] = diversity_data.T_50/diversity_data.T_50_B1
diversity_data['S_initial_int'] = round(diversity_data.S_initial, 1)
compl = diversity_data[diversity_data.Status>=0]
init_doc_list = list(compl.DOC_initial_int.unique())
activity_list = list(compl.activity.unique())

#%%
res = []
for doci in init_doc_list:
    for act in activity_list:
        sub = compl[(compl['DOC_initial_int']==doci) & (compl['activity']==act)].sort_values(by=['S_initial'])
        #csub = sub
        for c_n in [3,6,12,18]:
            csub = sub[sub.carbon_species == c_n]
            mean_df = csub.groupby("S_initial_int", as_index=False)["ratio_t_50"].median()
            bio_df = csub.groupby("S_initial_int", as_index=False)["Biomass_initial"].mean()
            y_mean = mean_df['ratio_t_50'].to_numpy(dtype=float)
            X = mean_df['S_initial_int'].to_numpy(dtype = float)
            b = bio_df["Biomass_initial"].to_numpy(dtype = float)[0]
            if y_mean.size>0:
                slope, intercept, rval, pval, std_err = stats.linregress(X,y_mean)
                y_pred = slope*X + intercept
                rval = 1 - (np.sum((y_mean - y_pred)**2)/np.sum((y_mean - np.mean(y_mean))**2))
                std_err = csub.ratio_t_50
            else:
                slope, intercept, rval, pval, std_err = [np.nan, np.nan, np.nan, np.nan]
            res.append([doci, act, c_n, b, slope, intercept, rval, pval, std_err])

resdf = pd.DataFrame.from_records(res, columns = ["DOC_initial", "activity", "carbon_species", "Biomass_initial", "slope", "intercept", "rval", "pval", "std_err"])
# Check how the plot looks like against the data points for specific cases
for doci in init_doc_list:
    for act in activity_list:
        sub = compl[(compl['DOC_initial_int']==doci) & (compl['activity']==act)].sort_values(by=['S_initial'])
        #csub = sub
        for c_n in [3,6,12,18]:
            csub = sub[sub.carbon_species == c_n]
            mean_df = csub.groupby("S_initial_int", as_index=False)["ratio_t_50"].median()
            bio_df = csub.groupby("S_initial_int", as_index=False)["Biomass_initial"].mean()
            y_mean = mean_df['ratio_t_50'].to_numpy(dtype=float)
            X = mean_df['S_initial_int'].to_numpy(dtype = float)
            b = bio_df["Biomass_initial"].to_numpy(dtype = float)
            ressub = resdf[(resdf['DOC_initial']==doci) & (resdf['activity']==act)].reset_index()
            ypred = ressub.slope[0]*X +ressub.intercept[0]
            plt.figure()
            plt.title(str(doci) + " " + str(act) + " " + str(c_n))
            plt.plot(X, ypred, 'r-', alpha = 0.6, label = 'linres')
            plt.scatter(sub.S_initial,sub.ratio_t_50, marker = '.', label = "all_data")
            plt.scatter(X,y_mean, marker = '^', c = 'purple', label = "res_data")
            plt.legend()

#%%
resdf["y_variable"] = resdf["DOC_initial"]
resdf["x_variable"] = resdf["Biomass_initial"]*resdf["activity"]/100

#%%

fig, axes = plt.subplots(2,2, figsize = (7,6), sharex = True, sharey = True)
sns.scatterplot(resdf["x_variable"], resdf["y_variable"], size = resdf["slope"], alpha = 0.5, ax= axes[0,0])
#axes[0,0].set_ylabel ("Slope")
sns.scatterplot(resdf["x_variable"], resdf["y_variable"], hue = resdf["intercept"], alpha = 0.6, ax= axes[0,1])
#axes[0,1].set_ylabel ("Intercept")
sns.scatterplot(resdf["x_variable"], resdf["y_variable"], hue = resdf["rval"], alpha = 0.6, ax= axes[1,0])
#axes[1,0].set_ylabel ("Score")
sns.scatterplot(resdf["x_variable"], resdf["y_variable"], hue = resdf["std_err"], alpha = 0.6, ax= axes[1,1])
#axes[1,1].set_ylabel ("Slope_stderror")
#axes[1,1].set_yscale("log")
#for a in axes[1,:]:
#    a.set_xlabel ('Average C available\nfor active biomass (-)')
fig.tight_layout()
#plt.savefig(os.path.join(figures_dir, "adaptation_regression_performance_char_rxn_tim_ratio.png"), dpi = 300)

#%%
#%%
res = []
for doci in init_doc_list:
    for act in activity_list:
        sub = compl[(compl['DOC_initial_int']==doci) & (compl['activity']==act)].sort_values(by=['S_initial'])
        for c_n in [3,6,12,18]:
            csub = sub[sub.carbon_species == c_n]
            y = np.log(csub['ratio_t_50'].to_numpy(dtype=float))
            X = csub['S_initial_int'].to_numpy(dtype = float)
            b = csub["Biomass_initial"].to_numpy(dtype = float)[0]
            if y.size>0:
                slope, intercept, rval, pval, std_err = stats.linregress(X,y)
            else:
                slope, intercept, rval, pval, std_err = [np.nan, np.nan, np.nan, np.nan]
            res.append([doci, act, c_n, b, slope, intercept, rval, pval, std_err])

resdf = pd.DataFrame.from_records(res, columns = ["DOC_initial", "activity", "carbon_species", "Biomass_initial", "slope", "intercept", "rval", "pval", "std_err"])
# Check how the plot looks like against the data points for specific cases
for doci in init_doc_list:
    for act in activity_list:
        sub = compl[(compl['DOC_initial_int']==doci) & (compl['activity']==act)].sort_values(by=['S_initial'])
        #csub = sub
        for c_n in [3,6,12,18]:
            csub = sub[sub.carbon_species == c_n]
            y = np.log(csub['ratio_t_50'].to_numpy(dtype=float))
            if y.size>0:
                X = csub['S_initial_int'].to_numpy(dtype = float)
                b = csub["Biomass_initial"].to_numpy(dtype = float)[0]
                ressub = resdf[(resdf['DOC_initial']==doci) & (resdf['activity']==act)].reset_index()
                ypred = np.e**(ressub.slope[0]*X +ressub.intercept[0])
                plt.figure()
                plt.title(str(doci) + " " + str(act) + " " + str(c_n))
                plt.plot(X, ypred, 'r-', alpha = 0.6, label = 'linres')
                plt.scatter(X,np.e**(y), marker = '^', c = 'purple', label = "res_data")
                plt.legend()
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
plt.savefig(os.path.join(figures_dir, "adaptation_regression_performance_char_rxn_tim_ratio.png"), dpi = 300)

#%%
res = []
for doci in init_doc_list:
    for act in activity_list:
        sub = diversity_data[(diversity_data['DOC_initial_int']==doci) & (diversity_data['activity']==act)]
        for c_n in [3,6,12,18]:
            X = sub[sub.carbon_species == c_n]['S_initial'].to_numpy(dtype = float)
            y = sub[sub.carbon_species == c_n]['ratio_t_50'].to_numpy(dtype = float)
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
plt.savefig(os.path.join(figures_dir, "adaptation_regression_performance_cnum_diff_char_rxn_tim_ratio.png"), dpi = 300)

#%%
# Plotting various possibilities of the x variable
resdf["y_variable"] = resdf["DOC_initial"]/resdf["carbon_species"]
resdf["x_variable"] = resdf["Biomass_initial"]*resdf["activity"]/100

fig, axes = plt.subplots(2,2, figsize = (7,6), sharex = True, sharey = True)
sns.scatterplot(resdf["x_variable"], resdf["y_variable"], alpha = 0.6, hue = resdf["slope"], ax = axes[0,0])
#axes[0,0].set_ylabel ("Slope")
sns.scatterplot(resdf["x_variable"], resdf["y_variable"], alpha = 0.6, hue= resdf["intercept"],ax = axes[0,1])
#axes[0,1].set_ylabel ("Intercept")
sns.scatterplot(resdf["x_variable"], resdf["y_variable"], alpha = 0.6,  hue= resdf["score"],size = 1/resdf["pval"],ax = axes[1,0])
#axes[1,0].set_ylabel ("Score")
sns.scatterplot(resdf["x_variable"], resdf["y_variable"], alpha = 0.6, hue= resdf["slope_error"], ax = axes[1,1])
#axes[1,1].set_ylabel ("Slope_stderror")
#for a in axes[1,:]:
#    a.set_xlabel ('Average DOC concentration per chemical\navailable for active biomass (-)')
#for a in axes.flatten():
#    a.set_xscale("log")
#    a.set_xticks([0.4,1,2,5,10])
#    a.set_xticklabels([0.4,1,2,5,10])
#for a in axes.flatten():
#    a.get_legend().remove()
#axes[1,0].legend(title = "Carbon species", ncol = 4, bbox_to_anchor = (1.1,-0.3))
fig.tight_layout()
#plt.savefig(os.path.join(figures_dir, "adaptation_regression_performance_cnum_diff_char_rxn_tim_ratio.png"), dpi = 300)