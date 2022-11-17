#%%
import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns
import scipy

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
all_data = pd.read_csv(os.path.join(project_dir,"simulation_results_diversity_parameters.csv"))
print(all_data.shape)
print(all_data.dtypes)
#%%
compl = all_data.dropna(subset = ['ratio_t_50'])
act_full = all_data[all_data.Sim_series=="b_1_all_"]
act_off = all_data[all_data.Sim_series!="b_1_all_"]
compl_act_off = compl[compl.activity<100]
no_na = all_data.dropna()
#%%
all_features = ['Variance','DOC_initial', 'carbon_species', 'NOSC_initial','biomass_species', 'activity', 'S_initial', 'FD_initial', 'FD_cov', 'vmax_mean', 'm_mean', 'k_mean', 'exo_enz_mean','Decay_constant', 'Norm_b1_decay_const','Norm_b2_decay_const']
all_corr = all_data[all_features].corr()
mask = np.triu(np.ones_like(all_corr, dtype=bool))
sns.heatmap(all_corr, vmin = -1, vmax = 1, center = 0, mask = mask, cmap='vlag')
#%%
features_to_test = ['Variance','DOC_initial', 'carbon_species', 'biomass_species', 'activity', 'S_initial', 'FD_initial', 'vmax_mean', 'Decay_constant', 'Norm_b1_decay_const','Norm_b2_decay_const']
feature_data = all_data[features_to_test]
print(feature_data.shape)
#%%
corr = feature_data.corr()
mask = np.triu(np.ones_like(corr, dtype=bool))
sns.heatmap(corr, vmin = -1, vmax = 1, center = 0, mask = mask, cmap='vlag')
#Limited linear relationships
#%%
subset = act_full[act_full.biomass_species <12]
g = sns.displot(data= act_full, x = "Decay_constant", hue = "DOC_initial_int", row = "carbon_species", col = "biomass_species", kind = 'kde')
g.set(xscale="log")
#%%
#%%
subset = act_full[act_full.biomass_species <8]
g = sns.jointplot(data= subset, x = "FD_initial", y = "Decay_constant", hue = "DOC_initial_int")
g.set(xscale="log")
#%%
import statsmodels.api as sm
import statsmodels.formula.api as smf
#%%
plt.figure()
sm.qqplot(np.log10(act_full['FD_initial']), line = 's')
plt.title("Initial functional diversity")
#%%
#plt.figure()
#sm.qqplot(np.log10(act_full['Decay_constant']), line = '45')
#plt.title("Decay constant")
plt.figure()
sm.qqplot(np.log10(act_full['Decay_constant']), line = 's')
plt.title("Decay constant")
#%%
#Log transformation of the dataframe
lognona = no_na.apply(lambda x: np.log10(x) if np.issubdtype(x.dtype, np.number) else x)
#%%
### Tests for normality ###
## Initial functional diversity
print("Shapiro test: ", scipy.stats.shapiro(lognona['FD_initial']))
print("KS test: ",scipy.stats.kstest(lognona['FD_initial'], scipy.stats.norm.cdf))
#%%
# Decay constant
print("Shapiro test: ", scipy.stats.shapiro(lognona['Decay_constant']))
print("KS test: ",scipy.stats.kstest(lognona['Decay_constant'], scipy.stats.norm.cdf))
#%%
### Linear Mixed Effects Models ###
results_dic = {}
#%%
#Given information is number of carbon groups, number of biomass groups,
#and level of activity%
#Given this, can we predict decay constant if we know available carbon concentration?
#In this, carbon availability has different intercepts.
cab_docgroup = smf.mixedlm('Decay_constant ~ carbon_species + biomass_species + activity', lognona, groups = lognona['DOC_initial_int'])
cab_docgroup = cab_docgroup.fit(reml=False)
print(cab_docgroup.summary())
results_dic.update({'cab_docgroup_aic': cab_docgroup.aic})
#%%
#Given information is number of carbon groups, number of biomass groups,
#and level of activity%
#Given this, can we predict decay constant if we know available carbon concentration?
#In this, carbon availability has different slope associated with the size of the carbon pool.
cab_docslopec = smf.mixedlm('Decay_constant ~ carbon_species + biomass_species + activity', lognona, groups = lognona['DOC_initial_int'], re_formula = "~ 0 + carbon_species")
cab_docslopec = cab_docslopec.fit(reml=False)
print(cab_docslopec.summary())
results_dic.update({'cab_docslopec_aic': cab_docslopec.aic})
#%%
#Given information is number of carbon groups, number of biomass groups,
#and level of activity%
#Given this, can we predict decay constant if we know available carbon concentration?
#In this, carbon availability has different slope associated with the activity%.
cab_docslopea = smf.mixedlm('Decay_constant ~ carbon_species + biomass_species + activity', lognona, groups = lognona['DOC_initial_int'], re_formula = "~ 0 + activity")
cab_docslopea = cab_docslopea.fit(reml=False)
print(cab_docslopea.summary())
results_dic.update({'cab_docslopea_aic': cab_docslopea.aic})
#%%
#Given information is number of carbon groups, number of biomass groups,
#and level of activity%
#Given this, can we predict decay constant if we know available carbon concentration?
#In this, carbon availability has different slope associated with the biomass_species.
cab_docslopeb = smf.mixedlm('Decay_constant ~ carbon_species + biomass_species + activity', lognona, groups = lognona['DOC_initial_int'], re_formula = "~ 0 + biomass_species")
cab_docslopeb = cab_docslopeb.fit(reml=False)
print(cab_docslopeb.summary())
results_dic.update({'cab_docslopeb_aic': cab_docslopeb.aic})
#%%
#Given information is number of carbon groups, number of biomass groups,
#and level of activity%
#Given this, can we predict decay constant if we know available carbon concentration?
#In this, carbon availability has different intercept and slope associated with the carbon species.
cab_docintslopeb = smf.mixedlm('Decay_constant ~ carbon_species + biomass_species + activity', lognona, groups = lognona['DOC_initial_int'], re_formula = "~biomass_species")
cab_docintslopeb = cab_docintslopeb.fit(reml=False)
print(cab_docintslopeb.summary())
results_dic.update({'cab_docintslopeb_aic': cab_docintslopeb.aic})
#%%
#Given information is number of carbon groups, number of biomass groups,
#and level of activity%
#Given this, can we predict decay constant if we know available carbon concentration?
#In this, carbon availability has different intercept and slope associated with the carbon species.
cab_docintslopec = smf.mixedlm('Decay_constant ~ carbon_species + biomass_species + activity', lognona, groups = lognona['DOC_initial_int'], re_formula = "~carbon_species")
cab_docintslopec = cab_docintslopec.fit(reml=False)
print(cab_docintslopec.summary())
results_dic.update({'cab_docintslopec_aic': cab_docintslopec.aic})
#%%
#Given information is number of carbon groups, number of biomass groups,
#and level of activity%
#Given this, can we predict decay constant if we know available carbon concentration?
#In this, carbon availability has different intercept and slope associated with the activity.
cab_docintslopea = smf.mixedlm('Decay_constant ~ carbon_species + biomass_species + activity', lognona, groups = lognona['DOC_initial_int'], re_formula = "~activity")
cab_docintslopea = cab_docintslopea.fit(reml=False)
print(cab_docintslopea.summary())
results_dic.update({'cab_docintslopea_aic': cab_docintslopea.aic})
#%%
#Given information is number of carbon groups, number of biomass groups,
#and level of activity%
#Given this, can we predict decay constant if we know available carbon concentration?
#In this, carbon availability has different intercept and slope associated with all the fixed factors together.
cab_docintslopeabc = smf.mixedlm('Decay_constant ~ carbon_species + biomass_species + activity', lognona, groups = lognona['DOC_initial_int'], re_formula = "~carbon_species+biomass_species+activity")
cab_docintslopeabc = cab_docintslopeabc.fit(reml=False)
print(cab_docintslopeabc.summary())
results_dic.update({'cab_docintslopeabc_aic': cab_docintslopeabc.aic})
#Don't take this forward
#%%
#Given information is number of carbon groups, number of biomass groups,
#and level of activity%
#Given this, can we predict decay constant if we know available carbon concentration?
#In this, carbon availability has different slope associated with all the fixed factors together.
cab_docslopeabc = smf.mixedlm('Decay_constant ~ carbon_species + biomass_species + activity', lognona, groups = lognona['DOC_initial_int'], re_formula = "~0+carbon_species+biomass_species+activity")
cab_docslopeabc = cab_docslopeabc.fit(reml=False)
print(cab_docslopeabc.summary())
results_dic.update({'cab_docslopeabc_aic': cab_docslopeabc.aic})
#%%
### Additionally, if we also know the community structure, would this be different
### as groups or fixed effects?
print("Variance as fixed effect")
cabv_docslopea = smf.mixedlm('Decay_constant ~ carbon_species + biomass_species + activity + Variance', lognona, groups = lognona['DOC_initial_int'], re_formula = "~0+activity")
cabv_docslopea = cabv_docslopea.fit(reml=False)
print(cabv_docslopea.summary())
results_dic.update({'cabv_docslopea_aic': cabv_docslopea.aic})
#%%
print("Variance as fixed effect")
cabv_docslopev = smf.mixedlm('Decay_constant ~ carbon_species + biomass_species + activity + Variance', lognona, groups = lognona['DOC_initial_int'], re_formula = "~0+Variance")
cabv_docslopev = cabv_docslopev.fit(reml=False)
print(cabv_docslopev.summary())
results_dic.update({'cabv_docslopev_aic': cabv_docslopev.aic})
#%%
cabv_docintslopev = smf.mixedlm('Decay_constant ~ carbon_species + biomass_species + activity + Variance', lognona, groups = lognona['DOC_initial_int'], re_formula = "~Variance")
cabv_docintslopev = cabv_docintslopev.fit(reml=False)
print(cabv_docintslopev.summary())
results_dic.update({'cabv_docintslopev_aic': cabv_docintslopev.aic})
#%%
cab_docintslopev = smf.mixedlm('Decay_constant ~ carbon_species + biomass_species + activity', lognona, groups = lognona['DOC_initial_int'], re_formula = "~Variance")
cab_docintslopev = cab_docintslopev.fit(reml=False)
print(cab_docintslopev.summary())
results_dic.update({'cab_docintslopev_aic': cab_docintslopev.aic})
#%%
cab_docslopev = smf.mixedlm('Decay_constant ~ carbon_species + biomass_species + activity', lognona, groups = lognona['DOC_initial_int'], re_formula = "~0+Variance")
cab_docslopev = cab_docslopev.fit(reml=False)
print(cab_docslopev.summary())
results_dic.update({'cab_docslopev_aic': cab_docslopev.aic})
#%%
print("Variance as grouped effect")
cab_vgroup = smf.mixedlm('Decay_constant ~ carbon_species + biomass_species + activity', lognona, groups = lognona['Variance'])
cab_vgroup = cab_vgroup.fit(reml=False)
print(cab_vgroup.summary())
results_dic.update({'cab_vgroup_aic': cab_vgroup.aic})
#%%
fd_docgroup = smf.mixedlm('Decay_constant ~ FD_initial', lognona, groups = lognona['DOC_initial_int'])
fd_docgroup = fd_docgroup.fit(reml=False)
print(fd_docgroup.summary())
results_dic.update({'fd_docgroup_aic': fd_docgroup.aic})
#Take this forward
#%%
fd_docintslopfd = smf.mixedlm('Decay_constant ~ FD_initial', lognona, groups = lognona['DOC_initial_int'], re_formula = "~FD_initial")
fd_docintslopfd = fd_docintslopfd.fit(reml=False)
print(fd_docintslopfd.summary())
results_dic.update({'fd_docintslopfd_aic': fd_docintslopfd.aic})
#Don't take this forward. Not converged.
#%%
fd_docslopefd = smf.mixedlm('Decay_constant ~ FD_initial', lognona, groups = lognona['DOC_initial_int'], re_formula = "~0+FD_initial")
fd_docslopefd = fd_docslopefd.fit(reml=False)
print(fd_docslopefd.summary())
results_dic.update({'fd_docslopefd_aic': fd_docslopefd.aic})

### ORDINARY LINEAR REGRESSION ###
#%%
fd = smf.ols('Decay_constant ~ FD_initial', lognona)
fd = fd.fit()
print(fd.summary())
results_dic.update({'fd_aic': fd.aic})
#%%
fd_doc = smf.ols('Decay_constant ~ FD_initial + DOC_initial_int', lognona)
fd_doc = fd_doc.fit()
print(fd_doc.summary())
results_dic.update({'fd_doc_aic': fd_doc.aic})
#%%
doc = smf.ols('Decay_constant ~ DOC_initial_int', lognona)
doc = doc.fit()
print(doc.summary())
results_dic.update({'doc_aic': doc.aic})
#%%
fd_docfd = smf.ols('Decay_constant ~ FD_initial + DOC_initial_int:FD_initial', lognona)
fd_docfd = fd_docfd.fit()
print(fd_docfd.summary())
results_dic.update({'fd_docfd_aic': fd_docfd.aic})
#%%
fd_doc_docfd = smf.ols('Decay_constant ~ FD_initial + DOC_initial_int + DOC_initial_int:FD_initial', lognona)
fd_doc_docfd = fd_doc_docfd.fit()
print(fd_doc_docfd.summary())
results_dic.update({'fd_doc_docfd_aic': fd_doc_docfd.aic})
#%%
#%%
cabvdoc = smf.ols('Decay_constant ~ carbon_species + biomass_species + activity + Variance + DOC_initial_int', lognona)
cabvdoc = cabvdoc.fit()
print(cabvdoc.summary())
results_dic.update({'cabvdoc_aic': cabvdoc.aic})
#%%
cabv_doca = smf.ols('Decay_constant ~ carbon_species + biomass_species + activity + Variance + DOC_initial_int:activity', lognona)
cabv_doca = cabv_doca.fit()
print(cabv_doca.summary())
results_dic.update({'cabv_doca_aic': cabv_doca.aic})
#%%
cabvdoc_doca = smf.ols('Decay_constant ~ carbon_species + biomass_species + activity + Variance + DOC_initial_int + DOC_initial_int:activity', lognona)
cabvdoc_doca = cabvdoc_doca.fit()
print(cabvdoc_doca.summary())
results_dic.update({'cabvdoc_doca_aic': cabvdoc_doca.aic})
#%%
cabdoc_docv = smf.ols('Decay_constant ~ carbon_species + biomass_species + activity + DOC_initial_int + DOC_initial_int:Variance', lognona)
cabdoc_docv = cabdoc_docv.fit()
print(cabdoc_docv.summary())
results_dic.update({'cabdoc_docv_aic': cabdoc_docv.aic})
#%%
cab_docv = smf.ols('Decay_constant ~ carbon_species + biomass_species + activity + DOC_initial_int:Variance', lognona)
cab_docv = cab_docv.fit()
print(cab_docv.summary())
results_dic.update({'cab_docv_aic': cab_docv.aic})
#%%
cab_doca = smf.ols('Decay_constant ~ carbon_species + biomass_species + activity + DOC_initial_int:activity', lognona)
cab_doca = cab_doca.fit()
print(cab_doca.summary())
results_dic.update({'cab_doca_aic': cab_doca.aic})
#%%
cabdoc = smf.ols('Decay_constant ~ carbon_species + biomass_species + activity + DOC_initial_int', lognona)
cabdoc = cabdoc.fit()
print(cabdoc.summary())
results_dic.update({'cabdoc_aic': cabdoc.aic})
#%%
cabdoc_doccab = smf.ols('Decay_constant ~ carbon_species + biomass_species + activity + DOC_initial_int + DOC_initial_int:(activity+carbon_species+biomass_species)', lognona)
cabdoc_doccab = cabdoc_doccab.fit()
print(cabdoc_doccab.summary())
results_dic.update({'cabdoc_doccab_aic': cabdoc_doccab.aic})
#%%
### Does Carbon amplify the effect of Variance? Are they interacting?
cab_docgroup_vinter = smf.mixedlm('Decay_constant ~ carbon_species + biomass_species + activity', lognona, groups = lognona['DOC_initial_int'], re_formula = "~Variance")
cab_docgroup_vinter = cab_docgroup_vinter.fit(reml=False)
print(cab_docgroup_vinter.summary())
#%%
### If I only take 100% activity data:
act_full_log = act_full.apply(lambda x: np.log10(x) if np.issubdtype(x.dtype, np.number) else x).dropna(subset=['Decay_constant'])
#%%
full_cab_docgroup = smf.mixedlm('Decay_constant ~ carbon_species + biomass_species', act_full_log, groups = act_full_log['DOC_initial_int'])
full_cab_docgroup = full_cab_docgroup.fit(reml=False)
print(full_cab_docgroup.summary())
#%%
full_cab_vgroup = smf.mixedlm('Decay_constant ~ carbon_species + biomass_species', act_full_log, groups = act_full_log['Variance'])
full_cab_vgroup = full_cab_vgroup.fit(reml=False)
print(full_cab_vgroup.summary())
#%%
full_cabv_docgroup = smf.mixedlm('Decay_constant ~ carbon_species + biomass_species + Variance', act_full_log, groups = act_full_log['DOC_initial_int'])
full_cabv_docgroup = full_cabv_docgroup.fit(reml=False)
print(full_cabv_docgroup.summary())
#%%
full_cabvdoc = smf.ols('Decay_constant ~ carbon_species + biomass_species + Variance + DOC_initial_int', act_full_log)
full_cabvdoc = full_cabvdoc.fit()
print(full_cabvdoc.summary())
#%%
full_cabdoc_vgroup = smf.mixedlm('Decay_constant ~ carbon_species + biomass_species + DOC_initial_int', act_full_log, groups = act_full_log["Variance"])
full_cabdoc_vgroup = full_cabdoc_vgroup.fit()
print(full_cabdoc_vgroup.summary())

#%%
mdf_4 = smf.ols('Decay_constant ~ FD_initial', act_full_log)
mdf_4 = mdf_4.fit()
print(mdf_4.summary())
#%%
mdf_5 = smf.ols('Decay_constant ~ FD_initial + DOC_initial_int', act_full_log)
mdf_5 = mdf_5.fit()
print(mdf_5.summary())
#%%
mdf_6 = smf.ols('Decay_constant ~ DOC_initial_int', act_full_log)
mdf_6 = mdf_6.fit()
print(mdf_6.summary())
#%%
mdf_7 = smf.ols('Decay_constant ~ FD_initial + DOC_initial_int:FD_initial', act_full_log)
mdf_7 = mdf_7.fit()
print(mdf_7.summary())
#%%
#%%
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.linear_model import Ridge
from sklearn.inspection import permutation_importance
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_pinball_loss, mean_squared_error
from sklearn.svm import SVR
pipe_ridge = make_pipeline(StandardScaler(), Ridge())
common_params = dict(
    learning_rate=0.05,
    n_estimators=200,
    max_depth=2,
    min_samples_leaf=9,
    min_samples_split=9,
)
gbr = GradientBoostingRegressor(loss="quantile", alpha=0.5, **common_params)
pipe_gbrt = make_pipeline(StandardScaler(), gbr)
def run_ml (X_set, y_set, features_to_test, pipe_to_use):
    X = np.array(X_set)
    y = np.array(y_set)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    pipe_to_use.fit(X_train, y_train)  # apply scaling on training data
    y_pred = pipe_to_use.predict(X_test)
    print("Model performances")
    print("Score: ", f"{pipe_to_use.score(X_test, y_test):.3f}")
    print("MSE: ", f"{mean_squared_error(y_test,y_pred):.5f}")
    print("MSE%: ", f"{mean_squared_error(y_test,y_pred)/np.mean(y_test):.1f}")
    r = permutation_importance(pipe_to_use, X_test, y_test,
                                n_repeats=30,
                                random_state=0)
    print("Feature importances")
    for i in r.importances_mean.argsort()[::-1]:
        if r.importances_mean[i] - 2 * r.importances_std[i] > 0:
            print(f"{features_to_test[i]:<8}"
                f" {r.importances_mean[i]:.3f}"
                f" +/- {r.importances_std[i]:.3f}")
    return None
gbrt_features = ['DOC_initial', 'carbon_species', 'biomass_species', 'FD_initial', 'vmax_mean']
gbrt_features_2 = ['DOC_initial', 'carbon_species', 'biomass_species', 'activity', 'Variance']
gbrt_features_3 = ['DOC_initial', 'carbon_species', 'biomass_species', 'activity', 'Variance', 'vmax_mean']
gbrt_features_4 = ['DOC_initial', 'activity', 'FD_initial', 'vmax_mean', 'S_initial']
#%%
X = np.array(no_na[gbrt_features])
y = np.array(no_na[['Decay_constant']])
print("GBRT: ")
run_ml(X,y,gbrt_features, pipe_gbrt)
print("Ridge: ")
run_ml(X,y,gbrt_features, pipe_ridge)
#%%
X = np.array(no_na[gbrt_features_2])
y = np.array(no_na[['Decay_constant']])
print("GBRT: ")
run_ml(X,y,gbrt_features_2, pipe_gbrt)
print("Ridge: ")
run_ml(X,y,gbrt_features_2, pipe_ridge)
#%%
#%%
X = np.array(no_na[gbrt_features_3])
y = np.array(no_na[['Decay_constant']])
print("GBRT: ")
run_ml(X,y,gbrt_features_3, pipe_gbrt)
print("Ridge: ")
run_ml(X,y,gbrt_features_3, pipe_ridge)
#%%
#%%
X = np.array(no_na[gbrt_features_4])
y = np.array(no_na[['Decay_constant']])
print("GBRT: ")
run_ml(X,y,gbrt_features_4, pipe_gbrt)
print("Ridge: ")
run_ml(X,y,gbrt_features_4, pipe_ridge)
print("SVR: ")
run_ml(X,y,gbrt_features_4, pipe_svr)
plot_true_vs_pred(y,pipe_svr.predict(X))
#%%
X = np.array(no_na[gbrt_features_4])
y = np.array(no_na[['Decay_constant']])
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
pipe_gbrt.fit(X_train, y_train)  # apply scaling on training data
y_pred = pipe_gbrt.predict(X_test)
print("Model performances")
print("Score: ", f"{pipe_gbrt.score(X_test, y_test):.3f}")
print("MSE: ", f"{mean_squared_error(y_test,y_pred):.5f}")
print("MSE%: ", f"{mean_squared_error(y_test,y_pred)/np.mean(y_test):.1f}")
#%%
def plot_true_vs_pred(y, y_pred_all):
    labelsize = 16
    titles = 18
    ticksize = 14
    plt.figure(figsize = (6,6))
    plt.scatter(y, y_pred_all, alpha = 0.5, label = "Data")
    plt.plot(y,y, 'r-', label = "1:1 line")
    plt.xticks(fontsize = ticksize)
    plt.yticks(fontsize = ticksize)
    plt.xlabel("True values", fontsize = ticksize)
    plt.ylabel ("Predictions", fontsize = labelsize)
    plt.legend(title_fontsize = titles, fontsize = labelsize)

    return None
#%%
plot_true_vs_pred(y,pipe_gbrt.predict(X))
plot_true_vs_pred(y,pipe_ridge.predict(X))
#%%
svr = SVR(kernel = 'linear')
pipe_svr = make_pipeline(StandardScaler(), svr)
X = np.array(no_na[gbrt_features_4])
y = np.array(no_na[['Decay_constant']])
print("SVR: ")
run_ml(X,y,gbrt_features_4, pipe_svr)
plot_true_vs_pred(y,pipe_svr.predict(X))
#%%
fig, axes = plt.subplots(2,4, sharey = True, figsize=(14,6))
axflat = axes.flatten()
axflat[0].scatter(no_na['carbon_species'],no_na['Decay_constant'], alpha = 0.7, )
axflat[0].set_xlabel("Carbon species")
axflat[0].set_ylabel("Decay constant")

axflat[1].scatter(no_na['biomass_species'],no_na['Decay_constant'], alpha = 0.7)
axflat[1].set_xlabel("Biomass species")

axflat[2].scatter(no_na['DOC_initial_int'],no_na['Decay_constant'], alpha = 0.7)
axflat[2].set_xlabel("Carbon availability")

axflat[3].scatter(no_na['activity'],no_na['Decay_constant'], alpha = 0.7)
axflat[3].set_xlabel("Activity (%)")

axflat[4].scatter(no_na['Variance'],no_na['Decay_constant'], alpha = 0.7)
axflat[4].set_xlabel("Variance")
axflat[4].set_ylabel("Decay constant")

axflat[5].scatter(no_na['FD_initial'],no_na['Decay_constant'], alpha = 0.7)
axflat[5].set_xscale("log")
axflat[5].set_xlabel("Initial functional diversity: Variance")

axflat[6].scatter(no_na['FD_cov'],no_na['Decay_constant'], alpha = 0.7)
axflat[6].set_xscale("log")
axflat[6].set_xlabel("Initial functional diversity: CoV")

#%%


