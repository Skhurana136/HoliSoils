#%%
import os
import pandas as pd

import scipy.stats as stats

import seaborn as sns
import matplotlib.pyplot as plt

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
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
null_data['S_initial_rnd'] = round(null_data.S_initial, 1)
null_data['ratio_t_50'] = null_data.T_50/null_data.T_50_B1
adapt_data = adapt_data.drop_duplicates()
adapt_data['DOC_initial_int'] = round(adapt_data.DOC_initial, -3)
adapt_data['S_initial_rnd'] = round(adapt_data.S_initial, 1)
adapt_data['ratio_t_50'] = adapt_data.T_50/adapt_data.T_50_B1
#%%
init_doc_list = list(null_data.DOC_initial_int.unique())
activity_list = list(null_data.activity.unique())
s_list = list(null_data.S_initial_rnd.unique())

res=[]
for doci in init_doc_list:
    for act in activity_list:
        for s in s_list:
            null = null_data[(null_data['DOC_initial_int']==doci) & (null_data['activity']==act) & (null_data['S_initial_rnd']==s)].dropna()
            adapt = adapt_data[(adapt_data['DOC_initial_int']==doci) & (adapt_data['activity']==act) & (adapt_data['S_initial_rnd']==s)].dropna()
            
            if (null.shape[0]>0) and (adapt.shape[0]>0):
                fvalue, pvalue = stats.f_oneway(null.ratio_t_50,adapt.ratio_t_50)
                res.append([doci, act, s, "ratio_t_50", fvalue, pvalue])

                fvalue, pvalue = stats.f_oneway(null.DOC_removal,adapt.DOC_removal)
                res.append([doci, act, s, "doc_removal", fvalue, pvalue])

                fvalue, pvalue = stats.f_oneway(null.t_50_days,adapt.t_50_days)
                res.append([doci, act, s, "t_50", fvalue, pvalue])

                fvalue, pvalue = stats.f_oneway(null.DOC_removal_mid,adapt.DOC_removal_mid)
                res.append([doci, act, s, "doc_removal_mid", fvalue, pvalue])
            else:
                pass


resdf = pd.DataFrame.from_records(res, columns = ["DOC_initial", "activity", "S_initial", "Metric", "fval", "pval"])
resdf.to_csv(os.path.join(results_dir, "anova_null_adapt.csv"))


#%%
resdf.loc[resdf["pval"]<0.1, "significance"] = 1
resdf.loc[resdf["pval"]>=0.1, "significance"] = 0
metric_to_plot = list(resdf.Metric.unique())
title_list = ["Impact on time", "Carbon removed", "Time scale", "Carbon removed at midpoint"]
resdf["x_var"] = resdf.S_initial*resdf.activity/100

fig, axes = plt.subplots(2,2, sharex = True, sharey = True)
for y, t in zip(metric_to_plot, title_list):
    nullsubset = resdf[(resdf.Metric==y)&(resdf.significance==0)]
    subset = resdf[(resdf.Metric==y)&(resdf.significance==1)]
    axes.flatten()[metric_to_plot.index(y)].scatter("x_var", "DOC_initial", alpha = 0.5, data = nullsubset, c= "grey")
    axes.flatten()[metric_to_plot.index(y)].scatter("x_var", "DOC_initial", alpha = 0.7, data = subset)
    axes.flatten()[metric_to_plot.index(y)].set_title(t)


for a in axes[-1,:]:
    a.set_xlabel("% active initial diversity")
for a in axes[:,0]:
    a.set_ylabel ("Initial DOC (uM)")

plt.savefig(os.path.join(figures_dir, "anova_test.png"), dpi = 300)

#%%
init_doc_list = list(null_data.DOC_initial_int.unique())
activity_list = list(null_data.activity.unique())
s_list = list(null_data.S_initial_rnd.unique())

res=[]
for doci in init_doc_list:
    for act in activity_list:
        for s in s_list:
            null = null_data[(null_data['DOC_initial_int']==doci) & (null_data['activity']==act) & (null_data['S_initial_rnd']==s)]#.dropna()
            adapt = adapt_data[(adapt_data['DOC_initial_int']==doci) & (adapt_data['activity']==act) & (adapt_data['S_initial_rnd']==s)]#.dropna()
            print(doci, act, s)
            if (null.shape[0]>0) and (adapt.shape[0]>0):
                fvalue, pvalue = stats.mannwhitneyu(null.DOC_removal,adapt.DOC_removal, alternative = "greater")#, method = "exact")#, nan_policy = "omit")
                res.append([doci, act, s, "doc_removal", fvalue, pvalue])

                fvalue, pvalue = stats.mannwhitneyu(null.t_50_days,adapt.t_50_days, alternative = "greater")#, method = "exact")#, nan_policy = "omit")
                res.append([doci, act, s, "t_50", fvalue, pvalue])

                fvalue, pvalue = stats.mannwhitneyu(null.DOC_removal_mid,adapt.DOC_removal_mid, alternative = "greater")#, method = "exact")#, nan_policy = "omit")
                res.append([doci, act, s, "doc_removal_mid", fvalue, pvalue])
            else:
                pass


resdf = pd.DataFrame.from_records(res, columns = ["DOC_initial", "activity", "S_initial", "Metric", "fval", "pval"])
resdf.to_csv(os.path.join(results_dir, "mannwhitney_null_adapt.csv"))
#%%
resdf.loc[resdf["pval"]<0.1, "significance"] = 1
resdf.loc[resdf["pval"]>=0.1, "significance"] = 0
metric_to_plot = list(resdf.Metric.unique())
title_list = ["Carbon removed", "Time scale", "Carbon removed at midpoint"]
resdf["x_var"] = resdf.S_initial*resdf.activity/100

fig, axes = plt.subplots(3,1, figsize = (3,7), sharex = True, sharey = True)
for y, t in zip(metric_to_plot, title_list):
    nullsubset = resdf[(resdf.Metric==y)&(resdf.significance==0)]
    subset = resdf[(resdf.Metric==y)&(resdf.significance==1)]
    axes.flatten()[metric_to_plot.index(y)].scatter("x_var", "DOC_initial", alpha = 0.5, data = nullsubset, c= "grey")
    axes.flatten()[metric_to_plot.index(y)].scatter("x_var", "DOC_initial", alpha = 0.7, data = subset)
    axes.flatten()[metric_to_plot.index(y)].set_title(t)

axes[-1].set_xlabel("% active initial diversity")
for a in axes[:]:
    a.set_ylabel ("Initial DOC (uM)")

plt.savefig(os.path.join(figures_dir, "mannwhitney_test.png"), dpi = 300)

#%%
init_doc_list = list(null_data.DOC_initial_int.unique())
activity_list = list(null_data.activity.unique())
s_list = list(null_data.S_initial_rnd.unique())

res=[]
for doci in init_doc_list:
    for act in activity_list:
        for s in s_list:
            null = null_data[(null_data['DOC_initial_int']==doci) & (null_data['activity']==act) & (null_data['S_initial_rnd']==s)]#.dropna()
            adapt = adapt_data[(adapt_data['DOC_initial_int']==doci) & (adapt_data['activity']==act) & (adapt_data['S_initial_rnd']==s)]#.dropna()
            print(doci, act, s)
            if (null.shape[0]>0) and (adapt.shape[0]>0):
                fvalue, pvalue = stats.ttest_ind(null.ratio_t_50,adapt.ratio_t_50, alternative = "greater")#, method = "exact")#nan_policy = "omit")
                res.append([doci, act, s, "ratio_t_50", fvalue, pvalue])

                fvalue, pvalue = stats.ttest_ind(null.DOC_removal,adapt.DOC_removal, alternative = "greater")#, method = "exact")#, nan_policy = "omit")
                res.append([doci, act, s, "doc_removal", fvalue, pvalue])

                fvalue, pvalue = stats.ttest_ind(null.t_50_days,adapt.t_50_days, alternative = "greater")#, method = "exact")#, nan_policy = "omit")
                res.append([doci, act, s, "t_50", fvalue, pvalue])

                fvalue, pvalue = stats.ttest_ind(null.DOC_removal_mid,adapt.DOC_removal_mid, alternative = "greater")#, method = "exact")#, nan_policy = "omit")
                res.append([doci, act, s, "doc_removal_mid", fvalue, pvalue])
            else:
                pass


resdf = pd.DataFrame.from_records(res, columns = ["DOC_initial", "activity", "S_initial", "Metric", "fval", "pval"])
resdf.to_csv(os.path.join(results_dir, "ttest_null_adapt.csv"))

#%%
resdf.loc[resdf["pval"]<0.1, "significance"] = 1
resdf.loc[resdf["pval"]>=0.1, "significance"] = 0
metric_to_plot = list(resdf.Metric.unique())
title_list = ["Impact on time", "Carbon removed", "Time scale", "Carbon removed at midpoint"]
resdf["x_var"] = resdf.S_initial*resdf.activity/100

fig, axes = plt.subplots(2,2, sharex = True, sharey = True)
for y, t in zip(metric_to_plot, title_list):
    nullsubset = resdf[(resdf.Metric==y)&(resdf.significance==0)]
    subset = resdf[(resdf.Metric==y)&(resdf.significance==1)]
    axes.flatten()[metric_to_plot.index(y)].scatter("x_var", "DOC_initial", alpha = 0.5, data = nullsubset, c= "grey")
    axes.flatten()[metric_to_plot.index(y)].scatter("x_var", "DOC_initial", alpha = 0.7, data = subset)
    axes.flatten()[metric_to_plot.index(y)].set_title(t)


for a in axes[-1,:]:
    a.set_xlabel("% active initial diversity")
for a in axes[:,0]:
    a.set_ylabel ("Initial DOC (uM)")

#%%
null_data = null_data.dropna()
adapt_data = adapt_data.dropna()

fvalue, pvalue = stats.mannwhitneyu(null_data.ratio_t_50,adapt_data.ratio_t_50, alternative = "less")
print("ratio_t_50", fvalue, pvalue)

fvalue, pvalue = stats.mannwhitneyu(null_data.DOC_removal,adapt_data.DOC_removal, alternative = "less")
print("DOC_removal", fvalue, pvalue)

fvalue, pvalue = stats.mannwhitneyu(null_data.t_50_days,adapt_data.t_50_days, alternative = "less")
print("Time scale", fvalue, pvalue)
