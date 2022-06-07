#%%
import os
from tkinter.ttk import Separator
import pandas as pd

import scipy.stats as stats

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
simulations_dir = os.path.join(project_dir, "simulations")
results_dir = os.path.join(project_dir, "results")
figures_dir = os.path.join(project_dir, "figures")

filename = os.path.join(results_dir, "null_cont_combined_dataset.pkl")
null_data = pd.read_pickle(filename)
null_data = null_data.drop_duplicates()
null_data['DOC_initial_int'] = round(null_data.DOC_initial, -3)
null_data['S_initial_rnd'] = round(null_data.S_initial, 1)
null_data['ratio_t_50'] = null_data.T_50/null_data.T_50_B1

filename = os.path.join(results_dir, "adaptation_cont_combined_dataset.pkl")
adapt_data = pd.read_pickle(filename)
adapt_data = adapt_data.drop_duplicates()
adapt_data['DOC_initial_int'] = round(adapt_data.DOC_initial, -3)
adapt_data['S_initial_rnd'] = round(adapt_data.S_initial, 1)
adapt_data['ratio_t_50'] = adapt_data.T_50/adapt_data.T_50_B1

#%%
init_doc_list = list(compl.DOC_initial_int.unique())
activity_list = list(compl.activity.unique())

docres = []
doct50b1 = []
docmidres = []
for doci in init_doc_list:
    for act in activity_list:
        subset = diversity_data[(diversity_data['DOC_initial_int']==doci) & (diversity_data['activity']==act)]
        subset_table1 = pd.pivot_table(subset, values = 'DOC_removal', index = ['carbon_species', 'Sim_series'], columns = ['S_initial_rnd'])
        fvalue, pvalue = stats.f_oneway(subset_table1[1.4],subset_table1[2.1],subset_table1[2.8],subset_table1[3.5])
        docres.append([doci, act, fvalue, pvalue])

        subset_table2 = pd.pivot_table(subset, values = 'DOC_removal_mid', index = ['carbon_species', 'Sim_series'], columns = ['S_initial_rnd'])
        fvalue, pvalue = stats.f_oneway(subset_table1[1.4],subset_table1[2.1],subset_table1[2.8],subset_table1[3.5])
        docmidres.append([doci, act, fvalue, pvalue])

        subset_table3 = pd.pivot_table(subset, values = 'DOC_T50_B1', index = ['carbon_species', 'Sim_series'], columns = ['S_initial_rnd'])
        fvalue, pvalue = stats.f_oneway(subset_table3[1.4],subset_table3[2.1],subset_table3[2.8],subset_table3[3.5])
        doct50b1.append([doci, act, fvalue, pvalue])

resdf = pd.DataFrame.from_records(docres, columns = ["DOC_initial", "activity", "fval", "pval"])
resdf.to_csv(os.path.join(results_dir, "anova_doc_removal.csv"))

resdf = pd.DataFrame.from_records(docmidres, columns = ["DOC_initial", "activity", "fval", "pval"])
resdf.to_csv(os.path.join(results_dir, "anova_doc_removal_mid.csv"))

resdf = pd.DataFrame.from_records(doct50b1, columns = ["DOC_initial", "activity", "fval", "pval"])
resdf.to_csv(os.path.join(results_dir, "anova_doc_t50_b1.csv"))


#%%
res = []
subst_1 = diversity_data#.dropna(subset = ['t_50_days']).reset_index()
init_doc_list = list(subst_1.DOC_initial_int.unique())
activity_list = list(subst_1.activity.unique())

for doci in init_doc_list:
    for act in activity_list:
        subset = subst_1[(subst_1['DOC_initial_int']==doci) & (subst_1['activity']==act)]
        subset_table = pd.pivot_table(subset, values = 't_50_days', index = ['carbon_species', 'Sim_series'], columns = ['S_initial_rnd'])
        if subset_table.shape[0]==0:
            pass
        else:
            fvalue, pvalue = stats.f_oneway(subset_table3[1.4],subset_table3[2.1],subset_table3[2.8],subset_table3[3.5])
            res.append([doci, act, fvalue, pvalue])
resdf = pd.DataFrame.from_records(res, columns = ["DOC_initial", "activity", "fval", "pval"])
resdf.to_csv(os.path.join(results_dir, "anova_char_rxn_time_scale.csv"))

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
import seaborn as sns
import matplotlib.pyplot as plt

#%%
resdf.loc[resdf["pval"]<0.1, "significance"] = 1
resdf.loc[resdf["pval"]>=0.1, "significance"] = 0
metric_to_plot = list(resdf.Metric.unique())
title_list = ["Impact on time", "Carbon removed", "Time scale", "Carbon removed at midpoint"]
resdf["x_var"] = resdf.S_initial*resdf.activity/100

fig, axes = plt.subplots(2,2, sharex = True, sharey = True)
for y, t in zip(metric_to_plot, title_list):
    subset = resdf[(resdf.Metric==y)&(resdf.significance==1)]
    axes.flatten()[metric_to_plot.index(y)].scatter("x_var", "DOC_initial", alpha = 0.7, data = subset)
    axes.flatten()[metric_to_plot.index(y)].set_title(t)


for a in axes[-1,:]:
    a.set_xlabel("% active initial diversity")
for a in axes[:,0]:
    a.set_ylabel ("Initial DOC (uM)")

plt.savefig(os.path.join(figures_dir, "anova_test.png"), dpi = 300)