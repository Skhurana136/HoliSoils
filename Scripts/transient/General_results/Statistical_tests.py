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

filename = os.path.join(results_dir, "combined_dataset.pkl")
diversity_data = pd.read_pickle(filename)
diversity_data = diversity_data.drop_duplicates()
diversity_data['DOC_initial_int'] = round(diversity_data.DOC_initial, -3)
diversity_data['S_initial_rnd'] = round(diversity_data.S_initial, 1)
diversity_data['ratio_t_50'] = diversity_data.T_50/diversity_data.T_50_B1
compl = diversity_data[diversity_data.Status>=0]
#%%
init_doc_list = list(compl.DOC_initial_int.unique())
activity_list = list(compl.activity.unique())

docres = []
doct50b1 = []
for doci in init_doc_list:
    for act in activity_list:
        subset = diversity_data[(diversity_data['DOC_initial_int']==doci) & (diversity_data['activity']==act)]
        subset_table1 = pd.pivot_table(subset, values = 'DOC_removal', index = ['carbon_species', 'Sim_series'], columns = ['S_initial_rnd'])
        fvalue, pvalue = stats.f_oneway(subset_table1[1.4],subset_table1[2.1],subset_table1[2.8],subset_table1[3.5])
        docres.append([doci, act, fvalue, pvalue])

        subset_table3 = pd.pivot_table(subset, values = 'DOC_T50_B1', index = ['carbon_species', 'Sim_series'], columns = ['S_initial_rnd'])
        fvalue, pvalue = stats.f_oneway(subset_table3[1.4],subset_table3[2.1],subset_table3[2.8],subset_table3[3.5])
        doct50b1.append([doci, act, fvalue, pvalue])

resdf = pd.DataFrame.from_records(docres, columns = ["DOC_initial", "activity", "fval", "pval"])
resdf.to_csv(os.path.join(results_dir, "anova_doc_removal.csv"))

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
