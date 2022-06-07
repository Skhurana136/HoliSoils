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