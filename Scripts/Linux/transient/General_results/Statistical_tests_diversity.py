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

a_n_data = [null_data, adapt_data]
all_data = pd.concat(a_n_data)
#%%
compl = null_data
init_doc_list = list(compl.DOC_initial_int.unique())
activity_list = list(compl.activity.unique())

#%%
docres = []
doct50b1 = []
docmidres = []
t50res = []
for doci in init_doc_list:
    for act in activity_list:
        subset = compl[(compl['DOC_initial_int']==doci) & (compl['activity']==act)]
        subset_table1 = pd.pivot_table(subset, values = 'DOC_removal', index = ['carbon_species', 'Sim_series'], columns = ['S_initial_rnd'])
        fvalue, pvalue = stats.f_oneway(subset_table1[1.4],subset_table1[2.1],subset_table1[2.5], subset_table1[2.8],subset_table1[3.0],subset_table1[3.2],subset_table1[3.3],subset_table1[3.5])
        docres.append([doci, act, fvalue, pvalue])

        subset_table2 = pd.pivot_table(subset, values = 'DOC_removal_mid', index = ['carbon_species', 'Sim_series'], columns = ['S_initial_rnd'])
        fvalue, pvalue = stats.f_oneway(subset_table1[1.4],subset_table1[2.1],subset_table1[2.5], subset_table1[2.8],subset_table1[3.0],subset_table1[3.2],subset_table1[3.3],subset_table1[3.5])
        docmidres.append([doci, act, fvalue, pvalue])

        subset_table3 = pd.pivot_table(subset, values = 'DOC_T50_B1', index = ['carbon_species', 'Sim_series'], columns = ['S_initial_rnd'])
        fvalue, pvalue = stats.f_oneway(subset_table1[1.4],subset_table1[2.1],subset_table1[2.5], subset_table1[2.8],subset_table1[3.0],subset_table1[3.2],subset_table1[3.3],subset_table1[3.5])
        doct50b1.append([doci, act, fvalue, pvalue])

        subset_table4 = pd.pivot_table(subset, values = 'T_50', index = ['carbon_species', 'Sim_series'], columns = ['S_initial_rnd'])
        fvalue, pvalue = stats.f_oneway(subset_table1[1.4],subset_table1[2.1],subset_table1[2.5], subset_table1[2.8],subset_table1[3.0],subset_table1[3.2],subset_table1[3.3],subset_table1[3.5])
        t50res.append([doci, act, fvalue, pvalue])

resdf = pd.DataFrame.from_records(docres, columns = ["DOC_initial", "activity", "fval", "pval"])
resdf.to_csv(os.path.join(results_dir, "anova_doc_removal_null.csv"))

resdf = pd.DataFrame.from_records(docmidres, columns = ["DOC_initial", "activity", "fval", "pval"])
resdf.to_csv(os.path.join(results_dir, "anova_doc_removal_mid_null.csv"))

resdf = pd.DataFrame.from_records(doct50b1, columns = ["DOC_initial", "activity", "fval", "pval"])
resdf.to_csv(os.path.join(results_dir, "anova_doc_t50_b1_null.csv"))

resdf = pd.DataFrame.from_records(t50res, columns = ["DOC_initial", "activity", "fval", "pval"])
resdf.to_csv(os.path.join(results_dir, "anova_t_50_null.csv"))