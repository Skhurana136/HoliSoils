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

all_data = pd.read_pickle(os.path.join(results_dir, "cont_carbon_sw_off_combined_dataset.pkl"))
all_data['DOC_initial_int'] = round(all_data.DOC_initial, -3)
all_data['S_initial_rnd'] = round(all_data.S_initial, 1)
all_data['ratio_t_50'] = all_data.T_50/all_data.T_50_B1
all_data ['decay_const'] = 1/all_data.T_50
compl = all_data
init_doc_list = list(compl.DOC_initial_int.unique())
activity_list = list(compl.activity.unique())

#%%
docres = []
doct50b1 = []
docmidres = []
t50res = []
dec_const = []
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
        
        subset_table4 = pd.pivot_table(subset, values = 'decay_const', index = ['carbon_species', 'Sim_series'], columns = ['S_initial_rnd'])
        fvalue, pvalue = stats.f_oneway(subset_table1[1.4],subset_table1[2.1],subset_table1[2.5], subset_table1[2.8],subset_table1[3.0],subset_table1[3.2],subset_table1[3.3],subset_table1[3.5])
        dec_const.append([doci, act, fvalue, pvalue])

resdf = pd.DataFrame.from_records(docres, columns = ["DOC_initial", "activity", "fval", "pval"])
resdf.to_csv(os.path.join(results_dir, "sw_off_anova_doc_removal_null.csv"))

resdf = pd.DataFrame.from_records(docmidres, columns = ["DOC_initial", "activity", "fval", "pval"])
resdf.to_csv(os.path.join(results_dir, "sw_off_anova_doc_removal_mid_null.csv"))

resdf = pd.DataFrame.from_records(doct50b1, columns = ["DOC_initial", "activity", "fval", "pval"])
resdf.to_csv(os.path.join(results_dir, "sw_off_anova_doc_t50_b1_null.csv"))

resdf = pd.DataFrame.from_records(t50res, columns = ["DOC_initial", "activity", "fval", "pval"])
resdf.to_csv(os.path.join(results_dir, "sw_off_anova_t_50_null.csv"))

resdf = pd.DataFrame.from_records(dec_const, columns = ["DOC_initial", "activity", "fval", "pval"])
resdf.to_csv(os.path.join(results_dir, "sw_off_anova_decay_const_null.csv"))