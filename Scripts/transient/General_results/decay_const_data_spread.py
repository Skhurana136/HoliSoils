#%%
import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
simulations_dir = os.path.join(project_dir, "simulations")
results_dir = os.path.join(project_dir, "results")
figures_dir = os.path.join(project_dir, "figures")

filename = os.path.join(results_dir, "null_cont_combined_dataset.pkl")
diversity_data = pd.read_pickle(filename)
diversity_data = diversity_data.drop_duplicates()
diversity_data['DOC_initial_int'] = round(diversity_data.DOC_initial, -3)
diversity_data['decay_const_base'] = 1/diversity_data.T_50_B1
diversity_data['decay_const'] = 1/diversity_data.T_50
diversity_data['decay_ratio'] = diversity_data.decay_const/diversity_data.decay_const_base
diversity_data['S_initial_int'] = round(diversity_data.S_initial, 1)
compl = diversity_data.dropna(subset = ['decay_ratio'])
init_doc_list = list(compl.DOC_initial_int.unique())
activity_list = list(compl.activity.unique())
s_initial_list = list(compl.S_initial_int.unique())

#%%
incompl = diversity_data[diversity_data.decay_ratio.isna()]
incompl["sim_folder"] = "cont_carbon_" + incompl.carbon_species.astype(str) + "_" + incompl.Seed.astype(str) + "_ip_0"
sim_fol_list = list(incompl.sim_folder.unique())
incompl_vmax_working = {}
incompl_sp_num = {}
incompl_ks_working = {}
for folder in sim_fol_list:
    filename = os.path.join(simulations_dir, folder, "seeds_randoms.pkl")
    seed_details = pd.read_pickle(filename)
    
    subset = incompl[incompl.sim_folder == folder]
    bio_sp = list(subset.biomass_species)
    seed_list = list(subset.Seed)
    c = list(subset.carbon_species.unique())[0]
    
    for see in seed_list:
        seed_set = subset[subset.Seed==see]
        bio_sp = list(seed_set.biomass_species)
        sim_list = list(seed_set.Sim_series)
        doc_list = list(seed_set.DOC_initial_int)

        for b, sim_scenario, doc in zip(bio_sp, sim_list, doc_list):
            dic_key = sim_scenario + "/bio_n_" + str(b) + "/dom_initial_" + str(int(doc)) + "/seed_" + str(see)
            base_dic_key = "b_1_all_/bio_n_" + str(b) + "/dom_initial_" + str(int(doc)) + "/seed_" + str(see)
            vmax = seed_details[dic_key]['max_rate_parameters']
            ks = seed_details[dic_key]['sat_const_parameters']
            base_vmax = seed_details[base_dic_key]['max_rate_parameters']
            base_ks = seed_details[base_dic_key]['sat_const_parameters']
            switched_on = np.where((vmax- base_vmax)==0)
            vmax_on = base_vmax[switched_on]
            ks_on = base_ks[switched_on]
            base_v_c_max = np.max(base_vmax, axis = 0)
            base_v_c_min = np.min(base_vmax, axis = 0)
            v_c_max = np.max(vmax_on, axis = 0)
            v_c_min = np.min(vmax_on, axis = 0)
            ks_on_c_max = np.max(ks_on, axis = 0)
            ks_on_c_min = np.min(ks_on, axis = 0)    
            incompl_ks_working[dic_key+"/carbon_" + str(c)] = {'max': ks_on_c_max, 'min' : ks_on_c_min}
            incompl_vmax_working[dic_key+"/carbon_" + str(c)] = {'max': (base_v_c_max-v_c_max)/base_v_c_max}#, 'min' : }
            incompl_sp_num[dic_key+"/carbon_" + str(c)] = {'bio': b, 'carb' : c}

compl["sim_folder"] = "cont_carbon_" + compl.carbon_species.astype(str) + "_" + compl.Seed.astype(str) + "_ip_0"
work_sim_fol_list = list(compl.sim_folder.unique())
compl_vmax_working = {}
compl_sp_num = {}
compl_ks_working = {}
for folder in work_sim_fol_list:
    filename = os.path.join(simulations_dir, folder, "seeds_randoms.pkl")
    seed_details = pd.read_pickle(filename)
    
    subset = compl[compl.sim_folder == folder]
    bio_sp = list(subset.biomass_species)
    seed_list = list(subset.Seed)
    c = list(subset.carbon_species.unique())[0]
    
    for see in seed_list:
        seed_set = subset[subset.Seed==see]
        bio_sp = list(seed_set.biomass_species)
        sim_list = list(seed_set.Sim_series)
        doc_list = list(seed_set.DOC_initial_int)

        for b, sim_scenario, doc in zip(bio_sp, sim_list, doc_list):
            dic_key = sim_scenario + "/bio_n_" + str(b) + "/dom_initial_" + str(int(doc)) + "/seed_" + str(see)
            base_dic_key = "b_1_all_/bio_n_" + str(b) + "/dom_initial_" + str(int(doc)) + "/seed_" + str(see)
            vmax = seed_details[dic_key]['max_rate_parameters']
            ks = seed_details[dic_key]['sat_const_parameters']
            base_vmax = seed_details[base_dic_key]['max_rate_parameters']
            base_ks = seed_details[base_dic_key]['sat_const_parameters']
            ox_state = seed_details[base_dic_key]['oxidation_state']
            switched_on = np.where((vmax- base_vmax)==0)
            vmax_on = base_vmax[switched_on]
            ks_on = base_ks[switched_on]
            base_v_c_max = np.max(base_vmax, axis = 0)
            base_v_c_min = np.min(base_vmax, axis = 0)
            v_c_max = np.max(vmax_on, axis = 0)
            v_c_min = np.min(vmax_on, axis = 0)
            ks_on_c_max = np.max(ks_on, axis = 0)
            ks_on_c_min = np.min(ks_on, axis = 0)    
            compl_ks_working[dic_key+"/carbon_" + str(c)] = {'max': ks_on_c_max, 'min' : ks_on_c_min}
            compl_vmax_working[dic_key+"/carbon_" + str(c)] = {'max': (base_v_c_max-v_c_max)/base_v_c_max}#, 'min' : }
            compl_sp_num[dic_key+"/carbon_" + str(c)] = {'bio': b, 'carb' : c}
            
    
#%%
nw_vm_list = []
for all_keys in list(incompl_vmax_working.keys()):
    nw_vm_list.append(incompl_vmax_working[all_keys]['max'])
w_vm_list = []
for all_keys in list(compl_vmax_working.keys()):
    w_vm_list.append(compl_vmax_working[all_keys]['max'])
sns.histplot(w_vm_list, color = 'orange', bins = 30, label = "worked")#, kde = True, fill = False)
#sns.histplot(nw_vm_list, bins = 30, label = "not worked")#, kde = True, fill = False)
#sns.kdeplot(w_vm_list, color = 'orange', label = "worked")
#sns.kdeplot(nw_vm_list, color = "steelblue", label = "not worked")
plt.legend()
    
#%%
nw_vm_list = []
w_vm_list = []
for all_keys in list(incompl_vmax_working.keys()):
    nw_vm_list.append(incompl_vmax_working[all_keys]['min'])
for all_keys in list(compl_vmax_working.keys()):
    w_vm_list.append(compl_vmax_working[all_keys]['min'])
#sns.histplot(w_vm_list, color = 'orange', bins = 30, label = "worked", kde = True, fill = False)
#sns.histplot(nw_vm_list, bins = 30, label = "not worked", kde = True, fill = False)
sns.kdeplot(w_vm_list, color = 'orange', label = "worked")
sns.kdeplot(nw_vm_list, color = "steelblue", label = "not worked")
plt.xlim(right = 4)
#plt.yscale("log")
plt.legend()

#%%
nw_c_list = []
w_c_list = []
for all_keys in list(incompl_sp_num.keys()):
    nw_c_list.append(incompl_sp_num[all_keys]['carb']/incompl_sp_num[all_keys]['bio'])
for all_keys in list(compl_sp_num.keys()):
    w_c_list.append(compl_sp_num[all_keys]['carb']/compl_sp_num[all_keys]['bio'])
sns.histplot(w_c_list, color = 'orange', bins = 30, label = "worked", kde = True, fill = False)
sns.histplot(nw_c_list, bins = 30, label = "not worked", kde = True, fill = False)
#sns.kdeplot(w_vm_list, color = 'orange', label = "worked")
#sns.kdeplot(nw_vm_list, color = "steelblue", label = "not worked")
plt.legend()

#%%
nw_ks_list = []
w_ks_list = []
for all_keys in list(incompl_ks_working.keys()):
    nw_ks_list.append(incompl_ks_working[all_keys]['max']/incompl_ks_working[all_keys]['min'])
for all_keys in list(compl_ks_working.keys()):
    w_ks_list.append(compl_ks_working[all_keys]['max']/compl_ks_working[all_keys]['min'])
#sns.histplot(w_ks_list, color = 'orange', bins = 30, label = "worked", kde = True, fill = False)
#sns.histplot(nw_ks_list, bins = 30, label = "not worked", kde = True, fill = False)
sns.kdeplot(w_ks_list, color = 'orange', label = "worked")
sns.kdeplot(nw_ks_list, color = "steelblue", label = "not worked")
plt.legend()