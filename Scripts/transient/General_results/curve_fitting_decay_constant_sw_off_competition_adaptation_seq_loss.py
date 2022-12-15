#%%
import os
import pandas as pd
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import h5py

def proc_datasets (df):
    print(df.columns)
    df['DOC_initial_int'] = round(df.DOC_initial, -3)
    df['S_initial_int'] = round(df.S_initial, 1)
    df['ratio_t_50'] = df.T_50/df.T_50_B1
    df["decay_const"] = 1/df.t_50_days
    return df

## LOAD RESULTS

project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient","seq_loss")
simulations_dir = os.path.join(project_dir, "simulations")
results_dir = os.path.join(project_dir, "results")
figures_dir = os.path.join(project_dir, "figures")

data_8sp = pd.read_pickle(os.path.join(results_dir, "competition_adaptation_seq_8_carbon__loss_0.9_combined_dataset.pkl"))
data_4sp = pd.read_pickle(os.path.join(results_dir, "competition_adaptation_seq_carbon__loss_0.9_combined_dataset.pkl"))

data_8sp = proc_datasets(data_8sp)
data_4sp = proc_datasets(data_4sp)
all_data = pd.concat([data_8sp, data_4sp])

compl = all_data.dropna(subset = ['ratio_t_50'])
init_doc_list = np.sort(list(compl.DOC_initial_int.unique()))
seed_list = list(compl.Seed.unique())
carbon_list = list(compl.carbon_species.unique())
biomass_list = list(compl.biomass_species.unique())
activity_list = np.sort(list(compl.activity.unique()))
s_initial_list = np.sort(list(compl.S_initial_int.unique()))
non_100 = compl[compl.activity<100]
generalist_act = compl[compl.activity==100].reset_index()

## PREDICT FUNCTION PARAMETERS FOR IMPACT ON DECAY CONSTANT ##
#%%
print(generalist_act.t_50_days.min(), generalist_act.t_50_days.max())
print(generalist_act.S_initial.min(), generalist_act.S_initial.max())
#%%
h = sns.violinplot(x = compl.DOC_initial_int, y = compl.decay_const, hue = compl.carbon_species)
plt.xlabel("Initial available carbon (uM)")
plt.ylabel("Decay constant")
#%%
h = sns.violinplot(x = compl.DOC_initial_int, y = compl.decay_const, hue = compl.biomass_species)
plt.xlabel("Initial available carbon (uM)")
plt.ylabel("Decay constant")

#%%
h = sns.violinplot(x = generalist_act.DOC_initial_int, y = generalist_act.decay_const, hue = generalist_act.carbon_species)
plt.xlabel("Initial available carbon (uM)")
plt.ylabel("Decay constant")
#%%
h = sns.violinplot(x = generalist_act.DOC_initial_int, y = generalist_act.decay_const, hue = generalist_act.biomass_species)
plt.xlabel("Initial available carbon (uM)")
plt.ylabel("Decay constant")


#%%
plt.figure(figsize = (8,4))
h = sns.boxplot(x = generalist_act.DOC_initial_int, y = generalist_act.decay_const, hue = generalist_act.S_initial_int)
plt.xlabel("Initial available carbon (uM)")
plt.ylabel("Decay constant")
#%%
h = sns.boxplot(x = generalist_act.carbon_species, y = generalist_act.decay_const, hue = generalist_act.DOC_initial_int)
plt.xlabel("Shannon diversity")
plt.ylabel("Decay constant")
#%%
h = sns.boxplot(x = generalist_act.biomass_species, y = generalist_act.decay_const, hue = generalist_act.DOC_initial_int)
plt.xlabel("Shannon diversity")
plt.ylabel("Decay constant")
#%%
plt.figure(figsize = (6,2))
h = sns.violinplot(x = generalist_act.carbon_species, y = generalist_act.decay_const,  hue = generalist_act.biomass_species)
plt.xlabel("Shannon diversity")
plt.ylabel("Decay constant")
#%%
def gather_paras(sim_data):
    dom_n, bio_n = sim_data['species_number']
    init_bio = sim_data['initial_conditions']['biomass']
    init_bio_sum = np.sum(init_bio)
    init_c = sim_data['initial_conditions']['dom']
    init_c_sum = np.sum(init_c)
    paras = sim_data['parameters']
    dying = np.tile(init_bio*np.asarray(paras['mortality_const'])/init_bio_sum, (dom_n, 1)).flatten().reshape(-1,1)
    exo_enz = np.tile(init_bio*np.asarray(paras['exo_enzyme_rate'])/init_bio_sum, (dom_n, 1)).flatten().reshape(-1,1)
    v_params = (init_bio*np.asarray((paras['max_rate']))/init_bio_sum).flatten().reshape(-1,1)*np.asarray(paras['carbon_uptake']).flatten().reshape(-1,1)*exo_enz
    k_params = np.asarray(paras['half_saturation']).flatten().reshape(-1,1)
    #ox = np.asarray(list([x*init_c[idx]/init_c_sum]*bio_n for x,idx in enumerate(paras['oxidation_state']))).flatten().reshape(-1,1)
    ox = np.asarray(list([x]*bio_n for x in paras['oxidation_state'])).flatten().reshape(-1,1)
    paras_arr = np.append(np.append(np.append(dying, exo_enz, axis=1),np.append(v_params, k_params, axis = 1), axis = 1), ox, axis = 1)
    return paras_arr

def run_gather_params(base_list, base_list_act):
    seed_sim_list = generalist_act.Seed.unique().tolist()
    cn_list = [3,6,12,18]
    bio_n_series = [4,8]
    init_dom_list = [1000,2000,5000,10000,15000]
    filestring = ["competition_adaptation_seq_carbon_","competition_adaptation_seq_8_carbon_"]
    all_para_arr = np.zeros((0,10))
    for c_n in cn_list:
        row = []
        c_b_row = []
        #results_filename = os.path.join(results_dir, filestring + str(c_n))

        for seed_sim in seed_sim_list:
            # Load all datasets and save their Shannon and diversity indices in a dataframe
            seed_all = 'seed_'+str(seed_sim)
            
            for b_n in bio_n_series:
                details_subfolder = filestring[bio_n_series.index(b_n)] + str(c_n) + '_'+str(seed_sim) + '_ip_0'
                simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
                hrw = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'r+')
                for t_dom_initial in init_dom_list:
                    c_b = "bio_n_"+ str(b_n)
                    dom_init = "dom_initial_" + str(t_dom_initial)
                    doc_input = (t_dom_initial)
                    seed_arr = np.zeros((c_n*b_n,1))+seed_sim
                    bio_arr = np.zeros((c_n*b_n,1))+b_n
                    seed_arr_bio = np.append(seed_arr, bio_arr, axis = 1) 
                    carbon_arr = np.zeros((c_n*b_n,1))+c_n
                    seed_carbon = np.append(seed_arr_bio, carbon_arr, axis = 1)
                    for base, act in zip(base_list, base_list_act):#zip(["b_1", "b_2", "b_3", "b_4", "b_5"], [1.0, 0.0, 0.25, 0.5, 0.75]):
                        act_arr = np.zeros((c_n*b_n,1))+act
                        seed_sim_arr = np.append(seed_carbon, act_arr, axis = 1)
                        if base == "b_1":
                            base_case = "b_1_all_"
                            level_arr = np.zeros((c_n*b_n,1))+0
                            sim_data = hrw[base_case][c_b][dom_init][seed_all]
                            seed_sim_level = np.append(seed_sim_arr, level_arr, axis = 1)
                            paras_sim = np.append(seed_sim_level, gather_paras(sim_data), axis = 1)
                            all_para_arr = np.append(all_para_arr, paras_sim, axis = 0)
                        else:
                            for label,label_id in zip(["a", "b", "c","d","e"], [1,2,3,4,5]):
                                sim = base + "_" + label + "_"
                                level_arr = np.zeros((c_n*b_n,1))+label_id
                                sim_data = hrw[sim][c_b][dom_init][seed_all]
                                seed_sim_level = np.append(seed_sim_arr, level_arr, axis = 1)
                                paras_sim = np.append(seed_sim_level, gather_paras(sim_data), axis = 1)
                                all_para_arr = np.append(all_para_arr, paras_sim, axis = 0)
            hrw.close()
    para_df = pd.DataFrame(all_para_arr, columns = ['Seed', 'biomass_species','carbon_species','Activity', 'level_id', 'm', 'exo_enz','vmax', 'k', 'oxidation_state'])
    para_df.loc[para_df["level_id"] == 0.0, "Sim_series"] = 'b_1_all_'
    para_df.loc[para_df["level_id"] == 1.0, "Sim_series"] = base_list[0]+'_a_'
    para_df.loc[para_df["level_id"] == 2.0, "Sim_series"] = base_list[0]+'_b_'
    para_df.loc[para_df["level_id"] == 3.0, "Sim_series"] = base_list[0]+'_c_'
    para_df.loc[para_df["level_id"] == 4.0, "Sim_series"] = base_list[0]+'_d_'
    para_df.loc[para_df["level_id"] == 5.0, "Sim_series"] = base_list[0]+'_e_'
    return para_df
#%%
base_para_df = run_gather_params(["b_1"], [1.0])
on1_para_df = run_gather_params(["b_2"], [0.1])
on2_para_df = run_gather_params(["b_3"], [0.25])
on3_para_df = run_gather_params(["b_4"], [0.5])
on4_para_df = run_gather_params(["b_5"], [0.75])

#%%
para_df = pd.concat([base_para_df,on1_para_df, on2_para_df, on3_para_df, on4_para_df])

#%%
func_div = para_df.groupby(['Seed', 'biomass_species', 'carbon_species', 'Sim_series', 'level_id'], as_index = False).agg({
    'vmax': ['mean','median','std','skew'],
    'k':['mean','median','std','skew'],
    'm':['mean','median','std','skew'],
    'exo_enz':['mean','median','std','skew']
}).reset_index()
func_div.columns = func_div.columns.map('_'.join)
func_div.rename(columns = {'Seed_':'Seed', 'carbon_species_':'carbon_species', 'biomass_species_': 'biomass_species', 'Sim_series_': 'Sim_series', 'level_id_':'level_id'}, inplace = True)
#%%
fig, axes = plt.subplots(3,3, sharex = 'col', figsize = (12,6))
sns.scatterplot(data = func_div, x = 'vmax_mean', y = 'vmax_median', hue = 'biomass_species', ax = axes[0,0])
axes[0,0].set_title("vmax")
sns.scatterplot(data = func_div, x = 'vmax_mean', y = 'vmax_std', hue = 'biomass_species', ax = axes[1,0])
sns.scatterplot(data = func_div, x = 'vmax_mean', y = 'vmax_skew', hue = 'biomass_species', ax = axes[2,0])
sns.scatterplot(data = func_div, x = 'k_mean', y = 'k_median', hue = 'biomass_species', ax = axes[0,1])
axes[0,1].set_title("Half-sat-const")
sns.scatterplot(data = func_div, x = 'k_mean', y = 'k_std', hue = 'biomass_species', ax = axes[1,1])
sns.scatterplot(data = func_div, x = 'k_mean', y = 'k_skew', hue = 'biomass_species', ax = axes[2,1])
sns.scatterplot(data = func_div, x = 'm_mean', y = 'm_median', hue = 'biomass_species', ax = axes[0,2])
axes[0,2].set_title("Mortality")
sns.scatterplot(data = func_div, x = 'm_mean', y = 'm_std', hue = 'biomass_species', ax = axes[1,2])
sns.scatterplot(data = func_div, x = 'm_mean', y = 'm_skew', hue = 'biomass_species', ax = axes[2,2])

for ax in axes.flatten()[:-1]:
    ax.legend([],[], frameon=False)
for ax in axes.flatten()[:]:
    ax.set_ylabel('')
    ax.set_xlabel('')
for ax in axes[2,:]:
    ax.set_xlabel("Mean")
axes[0,0].set_ylabel("Median")
axes[1,0].set_ylabel("Sdev")
axes[2,0].set_ylabel("Skew")
plt.figure()
sns.scatterplot(data = func_div, x = "vmax_mean", y= "k_mean", hue = "carbon_species")
#%%
#Compare baselines with switched off activities
for s in seed_list:
    for c in carbon_list:
        for b in biomass_list:
            sub_scb = func_div[(func_div["Seed"].astype(int)==s)&(func_div["carbon_species"].astype(int)==c)&(func_div["biomass_species"].astype(int)==b)]
            v_mean_base = sub_scb[sub_scb.Sim_series=='b_1_all_']['vmax_mean'].values[0]
            func_div.loc[(func_div["Seed"].astype(int)==s)&(func_div["carbon_species"].astype(int)==c)&(func_div["biomass_species"].astype(int)==b), "vmax_mean_base"]=v_mean_base

#%%
#seed_paras = pd.merge(generalist_act, func_div[["Seed", "biomass_species", "carbon_species", "Sim_series", "level_id","m_skew","exo_enz_skew", "vmax_skew", "k_skew","m_mean","exo_enz_mean", "vmax_mean", "k_mean"]], on = ["Seed", "biomass_species", "carbon_species", "Sim_series"])
#cols_to_merge = ["Seed", "biomass_species", "carbon_species", "Sim_series", "level_id","m_skew","exo_enz_skew", "vmax_skew", "k_skew","m_mean","exo_enz_mean", "vmax_mean", "k_mean"]
cols_to_merge = func_div.columns.to_list()
seed_paras = pd.merge(all_data, func_div[cols_to_merge], on = ["Seed", "biomass_species", "carbon_species", "Sim_series"])
#%%
seed_paras["docperc"] = seed_paras.DOC_initial/seed_paras.carbon_species
seed_paras["x_var"] = seed_paras.DOC_initial*seed_paras.carbon_species*seed_paras.vmax_mean##/seed_paras.k_mean
#%%
plt.figure()
sns.scatterplot(y = 'decay_const', x = "x_var", style = "carbon_species", data = seed_paras)
plt.xlabel("Available carbon x carbon species number x average vmax x exo_emzyme rate")
plt.legend( bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.figure()
sns.scatterplot(y = 'decay_const', x = "x_var", hue = "biomass_species", data = seed_paras)
plt.xlabel("Available carbon x carbon species number x average vmax x exo_emzyme rate")
plt.legend( bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
#%%
sns.scatterplot(data = seed_paras, x = "vmax_mean", y= "k_mean", hue = "biomass_species")
#%%
sns.scatterplot(data = seed_paras, x = "vmax_mean_base", y = "vmax_mean", hue = "activity", style= "S_initial_int")
plt.legend(bbox_to_anchor = (1.02,1))
#%%
seed_paras["reldelvmax"] = seed_paras.vmax_mean/seed_paras.vmax_mean_base
sub_off = seed_paras[seed_paras['Sim_series']!='b_1_all_']
sub_off["ratio_xvar"] = (1 - sub_off.reldelvmax)
sub_off["yplot"]  = 1/sub_off.ratio_t_50
sns.scatterplot(y = 'yplot', x = "ratio_xvar", hue = "DOC_initial_int", style = "activity",data = sub_off)
plt.xlabel("Reduction in decomposition potential (difference in vmean)\nnormalized by decomposition potential of base case(vmean)")
plt.ylabel("Decay constant normalized by\nthat of the base case")
plt.legend( bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
#%%
sns.scatterplot(y = 'decay_const', x = "x_var", hue = "DOC_initial_int", style = "biomass_species",data = seed_paras)
plt.xlabel("Available carbon x carbon species number x biomass weighted vmax")
plt.ylabel("Decay constant")
plt.legend( bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
#%%
sns.scatterplot(y = "ratio_xvar", x = "activity", hue = "DOC_initial_int", style = "carbon_species", data = sub_off)
plt.ylabel("Reduction in decomposition potential (difference in vmean)\nnormalized by decomposition potential of base case(vmean)")
plt.legend( bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)