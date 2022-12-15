import os
import pandas as pd
import numpy as np
import itertools

## LOAD RESULTS
#project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
project_dir = os.path.join('/proj', 'hs_micro_div_072022', 'Project_data', 'transient')
results_dir = os.path.join(project_dir, "results")
sim_suffixes = ["_0_01", "_0_1", "_0_5", "", "_1_5x"]
sim_suffixes_var = [0.01, 0.1, 0.5, 1, 1.5]

fd_file = os.path.join(results_dir,"competition_adaptation_carbon_cue_wremaining_c_data_combined_dataset.pkl")
tim_file = os.path.join(results_dir,"competition_adaptation_carbon__decay_const_c_pools_data_initial_conditions_combined_dataset.pkl")
para_file = os.path.join(results_dir,"competition_adaptation_parameters_combined_dataset.pkl")
fd_data = pd.read_pickle(fd_file)
tim_data = pd.read_pickle(tim_file)
para_data = pd.read_pickle(para_file)

#fd_tim_data = pd.merge(fd_data, tim_data, on = ["Seed", "Variance", "DOC_initial_int","biomass_species", "carbon_species", "Sim_series", "activity"])

#cols_to_merge = para_data.columns.to_list()[2:]
fd_tim_data = pd.merge(fd_data,tim_data, on = ["Seed", "Variance", "DOC_initial_int","biomass_species", "carbon_species", "Sim_series", "activity"], suffixes=('', '_y'))

fd_tim_data.drop(fd_tim_data.filter(regex='_y$').columns, axis=1, inplace=True)
all_data = pd.merge(fd_tim_data, para_data, on = ["Seed", "Variance", "biomass_species", "carbon_species", "Sim_series"],suffixes=('', '_y'))
all_data.drop(all_data.filter(regex='_y$').columns, axis=1, inplace=True)
all_data['FD_ratio'] = all_data.FD/all_data.FD_initial
all_data['Biomass_ratio'] = all_data.Biomass/all_data.Biomass_initial
all_data["active_H_c_connections"] = all_data['S_initial']*all_data['carbon_species']*all_data["activity"]/100
all_data["active_H"] = all_data['S_initial']*all_data["activity"]/100
print(all_data.columns)

init_doc_list = np.sort(list(all_data.DOC_initial_int.unique()))
#activity_list = np.sort(list(all_data.activity_x.unique()))
seed_sim_list = np.sort(list(all_data.Seed.unique()))
c_sp_list = list(all_data.carbon_species.unique())
bio_sp_list = list(all_data.biomass_species.unique())
sim_list = list(all_data.Sim_series.unique())

all_data['Decay_const'] = 1/all_data['T_loss']
initial_data = all_data[all_data['%C']==100.]
initial_data.rename(columns={'Decay_const':'Decay_const_initial'}, inplace=True)
tim_initial_data = pd.merge(all_data, initial_data, on = ["Seed", "Variance", "biomass_species", "carbon_species", "Sim_series", "C_pool"],suffixes=('', '_y'))
tim_initial_data.drop(tim_initial_data.filter(regex='_y$').columns, axis=1, inplace=True)
#tim_initial_data.to_csv(os.path.join("C:\Users\swkh9804\Documents\Project_data\HoliSoilssimulation_results_temporal_initial_conditions_decay_const_modified.csv"), index=False)

#for v,s,c,b,d,sim in itertools.product([0.01, 0.1, 0.5, 1., 1.5],seed_sim_list, c_sp_list,bio_sp_list,init_doc_list,sim_list):
#    sub_scb = all_data[(all_data["Variance"]==v)&(all_data["Seed"].astype(int)==s)&(all_data["carbon_species"].astype(int)==c)&(all_data["biomass_species"].astype(int)==b)&(all_data["DOC_initial_int"].astype(int)==d)&(all_data["Sim_series"]==sim)]
#    decay_initial = sub_scb[sub_scb['%C']==100.]['decay_const'].values[0]
#    all_data.loc[(all_data["Variance"]==v)&(all_data["Seed"].astype(int)==s)&(all_data["carbon_species"].astype(int)==c)&(all_data["biomass_species"].astype(int)==b)&(all_data["DOC_initial_int"].astype(int)==d)&(all_data["Sim_series"]==sim), "Decay_const_initial"]=decay_initial
#all_data["FD_cov"]= np.sqrt(all_data.FD_initial)/all_data.vmax_mean
#all_data['vmax_ratio'] = all_data.vmax_mean/all_data.vmax_mean_base
all_data['decay_ratio'] = all_data.decay_const/all_data.Decay_const_initial
all_data.to_csv(os.path.join(project_dir,"simulation_results_temporal_initial_conditions_decay_const.csv"),index=False)