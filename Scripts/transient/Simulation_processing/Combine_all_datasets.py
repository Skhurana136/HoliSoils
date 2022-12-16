#%%
import os
import pandas as pd
import numpy as np
import itertools

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
results_dir = os.path.join(project_dir, "results")
sim_suffixes = ["_0_01", "_0_1", "_0_5", "", "_1_5x"]
sim_suffixes_var = [0.01, 0.1, 0.5, 1, 1.5]

fd_file = os.path.join(results_dir,"competition_adaptation_carbon_cue_wremaining_c_data_combined_dataset.pkl")
tim_file = os.path.join(results_dir,"competition_adaptation_carbon__decay_const_c_pools_data_initial_conditions_combined_dataset.pkl")
para_file = os.path.join(results_dir,"competition_adaptation_parameters_combined_dataset.pkl")
fd_data = pd.read_pickle(fd_file)
tim_data = pd.read_pickle(tim_file)
para_data = pd.read_pickle(para_file)

fd_tim_data = pd.merge(fd_data,tim_data, on = ["Seed", "Variance", "DOC_initial_int","biomass_species", "carbon_species", "Sim_series", "activity"], suffixes=('', '_y'))
fd_tim_data.drop(fd_tim_data.filter(regex='_y$').columns, axis=1, inplace=True)
all_data = pd.merge(fd_tim_data, para_data, on = ["Seed", "Variance", "biomass_species", "carbon_species", "Sim_series"],suffixes=('', '_y'))
all_data.drop(all_data.filter(regex='_y$').columns, axis=1, inplace=True)
print(all_data.columns)

all_data["active_H_c_connections"] = all_data['S_initial']*all_data['carbon_species']*all_data["activity"]/100
all_data["active_H"] = all_data['S_initial']*all_data["activity"]/100
all_data['FD_ratio'] = all_data.FD/all_data.FD_initial
all_data['Biomass_ratio'] = all_data.Biomass/all_data.Biomass_initial
all_data['decay_ratio'] = all_data.decay_const/all_data.Decay_const_initial
print(all_data.columns)

all_data.to_csv(os.path.join(project_dir,"simulation_results_temporal_initial_conditions_decay_const.csv"),index=False)