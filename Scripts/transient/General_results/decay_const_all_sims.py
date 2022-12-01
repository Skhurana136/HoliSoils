#%%
import os
import pandas as pd
import numpy as np

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
sim_suffixes = ["_0_01", "_0_1", "_0_5", "", "_1_5x"]
sim_suffixes_var = [0.01, 0.1, 0.5, 1, 1.5]
#%%
files = []
for p,a in zip(sim_suffixes, sim_suffixes_var):
    filename = os.path.join(project_dir, "gen_spec_lognorm" + p, "results", "competition_adaptation_carbon__loss_0.9_cue_combined_dataset.pkl")
    data = pd.read_pickle(filename)
    var_arr = np.zeros((data.shape[0],))+a
    var_ser = pd.Series(var_arr, copy=False,name = "Variance")
    data_var = data.join(var_ser)
    files.append(data_var)
fd_data = pd.concat(files)
print(fd_data.columns)
fd_data['DOC_initial_int'] = round(fd_data.DOC_initial, -3)
#%%
### CHARACTERISTIC REACTION TIME SCALE
files = []
for p,a in zip(sim_suffixes, sim_suffixes_var):
    filename = os.path.join(project_dir, "gen_spec_lognorm" + p, "results", "competition_adaptation_carbon_loss_temporal_temporal_decay_const_c_pools_combined_dataset.pkl")
    data = pd.read_pickle(filename)
    var_arr = np.zeros((data.shape[0],))+a
    var_ser = pd.Series(var_arr, copy=False,name = "Variance")
    data_var = data.join(var_ser)
    files.append(data_var)
tim_data = pd.concat(files)
print(tim_data.columns)
tim_data['DOC_initial_int'] = round(tim_data.DOC_initial, -3)
#%%
cols_to_merge = ['Seed', 'Sim_series', 'carbon_species', 'biomass_species', 'C_pool', 'T10', 'T10_base', 'T20', 'T20_base', 'T30', 'T30_base', 'activity', 'Variance',  'DOC_initial_int']
fd_tim_data = pd.merge(fd_data, tim_data[cols_to_merge], on = ["Seed", "Variance", "DOC_initial_int","biomass_species", "carbon_species", "Sim_series", "activity"])
files = []
for p,a in zip(sim_suffixes, sim_suffixes_var):
    filename = os.path.join(project_dir, "gen_spec_lognorm" + p, "results", "gen_spec_lognorm" + p + "_parameters.csv")
    data = pd.read_csv(filename)
    var_arr = np.zeros((data.shape[0],))+a
    var_ser = pd.Series(var_arr, copy=False,name = "Variance")
    data_var = data.join(var_ser)
    files.append(data_var)
para_data = pd.concat(files)
print(para_data.columns)
cols_to_merge = para_data.columns.to_list()[2:]
all_data = pd.merge(fd_tim_data, para_data[cols_to_merge], on = ["Seed", "Variance", "biomass_species", "carbon_species", "Sim_series"])
all_data["active_H_c_connections"] = all_data['S_initial']*all_data['carbon_species']*all_data["activity"]/100
all_data["Decay_constant_10"] = 1/(all_data['T10']*5).to_numpy(dtype = float)    
all_data["Decay_constant_20"] = 1/(all_data['T20']*5).to_numpy(dtype = float)    
all_data["Decay_constant_30"] = 1/(all_data['T30']*5).to_numpy(dtype = float)
all_data["active_H"] = all_data['S_initial']*all_data["activity"]/100
all_data["FD_cov"]= np.sqrt(all_data.FD_initial)/all_data.vmax_mean
init_doc_list = np.sort(list(all_data.DOC_initial_int.unique()))
activity_list = np.sort(list(all_data.activity.unique()))
seed_sim_list = np.sort(list(all_data.Seed.unique()))
c_sp_list = list(all_data.carbon_species.unique())
bio_sp_list = list(all_data.biomass_species.unique())
for v in [0.01, 0.1, 0.5, 1., 1.5]:
    for s in seed_sim_list:
        for c in c_sp_list:
            for b in bio_sp_list:
                for d in init_doc_list:
                    sub_scb = all_data[(all_data["Variance"]==v)&(all_data["Seed"].astype(int)==s)&(all_data["carbon_species"].astype(int)==c)&(all_data["biomass_species"].astype(int)==b)&(all_data["DOC_initial_int"].astype(int)==d)]
                    fd_base = sub_scb[sub_scb.Sim_series=="b_1_all_"]['FD_initial'].values[0]
                    fd_cov_base = sub_scb[sub_scb.Sim_series=="b_1_all_"]['FD_cov'].values[0]
                    v_mean_base = sub_scb[sub_scb.Sim_series.isin(['b_1_all_'])]['vmax_mean'].mean()
                    all_data.loc[(all_data["Variance"]==v)&(all_data["Seed"].astype(int)==s)&(all_data["carbon_species"].astype(int)==c)&(all_data["biomass_species"].astype(int)==b)&(all_data["DOC_initial_int"].astype(int)==d), "vmax_mean_base"]=v_mean_base
                    all_data.loc[(all_data["Variance"]==v)&(all_data["Seed"].astype(int)==s)&(all_data["carbon_species"].astype(int)==c)&(all_data["biomass_species"].astype(int)==b)&(all_data["DOC_initial_int"].astype(int)==d), "FD_initial_base"]=fd_base
                    all_data.loc[(all_data["Variance"]==v)&(all_data["Seed"].astype(int)==s)&(all_data["carbon_species"].astype(int)==c)&(all_data["biomass_species"].astype(int)==b)&(all_data["DOC_initial_int"].astype(int)==d), "FD_cov_base"]=fd_cov_base
all_data['vmax_ratio'] = all_data.vmax_mean/all_data.vmax_mean_base
all_data['FD_ratio'] = all_data.FD_maxbio/all_data.FD_initial
all_data['Biomass_ratio'] = all_data.Biomass_maxbio/all_data.Biomass_initial
all_data["FD_initial_ratio"] = all_data.FD_initial/all_data.FD_initial_base                    
all_data["FD_cov_ratio"] = all_data.FD_cov/all_data.FD_cov_base
all_data.to_csv(os.path.join(project_dir,"simulation_results_temporal_decay_const.csv"),index=False)