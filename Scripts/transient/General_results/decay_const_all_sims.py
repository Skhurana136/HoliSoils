#%%
import os
import pandas as pd
import numpy as np

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
sim_suffixes = ["_0_01", "_0_1", "_0_5", "", "_1_5x"]
sim_suffixes_var = [0.01, 0.1, 0.5, 1, 1.5]

files = []
for p,a in zip(sim_suffixes, sim_suffixes_var):
    filename = os.path.join(project_dir, "gen_spec_lognorm" + p, "results", "competition_adaptation_carbon_cue_wremaining_c_data.pkl")
    data = pd.read_pickle(filename)
    var_arr = np.zeros((data.shape[0],))+a
    var_ser = pd.Series(var_arr, copy=False,name = "Variance")
    cases = list(data.Sim_series.unique())
    for c, a in zip (cases, [100, 10, 10, 10, 10, 10, 25, 25, 25, 25, 25, 50, 50, 50, 50, 50, 75, 75, 75, 75, 75]):        
        data.loc[data['Sim_series']==c, 'activity'] = a
    data_var = data.join(var_ser)
    files.append(data_var)
fd_data = pd.concat(files)
print(fd_data.columns)
fd_data['DOC_initial_int'] = round(fd_data.DOC_initial, -3)

### CHARACTERISTIC REACTION TIME SCALE
files = []
for p,a in zip(sim_suffixes, sim_suffixes_var):
    filename = os.path.join(project_dir, "gen_spec_lognorm" + p, "results", "competition_adaptation_carbon__decay_const_c_pools_data_initial_conditions.pkl")
    data = pd.read_pickle(filename)
    var_arr = np.zeros((data.shape[0],))+a
    var_ser = pd.Series(var_arr, copy=False,name = "Variance")
    cases = list(data.Sim_series.unique())
    for c, a in zip (cases, [100, 10, 10, 10, 10, 10, 25, 25, 25, 25, 25, 50, 50, 50, 50, 50, 75, 75, 75, 75, 75]):        
        data.loc[data['Sim_series']==c, 'activity'] = a
    data_var = data.join(var_ser)
    files.append(data_var)
tim_data = pd.concat(files)
print(tim_data.columns)
tim_data['DOC_initial_int'] = round(tim_data.DOC_initial, -3)
fd_tim_data = pd.merge(fd_data, tim_data, on = ["Seed", "Variance", "DOC_initial_int","biomass_species", "carbon_species", "Sim_series", "activity"])
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
#all_data["active_H_c_connections"] = all_data['S_initial']*all_data['carbon_species']*all_data["activity"]/100
#all_data["active_H"] = all_data['S_initial']*all_data["activity"]/100
#all_data["FD_cov"]= np.sqrt(all_data.FD_initial)/all_data.vmax_mean
init_doc_list = np.sort(list(all_data.DOC_initial_int.unique()))
activity_list = np.sort(list(all_data.activity.unique()))
seed_sim_list = np.sort(list(all_data.Seed.unique()))
c_sp_list = list(all_data.carbon_species.unique())
bio_sp_list = list(all_data.biomass_species.unique())
#T_cols = list(x for x in all_data.columns if 'T' in x)
#old_t = 'T0'
#for t_col in T_cols:
#    if t_col == 'T0':
#        all_data['dec0'] = 1/(all_data['T0']*5)
#    else:
#        #all_data['dec'+t_col[1:]] = 1/(all_data[t_col]*5) #execute when all calcs are for subsequent calcs
#        all_data['T_diff'+t_col[1:]] = abs((all_data[t_col]-all_data[old_t])) #execute when all calcs are wrt initial conditions
#        all_data['dec'+t_col[1:]] = 1/(all_data['T_diff'+t_col[1:]]*5) #execute when all calcs are wrt initial conditions
#    all_data['dec_ratio'+t_col[1:]]=all_data['dec'+t_col[1:]]/all_data['dec0']#compl = all_data.dropna(subset = ['ratio_t_50'])
#   old_t = t_col

#for v in [0.01, 0.1, 0.5, 1., 1.5]:
#    for s in seed_sim_list:
#        for c in c_sp_list:
#            for b in bio_sp_list:
#                for d in init_doc_list:
#                    sub_scb = all_data[(all_data["Variance"]==v)&(all_data["Seed"].astype(int)==s)&(all_data["carbon_species"].astype(int)==c)&(all_data["biomass_species"].astype(int)==b)&(all_data["DOC_initial_int"].astype(int)==d)]
#                    fd_base = sub_scb[sub_scb.Sim_series=="b_1_all_"]['FD_initial'].values[0]
#                    fd_cov_base = sub_scb[sub_scb.Sim_series=="b_1_all_"]['FD_cov'].values[0]
#                    v_mean_base = sub_scb[sub_scb.Sim_series.isin(['b_1_all_'])]['vmax_mean'].mean()
#                    all_data.loc[(all_data["Variance"]==v)&(all_data["Seed"].astype(int)==s)&(all_data["carbon_species"].astype(int)==c)&(all_data["biomass_species"].astype(int)==b)&(all_data["DOC_initial_int"].astype(int)==d), "vmax_mean_base"]=v_mean_base
#                    all_data.loc[(all_data["Variance"]==v)&(all_data["Seed"].astype(int)==s)&(all_data["carbon_species"].astype(int)==c)&(all_data["biomass_species"].astype(int)==b)&(all_data["DOC_initial_int"].astype(int)==d), "FD_initial_base"]=fd_base
#                    all_data.loc[(all_data["Variance"]==v)&(all_data["Seed"].astype(int)==s)&(all_data["carbon_species"].astype(int)==c)&(all_data["biomass_species"].astype(int)==b)&(all_data["DOC_initial_int"].astype(int)==d), "FD_cov_base"]=fd_cov_base
#all_data['vmax_ratio'] = all_data.vmax_mean/all_data.vmax_mean_base
#all_data['FD_ratio'] = all_data.FD_maxbio/all_data.FD_initial
#all_data['Biomass_ratio'] = all_data.Biomass_maxbio/all_data.Biomass_initial
#all_data["FD_initial_ratio"] = all_data.FD_initial/all_data.FD_initial_base                    
#all_data["FD_cov_ratio"] = all_data.FD_cov/all_data.FD_cov_base
all_data.to_csv(os.path.join(project_dir,"simulation_results_temporal_initial_conditions_decay_const.csv"),index=False)