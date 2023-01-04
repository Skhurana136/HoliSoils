## Import libraries
import os
import pandas as pd
import numpy as np

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
results_dir = os.path.join(project_dir, "results")
filestring = "competition_adaptation_carbon__decay_const_c_pools_data_initial_conditions"
results_filename = filestring + "_combined_dataset"
sim_suffixes = ["_0_01", "_0_1", "_0_5", "", "_1_5x"]
sim_suffixes_var = [0.01, 0.1, 0.5, 1, 1.5]

### CHARACTERISTIC REACTION TIME SCALE
files = []
for p,a in zip(sim_suffixes, sim_suffixes_var):
    filename = os.path.join(project_dir, "gen_spec_lognorm" + p, "results", filestring+".pkl")
    data = pd.read_pickle(filename)
    data['decay_const'] = 1/(data['T_loss']*5)
    data['DOC_initial_int'] = round(data.DOC_initial, -3)
    initial_data = data[data['%C']==100.]
    initial_data = initial_data[["Seed", "C_pool","DOC_initial_int","biomass_species", "carbon_species", "Sim_series","decay_const"]]
    initial_data.rename(columns={'decay_const':'decay_const_initial'}, inplace=True)
    print(initial_data.columns)
    data = pd.merge(data,initial_data, on = ["Seed", "C_pool", "DOC_initial_int","biomass_species", "carbon_species", "Sim_series"], suffixes=('', '_y'))
    data.drop(data.filter(regex='_y$').columns, axis=1, inplace=True)
    data['decay_ratio'] = data.decay_const/data.decay_const_initial
    var_arr = np.zeros((data.shape[0],))+a
    var_ser = pd.Series(var_arr, copy=False,name = "Variance")
    cases = list(data.Sim_series.unique())
    for c, a in zip (cases, [100, 10, 10, 10, 10, 10, 25, 25, 25, 25, 25, 50, 50, 50, 50, 50, 75, 75, 75, 75, 75]):        
        data.loc[data['Sim_series']==c, 'activity'] = a
    data = data.join(var_ser)
    data.replace({'decay_ratio': {np.inf:-1, -np.inf:-1, np.nan:-1}}, inplace=True)
    maxid = data.groupby(["Seed", "C_pool","DOC_initial_int","biomass_species", "carbon_species", "Sim_series"])['decay_ratio'].transform('idxmax').values
    minid = data.groupby(["Seed", "C_pool","DOC_initial_int","biomass_species", "carbon_species", "Sim_series"])['decay_ratio'].transform('idxmin').values
    data['decay_max'] = data.loc[maxid,'%C'].values
    data['decay_stop'] = data.loc[minid,'%C'].values
    files.append(data)

tim_data = pd.concat(files)

print(tim_data.columns)
print(tim_data.shape)
tim_data.drop_duplicates()
tim_data.to_csv(os.path.join(project_dir, "results", results_filename+".csv"))
tim_data.to_pickle(os.path.join(project_dir, "results", results_filename + ".pkl"))