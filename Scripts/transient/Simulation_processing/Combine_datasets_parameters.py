## Import libraries
import os
import pandas as pd
import numpy as np

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
results_dir = os.path.join(project_dir, "results")
results_filename = "competition_adaptation_parameters_combined_dataset"
sim_suffixes = ["_0_01", "_0_1", "_0_5", "", "_1_5x"]
sim_suffixes_var = [0.01, 0.1, 0.5, 1, 1.5]

### CHARACTERISTIC REACTION TIME SCALE
files = []
for p,a in zip(sim_suffixes, sim_suffixes_var):
    filename = os.path.join(project_dir, "gen_spec_lognorm" + p, "results", "gen_spec_lognorm" + p + "_parameters.csv")
    data = pd.read_csv(filename)
    var_arr = np.zeros((data.shape[0],))+a
    var_ser = pd.Series(var_arr, copy=False,name = "Variance")
    cases = list(data.Sim_series.unique())
    for c, a in zip (cases, [100, 10, 10, 10, 10, 10, 25, 25, 25, 25, 25, 50, 50, 50, 50, 50, 75, 75, 75, 75, 75]):        
        data.loc[data['Sim_series']==c, 'activity'] = a
    data_var = data.join(var_ser)
    files.append(data_var)
para_data = pd.concat(files)

tim_data = pd.concat(files)
tim_data['DOC_initial_int'] = round(tim_data.DOC_initial, -3)
tim_data = pd.concat(files)
print(tim_data.columns)
print(tim_data.shape)
tim_data.drop_duplicates()
tim_data.to_csv(os.path.join(project_dir, "results", results_filename+".csv"))
tim_data.to_pickle(os.path.join(project_dir, "results", results_filename + ".pkl"))