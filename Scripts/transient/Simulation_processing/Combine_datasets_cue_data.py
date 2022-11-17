## Import libraries
import os
import pandas as pd
import numpy as np
import sys

## LOAD RESULTS
#project_dir = os.path.join("C:/", "Users", "swkh9804", "Documents", "Project_data", "HoliSoils", "transient", sys.argv[1])
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient", sys.argv[1])
results_dir = os.path.join(project_dir, "results")
filestring = "competition_adaptation_carbon_"
filestring2 = "_loss_0.9"# + sys.argv[2]
results_filename = filestring + filestring2 + "_cue_combined_dataset"

seed_sim_list = [610229235, 983307757, 643338060, 714504443, 277077803, 898393994, 420,13012022,13061989]
cn_list = [3,6,12,18]

count = 0
for c_n in cn_list:
    row = []
    filename = os.path.join(results_dir, filestring + str(c_n) + filestring2 + "_cue_data.pkl")
    diversity_data = pd.read_pickle(filename)
    # Additional data processing
    diversity_data = diversity_data.replace('NA', np.nan)
    if count == 0:
        combined_data = diversity_data
    else:
        combined_data = pd.concat([combined_data, diversity_data]).reindex()
    count +=1

print(combined_data.shape)
combined_data.drop_duplicates()
print(combined_data.shape)

cases = list(diversity_data.Sim_series.unique())
for c, a in zip (cases, [100, 10, 10, 10, 10, 10, 25, 25, 25, 25, 25, 50, 50, 50, 50, 50, 75, 75, 75, 75, 75]):        
    combined_data.loc[combined_data['Sim_series']==c, 'activity'] = a

combined_data = combined_data.drop_duplicates()
combined_data.to_csv(os.path.join(project_dir, "results", results_filename+".csv"))
combined_data.to_pickle(os.path.join(project_dir, "results", results_filename + ".pkl"))