## Import libraries
import os
import pandas as pd
import numpy as np

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
results_dir = os.path.join(project_dir, "results")

seed_sim_list = [610229235, 983307757, 643338060, 714504443, 277077803, 898393994, 420,13012022,13061989]
cn_list = [3,6,12,18]

count = 0
for c_n in cn_list:
    row = []
    filename = os.path.join(results_dir, 'carbon_' + str(c_n) + "_diversity_data.pkl")
    diversity_data = pd.read_pickle(filename)
    # Additional data processing
    diversity_data['DOC_removal'] = (1 - diversity_data.DOC_end/diversity_data.DOC_input) * 100
    diversity_data['carbon_biomass'] = diversity_data.carbon_species * diversity_data.biomass_species
    diversity_data = diversity_data.replace('NA', np.nan)
    diversity_data['t_50_days']  = diversity_data.T_50
    diversity_data['t_50_b1_days']  = diversity_data.T_50_B1
    if count == 0:
        combined_data = diversity_data
    else:
        combined_data = pd.concat([combined_data, diversity_data]).reindex()
    count +=1

print(combined_data.shape)
combined_data.drop_duplicates()
print(combined_data.shape)

cases = list(diversity_data.Sim_series.unique())
combined_data['activity'] = combined_data.Sim_series
for c, a in zip (cases, [100, 25, 25, 25, 50, 50, 50, 75, 75, 75]):        
    combined_data.loc[combined_data['Sim_series']==c, 'activity'] = a

combined_data = combined_data.drop_duplicates()
combined_data.to_csv(os.path.join(project_dir, "results", "combined_dataset.csv"))
combined_data.to_pickle(os.path.join(project_dir, "results", "combined_dataset.pkl"))