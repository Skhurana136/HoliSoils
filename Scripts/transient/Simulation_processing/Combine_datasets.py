## Import libraries
import os
import pandas as pd
import numpy as np

## LOAD RESULTS
#project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
project_dir = "D:\Projects\HoliSoils\data"
results_dir = os.path.join(project_dir, "results")

seed_sim_list = [420,13012022,13061989]
ip_list = [0]#,1,2]
carbon_nums = [3,6,8,15,18]

count = 0
for c_case in carbon_nums:
    case = 'carbon_' + str(c_case)
    for seed_sim in seed_sim_list:
        for ip in ip_list:
            details_subfolder = case + '_' + str(seed_sim) + '_ip_' + str(ip)
            filename = os.path.join(results_dir, "carbon_"+str(c_case)+"_ip_"+str(ip)+"_diversity_data.pkl")
            diversity_data = pd.read_pickle(filename)
            # Additional data processing
            diversity_data['DOC_removal'] = (1 - diversity_data.DOC_end/diversity_data.DOC_input) * 100
            diversity_data['carbon_biomass'] = diversity_data.carbon_species * diversity_data.biomass_species
            diversity_data = diversity_data.replace('NA', np.nan)
            diversity_data['t_50_days']  = diversity_data.T_50/5
            diversity_data['t_50_b1_days']  = diversity_data.T_50_B1/5
            print(diversity_data.shape)
            print(diversity_data.columns)
            if count == 0:
                combined_data = diversity_data
            else:
                combined_data = pd.concat([combined_data, diversity_data]).reindex()
            count +=1

print(combined_data.shape)

cases = list(diversity_data.Sim_series.unique())
combined_data['activity'] = combined_data.Sim_series
for c, a in zip (cases, [100, 10, 10, 10, 10, 10, 30, 30, 30, 30, 30, 50, 50, 50, 50, 50, 70, 70, 70, 70, 70]):        
    combined_data.loc[combined_data['Sim_series']==c, 'activity'] = a

combined_data.to_csv(os.path.join(project_dir, "results", "combined_dataset.csv"))
combined_data.to_pickle(os.path.join(project_dir, "results", "combined_dataset.pkl"))