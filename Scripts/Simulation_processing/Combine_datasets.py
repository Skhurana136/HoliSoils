## Import libraries
import os
import pandas as pd
import numpy as np

## LOAD RESULTS
#project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
project_dir = "C:/Users/swkh9804/Documents/Projects/HoliSoils/data"
for case in ["carbon_3","carbon_6", "carbon_8", "carbon_15", "carbon_18"]:
    details_subfolder = case + '_seed_420_ip_0'
    simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
    results_dir = os.path.join(project_dir, "results", details_subfolder)
    figures_dir = os.path.join(project_dir, "figures", details_subfolder)
    filename = os.path.join(results_dir, "diversity_data.pkl")
    diversity_data = pd.read_pickle(filename)
    # Additional data processing
    diversity_data['DOC_removal'] = (1 - diversity_data.DOC_end/diversity_data.DOC_input) * 100
    diversity_data['carbon_biomass'] = diversity_data.carbon_species * diversity_data.biomass_species
    diversity_data = diversity_data.replace('NA', np.nan)
    diversity_data['t_50_days']  = diversity_data.T_50/5
    diversity_data['t_50_b1_days']  = diversity_data.T_50_B1/5
    print(diversity_data.shape)
    print(diversity_data.columns)
    if case == "carbon_3":
        combined_data = diversity_data
    else:
        combined_data = pd.concat([combined_data, diversity_data]).reindex()

print(combined_data.shape)

cases = list(diversity_data.Sim_series.unique())
combined_data['activity'] = combined_data.Sim_series
for c, a in zip (cases, [100, 10, 10, 10, 10, 10, 30, 30, 30, 30, 30, 50, 50, 50, 50, 50, 70, 70, 70, 70, 70]):        
    combined_data.loc[combined_data['Sim_series']==c, 'activity'] = a

combined_data.to_csv(os.path.join(project_dir, "results", "combined_dataset_420.csv"))
combined_data.to_pickle(os.path.join(project_dir, "results", "combined_dataset_420.pkl"))