## Import libraries
import os
import pandas as pd
import numpy as np
import sys
import itertools

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient", sys.argv[1])
results_dir = os.path.join(project_dir, "results")
filestring = "competition_adaptation_carbon_"
results_filename = filestring + "_decay_const_c_pools_data_initial_conditions_combined_dataset"
filename = os.path.join(results_dir, filestring + "_decay_const_c_pools_data_initial_conditions.pkl")
all_data = pd.read_pickle(filename)
all_data['DOC_initial']=all_data["DOC_initial"].astype(int)    
seed_sim_list = all_data.Seed.unique().tolist()#[610229235, 983307757, 643338060, 714504443, 277077803, 898393994, 420,13012022,13061989]
cn_list = all_data.carbon_species.unique().tolist()
bios = all_data.biomass_species.unique().tolist()
carb_types = ['DOC','reduced_C','necromass','oxidized_C']
docs = all_data.DOC_initial.unique().tolist()
sims = all_data.Sim_series.unique().tolist()
all_data.rename(columns={'DOC':'Tdoc', 'reduced_C':'Tredc', 'oxidized_C':'Toxc'}, inplace=True)
count = 0
for c, b, d, seed, s in itertools.product(cn_list,bios, docs, seed_sim_list, sims):
    subset = all_data[(all_data["Seed"].astype(int)==seed)&(all_data["carbon_species"].astype(int)==c)&(all_data["biomass_species"].astype(int)==b)&(all_data["DOC_initial"].astype(int)==d)&(all_data["Sim_series"]==s)]
    print(c,b,d,seed,s)
    print(subset[['%C','Tdoc']])
    #print(c_n, b_n, doci, seed, sim)
#for c_n in cn_list:
#    row = []
#    # Additional data processing
#    if count == 0:
#        combined_data = diversity_data
#    else:
#        combined_data = pd.concat([combined_data, diversity_data]).reindex()
#    count +=1

#print(combined_data.shape)
#combined_data.drop_duplicates()
#print(combined_data.shape)

#cases = list(diversity_data.Sim_series.unique())
#for c, a in zip (cases, [100, 10, 10, 10, 10, 10, 25, 25, 25, 25, 25, 50, 50, 50, 50, 50, 75, 75, 75, 75, 75]):        
#    combined_data.loc[combined_data['Sim_series']==c, 'activity'] = a

#combined_data = combined_data.drop_duplicates()
#combined_data.to_csv(os.path.join(project_dir, "results", results_filename+".csv"))
#combined_data.to_pickle(os.path.join(project_dir, "results", results_filename + ".pkl"))