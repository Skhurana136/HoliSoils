## Import libraries
import os
from re import S
import pandas as pd

## LOAD RESULTS
#project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
project_dir = "C:/Users/swkh9804/Documents/Projects/HoliSoils/data"
subfolders = ['carbon_3', 'carbon_6', 'carbon_8','carbon_15','carbon_18']
bio_n_series_master = [[2,3,5,6],[3,6,9,12],[4,6,8,16],[3,5,9,12,15,20,30],[3,9,12,15,18]]

for s, bio_n_series in zip(subfolders, bio_n_series_master):
    details_subfolder = "m_5_"+s + "_ip_0"
    print(details_subfolder)
    simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
    filename = os.path.join(simulations_dir, "seeds_randoms.pkl")
    seed_details = pd.read_pickle(filename)
    row = []
    for base in ["b_1"]:
        #print(base)
        for case in ["all"]:
            #print(case)
            for c_b_r in bio_n_series:
                #print(c_b_r)
                for dom_init in [1000,5000,10000,15000,20000]:
                    #print(dom_init)
                    sim_id = base + "_" + case + "_/bio_n_" + str(c_b_r) + "/dom_initial_" + str(dom_init) + "/seed_13061989"
                    status = seed_details[sim_id]["sim_status"]
                    if status == "failed":
                        print(base, case, c_b_r, dom_init)
    for base in ["b_2", "b_3", "b_4", "b_5"]:
        #print(base)
        for case in ["a", "b", "c", "d", "e"]:
            #print(case)
            for c_b_r in bio_n_series:
                #print(c_b_r)
                for dom_init in [1000,5000,10000,15000,20000]:
                    #print(dom_init)
                    sim_id = base + "_" + case + "_/bio_n_" + str(c_b_r) + "/dom_initial_" + str(dom_init) + "/seed_13061989"
                    status = seed_details[sim_id]["sim_status"]
                    if status == "failed":
                        print(base, case, c_b_r, dom_init)

