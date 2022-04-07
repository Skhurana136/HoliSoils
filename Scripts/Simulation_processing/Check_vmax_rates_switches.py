## Import libraries
import os
import pandas as pd
import sys

old_stdout = sys.stdout

## LOAD RESULTS
#project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
project_dir = "C:/Users/swkh9804/Documents/Projects/HoliSoils/data"
subfolders = ['carbon_3', 'carbon_6','carbon_9','carbon_15']

log_file = open(os.path.join(project_dir, "switches_of_failed_sims_ip_0.log"),"w")

sys.stdout = log_file

for details_subfolder in subfolders:
    print(details_subfolder)
    simulations_dir = os.path.join(project_dir, "simulations", details_subfolder+"_ip_0")
    filename = os.path.join(simulations_dir, "seeds_randoms.pkl")
    seed_details = pd.read_pickle(filename)
    row = []
    for base in ["b_1"]:
        #print(base)
        for case in ["all"]:
            #print(case)
            for c_b_r in [1/3, 1, 3]:
                #print(c_b_r)
                for dom_init in [1000,5000,10000,15000,20000]:
                    #print(dom_init)
                    sim_id = base + "_" + case + "_/c_b_ratio_" + str(c_b_r) + "/dom_initial_" + str(dom_init) + "/seed_13061989"
                    status = seed_details[sim_id]["sim_status"]
                    if status == "failed":
                        print(base, case, c_b_r, dom_init, status)
                        print(seed_details[sim_id]["max_rate_parameters"])
    for base in ["b_2", "b_3", "b_4", "b_5"]:
        #print(base)
        for case in ["a", "b", "c", "d", "e"]:
            #print(case)
            for c_b_r in [1/3, 1, 3]:
                #print(c_b_r)
                for dom_init in [1000,5000,10000,15000,20000]:
                    #print(dom_init)
                    sim_id = base + "_" + case + "_/c_b_ratio_" + str(c_b_r) + "/dom_initial_" + str(dom_init) + "/seed_13061989"
                    status = seed_details[sim_id]["sim_status"]
                    if status == "failed":
                        print(base, case, c_b_r, dom_init, status)
                        print(seed_details[sim_id]["max_rate_parameters"])

sys.stdout = old_stdout

log_file.close()