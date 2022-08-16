## Import libraries
import os
import numpy as np
import h5py
import sys
import pandas as pd

def gather_paras(sim_data):
    dom_n, bio_n = sim_data['species_number']
    paras = sim_data['parameters']
    v_params = np.asarray(paras['max_rate']).flatten().reshape(-1,1)
    k_params = np.asarray(paras['half_saturation']).flatten().reshape(-1,1)
    ox = np.asarray(list([x]*bio_n for x in paras['oxidation_state'])).flatten().reshape(-1,1)
    paras_arr = np.append(np.append(v_params, k_params, axis = 1), ox, axis = 1)
    return paras_arr

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient", "activity_loss_-02")
results_dir = os.path.join(project_dir, "results")
filestring = sys.argv[1] + "_carbon_"
seed_sim_list = [610229235, 983307757, 643338060, 714504443, 277077803, 898393994, 420,13012022,13061989]
cn_list = [3,6,12,18]
bio_n_series = [4,8,16,32]
init_dom_list = [1000,2000,5000,10000,15000]

all_para_arr = np.zeros((0,4))
for c_n in cn_list:
    row = []
    c_b_row = []
    results_filename = os.path.join(results_dir, filestring + str(c_n))

    for seed_sim in seed_sim_list:
        # Load all datasets and save their Shannon and diversity indices in a dataframe
        seed_all = 'seed_'+str(seed_sim)
        
        details_subfolder = filestring + str(c_n) + '_'+str(seed_sim) + '_ip_0'
        simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
        hrw = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'r+')

        for b_n in bio_n_series:
            for t_dom_initial in init_dom_list:
                c_b = "bio_n_"+ str(b_n)
                dom_init = "dom_initial_" + str(t_dom_initial)
                doc_input = (t_dom_initial)
                for base, act in zip(["b_1","b_2", "b_3", "b_4", "b_5"], [1, 0.1, 0.25, 0.5, 0.75]):
                    if base == "b_1":
                        base_case = "b_1_all_"
                        sim_data = hrw[base_case][c_b][dom_init][seed_all]
                        paras_sim = np.append(gather_paras(sim_data), np.zeros((c_n*b_n,1))+act, axis = 1)
                        all_para_arr = np.append(all_para_arr, paras_sim, axis = 0)
                    else:
                        for label in ["a", "b", "c","d","e"]:
                            sim = base + "_" + label + "_"
                            sim_data = hrw[sim][c_b][dom_init][seed_all]
                            paras_sim = np.append(gather_paras(sim_data), np.zeros((c_n*b_n,1))+act, axis = 1)
                            all_para_arr = np.append(all_para_arr, paras_sim, axis = 0)
        hrw.close()

import pandas as pd
df = pd.DataFrame(all_para_arr, columns = ['vmax','Ks','Oxidation_state','Activity'])
df.to_csv(results_dir, sys.argv[1] + "_parameters.csv")
