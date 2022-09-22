## Import libraries
import os
import numpy as np
import h5py
import sys
import pandas as pd

def gather_paras(sim_data):
    dom_n, bio_n = sim_data['species_number']
    init_bio = sim_data['initial_conditions']['biomass']
    init_bio_sum = np.sum(init_bio)
    init_c = sim_data['initial_conditions']['dom']
    init_c_sum = np.sum(init_c)
    paras = sim_data['parameters']
    dying = np.tile(init_bio*np.asarray(paras['mortality_const'])/init_bio_sum, (dom_n, 1)).flatten().reshape(-1,1)
    exo_enz = np.tile(init_bio*np.asarray(paras['exo_enzyme_rate'])/init_bio_sum, (dom_n, 1)).flatten().reshape(-1,1)
    v_params = (init_bio*np.asarray((paras['max_rate']))/init_bio_sum).flatten().reshape(-1,1)*np.asarray(paras['carbon_uptake']).flatten().reshape(-1,1)*exo_enz
    k_params = np.asarray(paras['half_saturation']).flatten().reshape(-1,1)
    #ox = np.asarray(list([x*init_c[idx]/init_c_sum]*bio_n for x,idx in enumerate(paras['oxidation_state']))).flatten().reshape(-1,1)
    ox = np.asarray(list([x]*bio_n for x in paras['oxidation_state'])).flatten().reshape(-1,1)
    paras_arr = np.append(np.append(np.append(dying, exo_enz, axis=1),np.append(v_params, k_params, axis = 1), axis = 1), ox, axis = 1)
    return paras_arr

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient", "activity_loss_-02")
results_dir = os.path.join(project_dir, "results")
filestring = sys.argv[1] + "_carbon_"
seed_sim_list = [610229235, 983307757, 643338060, 714504443, 277077803, 898393994, 420,13012022,13061989]
cn_list = [3,6,12,18]
bio_n_series = [4,8,16,32]
init_dom_list = [1000,2000,5000,10000,15000]

def run_gather_params(base_list, base_list_act):
    init_dom_list = [1000,2000,5000,10000,15000]
    all_para_arr = np.zeros((0,10))
    for c_n in cn_list:
        row = []
        c_b_row = []
        #results_filename = os.path.join(results_dir, filestring + str(c_n))

        for seed_sim in seed_sim_list:
            # Load all datasets and save their Shannon and diversity indices in a dataframe
            seed_all = 'seed_'+str(seed_sim)
            
            for b_n in bio_n_series:
                details_subfolder = filestring + str(c_n) + '_'+str(seed_sim) + '_ip_0'
                simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
                hrw = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'r+')
                for t_dom_initial in init_dom_list:
                    c_b = "bio_n_"+ str(b_n)
                    dom_init = "dom_initial_" + str(t_dom_initial)
                    doc_input = (t_dom_initial)
                    seed_arr = np.zeros((c_n*b_n,1))+seed_sim
                    bio_arr = np.zeros((c_n*b_n,1))+b_n
                    seed_arr_bio = np.append(seed_arr, bio_arr, axis = 1) 
                    carbon_arr = np.zeros((c_n*b_n,1))+c_n
                    seed_carbon = np.append(seed_arr_bio, carbon_arr, axis = 1)
                    for base, act in zip(base_list, base_list_act):#zip(["b_1", "b_2", "b_3", "b_4", "b_5"], [1.0, 0.0, 0.25, 0.5, 0.75]):
                        act_arr = np.zeros((c_n*b_n,1))+act
                        seed_sim_arr = np.append(seed_carbon, act_arr, axis = 1)
                        if base == "b_1":
                            base_case = "b_1_all_"
                            level_arr = np.zeros((c_n*b_n,1))+0
                            sim_data = hrw[base_case][c_b][dom_init][seed_all]
                            seed_sim_level = np.append(seed_sim_arr, level_arr, axis = 1)
                            paras_sim = np.append(seed_sim_level, gather_paras(sim_data), axis = 1)
                            all_para_arr = np.append(all_para_arr, paras_sim, axis = 0)
                        else:
                            for label,label_id in zip(["a", "b", "c","d","e"], [1,2,3,4,5]):
                                sim = base + "_" + label + "_"
                                level_arr = np.zeros((c_n*b_n,1))+label_id
                                sim_data = hrw[sim][c_b][dom_init][seed_all]
                                seed_sim_level = np.append(seed_sim_arr, level_arr, axis = 1)
                                paras_sim = np.append(seed_sim_level, gather_paras(sim_data), axis = 1)
                                all_para_arr = np.append(all_para_arr, paras_sim, axis = 0)
            hrw.close()
    para_df = pd.DataFrame(all_para_arr, columns = ['Seed', 'biomass_species','carbon_species','Activity', 'level_id', 'm', 'exo_enz','vmax', 'k', 'oxidation_state'])
    para_df.loc[para_df["level_id"] == 0.0, "Sim_series"] = 'b_1_all_'
    para_df.loc[para_df["level_id"] == 1.0, "Sim_series"] = base_list[0]+'_a_'
    para_df.loc[para_df["level_id"] == 2.0, "Sim_series"] = base_list[0]+'_b_'
    para_df.loc[para_df["level_id"] == 3.0, "Sim_series"] = base_list[0]+'_c_'
    para_df.loc[para_df["level_id"] == 4.0, "Sim_series"] = base_list[0]+'_d_'
    para_df.loc[para_df["level_id"] == 5.0, "Sim_series"] = base_list[0]+'_e_'
    return para_df
    
base_para_df = run_gather_params(["b_1"], [1.0])
on1_para_df = run_gather_params(["b_2"], [0.1])
on2_para_df = run_gather_params(["b_3"], [0.25])
on3_para_df = run_gather_params(["b_4"], [0.5])
on4_para_df = run_gather_params(["b_5"], [0.75])
para_df = pd.concat([base_para_df,on1_para_df, on2_para_df, on3_para_df, on4_para_df])
func_div = para_df.groupby(['Seed', 'biomass_species', 'carbon_species', 'Sim_series', 'level_id'], as_index = False).agg({
    'vmax': ['mean','median','std','skew'],
    'k':['mean','median','std','skew'],
    'm':['mean','median','std','skew'],
    'exo_enz':['mean','median','std','skew']
}).reset_index()
func_div.columns = func_div.columns.map('_'.join)
func_div.rename(columns = {'Seed_':'Seed', 'carbon_species_':'carbon_species', 'biomass_species_': 'biomass_species', 'Sim_series_': 'Sim_series', 'level_id_':'level_id'}, inplace = True)

func_div.to_csv(os.path.join(results_dir, sys.argv[1] + "_parameters.csv"))
