## Import libraries
import os
import pandas as pd
import numpy as np

import h5py 


## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
results_dir = os.path.join(project_dir, "results")

seed_sim_list = seed_sim_list = [610229235, 983307757, 643338060, 714504443, 277077803, 898393994, 420,13012022,13061989]

cn_list = [3,6,12,18]
bio_n_series = [4,8,12,16,20,24,28,32]
ip = 0
init_dom_list = [1000,2000,5000,10000,15000]
transient_switch = 0
input_factor = transient_switch*5/365
loss_crit = 0.37#0.63

for c_n in cn_list:
    row = []
    results_filename = os.path.join(results_dir, 'adapt_carbon_sw_off_' + str(c_n) + "_c_b_dynamics.pkl")

    for seed_sim in seed_sim_list:
        # Load all datasets and save their Shannon and diversity indices in a dataframe
        seed_all = 'seed_'+str(seed_sim)
        
        details_subfolder = 'adapt_carbon_sw_off_' + str(c_n) + '_'+str(seed_sim) + '_ip_' + str(ip)
        simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
        hr = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'r')

        filename = os.path.join(simulations_dir, "seeds_randoms.pkl")
        seed_details = pd.read_pickle(filename)

        for b_n in bio_n_series:
            for t_dom_initial in init_dom_list:
                base_case = "b_1_all_"
                c_b = "bio_n_"+ str(b_n)
                dom_init = "dom_initial_" + str(t_dom_initial)
                doc_input = (t_dom_initial) * input_factor
                if hr[base_case][c_b][dom_init][seed_all]:
                    sim_data = hr[base_case][c_b][dom_init][seed_all]
                    sim_seed_details = base_case+"/"+c_b+"/"+dom_init+"/"+seed_all
                    init = sim_data['initial_conditions']
                    paras = sim_data['parameters']
                    dom_n, bio_n = sim_data['species_number']
                    x = sim_data['solution']
                    #sim_status = seed_details[sim_seed_details]['sim_status']
                    C = np.asarray(x['dom'])
                    B = np.asarray(x['biomass'])
                    DOC = np.sum(C,axis=1)
                    Biomass = np.sum(B, axis = 1)
                    proportion = B/B.sum(axis=1, keepdims = True)
                    S = -np.sum(proportion*np.log(proportion), axis = 1)
                    S_i = S[0]
                    C_B_r = DOC/Biomass
                    C_B_r_i = C_B_r[0]
                    C_B_r_end = C_B_r[-1] 
                    C_B_r_min = np.min(C_B_r)
                    C_B_r_max = np.max(C_B_r)
                    C_B_r_median = np.median(C_B_r)
                    C_B_r_mean = np.mean(C_B_r)
                    C_i_max_o_s = seed_details[sim_seed_details]['oxidation_state'][np.argmax(C[0,:])]
                    C_end_max_o_s = seed_details[sim_seed_details]['oxidation_state'][np.argmax(C[-1,:])]
                    C_i_max_o_s_pc = np.max(C[0,:])/DOC[0]
                    C_end_max_o_s_pc = np.max(C[-1,:])/DOC[-1]
                    row.append([seed_sim,base_case, dom_n, bio_n, t_dom_initial, S_i, C_B_r_min, C_B_r_max, C_B_r_median, C_B_r_mean, C_B_r_i, C_B_r_end, C_i_max_o_s, C_end_max_o_s, C_i_max_o_s_pc, C_end_max_o_s_pc])
                for baseline in ["b_2", "b_3", "b_4","b_5"]:
                    for label in ["a", "b", "c","d","e"]:
                        sim = baseline + "_" + label + "_"
                        sim_seed_details = sim+"/"+c_b+"/"+dom_init+"/"+seed_all
                        if hr[sim][c_b][dom_init][seed_all]:
                            sim_data = hr[sim][c_b][dom_init][seed_all]
                            init = sim_data['initial_conditions']
                            paras = sim_data['parameters']
                            dom_n, bio_n = sim_data['species_number']
                            x = sim_data['solution']
                            #sim_status = seed_details[sim_seed_details]['sim_status']
                            C = np.asarray(x['dom'])
                            B = np.asarray(x['biomass'])
                            DOC = np.sum(C,axis=1)
                            Biomass = np.sum(B, axis = 1)
                            proportion = B/B.sum(axis=1, keepdims = True)
                            S = -np.sum(proportion*np.log(proportion), axis = 1)
                            S_i = S[0]
                            C_B_r = DOC/Biomass
                            C_B_r_i = C_B_r[0]
                            C_B_r_end = C_B_r[-1] 
                            C_B_r_min = np.min(C_B_r)
                            C_B_r_max = np.max(C_B_r)
                            C_B_r_median = np.median(C_B_r)
                            C_B_r_mean = np.mean(C_B_r)
                            C_i_max_o_s = seed_details[sim_seed_details]['oxidation_state'][np.argmax(C[0,:])]
                            C_end_max_o_s = seed_details[sim_seed_details]['oxidation_state'][np.argmax(C[-1,:])]
                            C_i_max_o_s_pc = np.max(C[0,:])/DOC[0]
                            C_end_max_o_s_pc = np.max(C[-1,:])/DOC[-1]
                            row.append([seed_sim, sim, dom_n, bio_n, t_dom_initial, S_i, C_B_r_min, C_B_r_max, C_B_r_median, C_B_r_mean, C_B_r_i, C_B_r_end, C_i_max_o_s, C_end_max_o_s, C_i_max_o_s_pc, C_end_max_o_s_pc])
        hr.close()

    c_b_dynamics = pd.DataFrame.from_records(row, columns = ["Seed", "Sim_series", "carbon_species", "biomass_species", "DOC_initial", "S_initial", "C_B_min", "C_B_max", "C_B_median", "C_B_mean", "C_B_initial", "C_B_end", "DOC_initial_max_os", "DOC_end_max_os", "DOC_initial_max_os_pc", "DOC_end_max_os_pc"])

    print("The shape of the dataframe is ", c_b_dynamics.shape)
    print("The dataframe contains the following data types ", c_b_dynamics.dtypes)

    filename = os.path.join(results_dir, "c_b_dynamics.pkl")
    c_b_dynamics.to_pickle(results_filename)
    print ("Diversity with carbon data is saved here ", results_filename)