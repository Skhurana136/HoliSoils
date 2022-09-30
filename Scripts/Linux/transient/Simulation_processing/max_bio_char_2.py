## Import libraries
import os
import pandas as pd
import numpy as np
import sys
import h5py 


## LOAD RESULTS
project_dir = os.path.join('/proj', 'hs_micro_div_072022', 'Project_data', 'transient', 'gen_spec_lognorm_1_5x')
results_dir = os.path.join(project_dir, "results")
filestring = sys.argv[1] + '_carbon_' #null
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
    results_filename = os.path.join(results_dir, filestring + str(c_n))

    for seed_sim in seed_sim_list:
        # Load all datasets and save their Shannon and diversity indices in a dataframe
        seed_all = 'seed_'+str(seed_sim)
        
        details_subfolder = filestring + str(c_n) + '_'+str(seed_sim) + '_ip_' + str(ip)
        simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
        hr = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'r')

        for b_n in bio_n_series:
            for t_dom_initial in init_dom_list:
                base_case = "b_1_all_"
                c_b = "bio_n_"+ str(b_n)
                dom_init = "dom_initial_" + str(t_dom_initial)
                doc_input = (t_dom_initial) * input_factor
                if hr[base_case][c_b][dom_init][seed_all]:
                    sim_data = hr[base_case][c_b][dom_init][seed_all]
                    paras = sim_data['parameters']
                    dom_n, bio_n = sim_data['species_number']
                    x = sim_data['solution']
                    C = np.asarray(x['dom'])
                    B = np.asarray(x['biomass'])
                    B_max_initial = np.max(B[0,:])
                    B_max_initial_which = np.argmax(B[0,:])
                    k_max_bmax_initial =  np.max(paras["half_saturation"][:,B_max_initial_which])
                    vmax_max_bmax_initial = np.max(paras["max_rate"][:,B_max_initial_which])
                    k_mean_bmax_initial = np.mean(paras["half_saturation"][:,B_max_initial_which])
                    vmax_mean_bmax_initial = np.mean(paras["max_rate"][:,B_max_initial_which])
                    B_max_all_time = np.max(B)
                    B_max_all_time_when = np.where(B==B_max_all_time)[0]
                    B_max_all_time_which = np.where(B==B_max_all_time)[1]
                    k_max_bmax_all_time =  np.max(paras["half_saturation"][:,B_max_all_time_which])
                    vmax_max_bmax_all_time = np.max(paras["max_rate"][:,B_max_all_time_which])
                    k_mean_bmax_all_time = np.mean(paras["half_saturation"][:,B_max_all_time_which])
                    vmax_mean_bmax_all_time = np.mean(paras["max_rate"][:,B_max_all_time_which])
                    B_max_end = np.max(B[-1,:])
                    B_max_end_which = np.argmax(B[-1,:])
                    k_max_bmax_end =  np.max(paras["half_saturation"][:,B_max_end_which])
                    vmax_max_bmax_end = np.max(paras["max_rate"][:,B_max_end_which])
                    k_mean_bmax_end = np.mean(paras["half_saturation"][:,B_max_end_which])
                    vmax_mean_bmax_end = np.mean(paras["max_rate"][:,B_max_end_which])
                    row.append([seed_sim,base_case, dom_n, bio_n, t_dom_initial, B_max_initial, B_max_initial_which, k_max_bmax_initial, k_mean_bmax_initial, vmax_max_bmax_initial, vmax_mean_bmax_initial, B_max_all_time, B_max_all_time_which, B_max_all_time_when, k_max_bmax_all_time, k_mean_bmax_all_time, vmax_max_bmax_all_time, vmax_mean_bmax_all_time, B_max_end, B_max_end_which, k_max_bmax_end, k_mean_bmax_end, vmax_max_bmax_end, vmax_mean_bmax_end])
                for baseline in ["b_2", "b_3", "b_4","b_5"]:
                    for label in ["a", "b", "c","d","e"]:
                        sim = baseline + "_" + label + "_"
                        if hr[sim][c_b][dom_init][seed_all]:
                            sim_data = hr[sim][c_b][dom_init][seed_all]
                            init = sim_data['initial_conditions']
                            paras = sim_data['parameters']
                            dom_n, bio_n = sim_data['species_number']
                            x = sim_data['solution']
                            C = np.asarray(x['dom'])
                            B = np.asarray(x['biomass'])
                            B_max_initial = np.max(B[0,:])
                            B_max_initial_which = np.argmax(B[0,:])
                            k_max_bmax_initial =  np.max(paras["half_saturation"][:,B_max_initial_which])
                            vmax_max_bmax_initial = np.max(paras["max_rate"][:,B_max_initial_which])
                            k_mean_bmax_initial = np.mean(paras["half_saturation"][:,B_max_initial_which])
                            vmax_mean_bmax_initial = np.mean(paras["max_rate"][:,B_max_initial_which])
                            B_max_all_time = np.max(B)
                            B_max_all_time_when = np.where(B==B_max_all_time)[0]
                            B_max_all_time_which = np.where(B==B_max_all_time)[1]
                            k_max_bmax_all_time =  np.max(paras["half_saturation"][:,B_max_all_time_which])
                            vmax_max_bmax_all_time = np.max(paras["max_rate"][:,B_max_all_time_which])
                            k_mean_bmax_all_time = np.mean(paras["half_saturation"][:,B_max_all_time_which])
                            vmax_mean_bmax_all_time = np.mean(paras["max_rate"][:,B_max_all_time_which])
                            B_max_end = np.max(B[-1,:])
                            B_max_end_which = np.argmax(B[-1,:])
                            k_max_bmax_end =  np.max(paras["half_saturation"][:,B_max_end_which])
                            vmax_max_bmax_end = np.max(paras["max_rate"][:,B_max_end_which])
                            k_mean_bmax_end = np.mean(paras["half_saturation"][:,B_max_end_which])
                            vmax_mean_bmax_end = np.mean(paras["max_rate"][:,B_max_end_which])
                            row.append([seed_sim,sim, dom_n, bio_n, t_dom_initial, B_max_initial, B_max_initial_which, k_max_bmax_initial, k_mean_bmax_initial, vmax_max_bmax_initial, vmax_mean_bmax_initial, B_max_all_time, B_max_all_time_which, B_max_all_time_when, k_max_bmax_all_time, k_mean_bmax_all_time, vmax_max_bmax_all_time, vmax_mean_bmax_all_time, B_max_end, B_max_end_which, k_max_bmax_end, k_mean_bmax_end, vmax_max_bmax_end, vmax_mean_bmax_end])
        hr.close()

    micro_char = pd.DataFrame.from_records(row, columns = ["Seed", "Sim_series", "carbon_species", "biomass_species", "DOC_initial", "B_initial_max_conc", "B_initial_max_idx", "B_initial_max_k_max", "B_initial_max_k_mean", "B_initial_max_vmax_max", "B_initial_max_vmax_mean", "B_max", "B_max_idx", "B_max_tim", "B_max_k_max", "B_max_k_mean", "B_max_vmax_max", "B_max_vmax_mean", "B_end_max", "B_end_max_idx", "B_end_max_k_max", "B_end_max_k_mean", "B_end_max_vmax_max", "B_end_max_vmax_mean"])

    print("The shape of the dataframe is ", micro_char.shape)
    print("The dataframe contains the following data types ", micro_char.dtypes)

    filename = os.path.join(results_dir, results_filename+"_micro_char.pkl")
    micro_char.to_pickle(filename)
    print ("Diversity with carbon data is saved here ", filename)