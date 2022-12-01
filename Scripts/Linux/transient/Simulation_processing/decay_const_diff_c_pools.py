## Import libraries
import os
import pandas as pd
import numpy as np
import sys
import h5py 


## LOAD RESULTS
#project_dir = os.path.join("D:\Projects", "HoliSoils","data","transient", sys.argv[1])
project_dir = os.path.join('/proj', 'hs_micro_div_072022', 'Project_data', 'transient', sys.argv[1])
results_dir = os.path.join(project_dir, "results")
filestring =  'competition_adaptation_carbon_' #null
seed_sim_list = [610229235, 983307757, 643338060, 714504443, 277077803, 898393994, 420,13012022,13061989]

cn_list = [3,6,12,18]
bio_n_series = [4,8,12,16,20,24,28,32]
ip = 0
init_dom_list = [1000,2000,5000,10000,15000]
transient_switch = 0
input_factor = transient_switch*5/365
loss_crit = 0.9#sys.argv[2]#0.37#0.63
result_fstring = "_loss_"+str(loss_crit)

all_results_dictionary={}
dictionary_iter = 0
for c_n in cn_list:
    row = []
    c_b_row = []
    results_filename = os.path.join(results_dir, filestring + str(c_n))

    for seed_sim in seed_sim_list:
        # Load all datasets and save their Shannon and diversity indices in a dataframe
        seed_all = 'seed_'+str(seed_sim)
        
        details_subfolder = filestring + str(c_n) + '_'+str(seed_sim) + '_ip_' + str(ip)
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
                #print(c_n, seed_sim, b_n, t_dom_initial)
                if hr[base_case][c_b][dom_init][seed_all]:
                    sim_data = hr[base_case][c_b][dom_init][seed_all]
                    sim_seed_details = base_case+"/"+c_b+"/"+dom_init+"/"+seed_all
                    init = sim_data['initial_conditions']
                    paras = sim_data['parameters']
                    os_i = np.asarray(paras["oxidation_state"])
                    c_os_less_0 = np.where(os_i<0.)
                    c_os_eq_0 = np.where(os_i==0.)
                    c_os_gr_0 = np.where(os_i>0.)
                    dom_n, bio_n = sim_data['species_number']
                    x = sim_data['solution']
                    C = np.asarray(x['dom'])
                    C_less_0 = np.sum(C[:,c_os_less_0], axis = 1)
                    C_eq_0 = np.sum(C[:,c_os_eq_0], axis = 1)
                    C_gr_0 = np.sum(C[:,c_os_gr_0], axis=1)
                    carbon_initial = np.asarray(init['dom'])
                    DOC_i = np.sum(carbon_initial)
                    C_less_0_i = np.sum(carbon_initial[c_os_less_0])
                    C_eq_0_i = np.sum(carbon_initial[c_os_eq_0])
                    C_gr_0_i = np.sum(carbon_initial[c_os_gr_0])
                    #Time taken for each C pool to be 90% of initial value
                    t10_c_less_0 = np.argwhere(np.round_(C_less_0/C_less_0_i, decimals = 2)==loss_crit)
                    t10_c_eq_0 = np.argwhere(np.round_(C_eq_0/C_eq_0_i, decimals = 2)==loss_crit)
                    t10_c_gr_0 = np.argwhere(np.round_(C_gr_0/C_gr_0_i, decimals = 2)==loss_crit)
                    for c_pool, t10_arr in zip(["reduced_C", "necromass", "oxidized_C"],[t10_c_less_0,t10_c_eq_0,t10_c_gr_0]):
                        if t10_arr.size > 0:
                            t10 = t10_arr[0][0]
                        else:
                            t10 = "NA"
                        t10_base = t10
                        all_results_dictionary[dictionary_iter]={"Seed":seed_sim, "Sim_series":base_case, "carbon_species":dom_n, "biomass_species":bio_n, "DOC_initial":DOC_i,"C_pool": c_pool, "T_10": t10, "T10_base":t10_base}
                        dictionary_iter+=1
                for baseline in ["b_2", "b_3", "b_4","b_5"]:
                    for label in ["a", "b", "c","d","e"]:
                        sim = baseline + "_" + label + "_"
                        if hr[sim][c_b][dom_init][seed_all]:
                            sim_data = hr[sim][c_b][dom_init][seed_all]
                            sim_seed_details = sim+"/"+c_b+"/"+dom_init+"/"+seed_all
                            init = sim_data['initial_conditions']
                            paras = sim_data['parameters']
                            os_i = np.asarray(paras["oxidation_state"])
                            c_os_less_0 = np.where(os_i<0.)
                            c_os_eq_0 = np.where(os_i==0.)
                            c_os_gr_0 = np.where(os_i>0.)
                            dom_n, bio_n = sim_data['species_number']
                            x = sim_data['solution']
                            C = np.asarray(x['dom'])
                            C_less_0 = np.sum(C[:,c_os_less_0], axis = 1)
                            C_eq_0 = np.sum(C[:,c_os_eq_0], axis = 1)
                            C_gr_0 = np.sum(C[:,c_os_gr_0], axis=1)
                            carbon_initial = np.asarray(init['dom'])
                            DOC_i = np.sum(carbon_initial)
                            C_less_0_i = np.sum(carbon_initial[c_os_less_0])
                            C_eq_0_i = np.sum(carbon_initial[c_os_eq_0])
                            C_gr_0_i = np.sum(carbon_initial[c_os_gr_0])
                            #Time taken for each C pool to be 90% of initial value
                            t10_c_less_0 = np.argwhere(np.round_(C_less_0/C_less_0_i, decimals = 2)==loss_crit)
                            t10_c_eq_0 = np.argwhere(np.round_(C_eq_0/C_eq_0_i, decimals = 2)==loss_crit)
                            t10_c_gr_0 = np.argwhere(np.round_(C_gr_0/C_gr_0_i, decimals = 2)==loss_crit)
                            for c_pool, t10_arr in zip(["reduced_C", "necromass", "oxidized_C"],[t10_c_less_0,t10_c_eq_0,t10_c_gr_0]):
                                if t10_arr.size > 0:
                                    t10 = t10_arr[0][0]
                                else:
                                    t10 = "NA"
                                all_results_dictionary[dictionary_iter]={"Seed":seed_sim, "Sim_series":sim, "carbon_species":dom_n, "biomass_species":bio_n, "DOC_initial":DOC_i,"C_pool": c_pool, "T_10": t10, "T10_base":t10_base}
                                dictionary_iter+=1
        hr.close()
    print(dictionary_iter)
    decay_const_c_pools_data = pd.DataFrame.from_dict(all_results_dictionary)#, columns = ["Seed", "Sim_series", "carbon_species", "biomass_species", "DOC_initial", "C_pool", "T_10", "T_10_B1"])
    filename = os.path.join(results_dir, results_filename+result_fstring+"_decay_const_c_pools_data.pkl")
    decay_const_c_pools_data.to_pickle(filename)
    print ("Decay constant for diff carbon pools are saved here ", filename)