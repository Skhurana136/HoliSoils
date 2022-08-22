## Import libraries
import os
import pandas as pd
import numpy as np
import sys
import h5py 


## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient", "seq_loss")
results_dir = os.path.join(project_dir, "results")
filestring = sys.argv[1] + '_carbon_' #null
seed_sim_list = [610229235, 983307757, 643338060, 714504443, 277077803]#, 898393994, 420,13012022,13061989]

cn_list = [3,6,12,18]
bio_n_series = [4]#[4,8,12,16,20,24,28,32]
ip = 0
init_dom_list = [1000,2000,5000,10000,15000]
transient_switch = 0
input_factor = transient_switch*5/365
loss_crit = 0.9#sys.argv[2]#0.37#0.63
result_fstring = "_loss_"+str(loss_crit)

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
                    dom_n, bio_n = sim_data['species_number']
                    x = sim_data['solution']
                    sim_status = seed_details[sim_seed_details]['sim_status']
                    C = np.asarray(x['dom'])
                    B = np.asarray(x['biomass'])
                    data = np.append(C,B,axis=1)
                    DOC = np.sum(C,axis=1)
                    proportion = B/B.sum(axis=1, keepdims = True)
                    S = -np.sum(proportion*np.log(proportion), axis = 1)
                    TOC = np.sum(C,axis=1) + np.sum(B, axis=1)
                    initial_data = sim_data['initial_conditions']
                    carbon_initial = np.asarray(initial_data['dom'])
                    #Carbon removal
                    #Carbon at end
                    DOC_end = DOC[-1]
                    DOC_mid = np.take(DOC, DOC.size//2)
                    #Initial carbon
                    DOC_i = np.sum(carbon_initial)
                    #Total input of DOC
                    doc_input_i = DOC_i + doc_input
                    #Time taken for DOC to be 50% of initial DOC
                    t50_arr = np.argwhere(np.round_(DOC/DOC_i, decimals = 2)==loss_crit)
                    if t50_arr.size > 0:
                        t50 = t50_arr[0][0]
                    else:
                        t50 = "NA"
                    t50_b1 = t50
                    #At baseline T50
                    #Biomass and Shannon
                    Biomass = np.sum(B, axis = 1)
                    biomass_initial = np.asarray(initial_data['biomass'])
                    bio_i = np.sum(biomass_initial)
                    bio_end = Biomass[-1]
                    bio_mid = np.take(Biomass, DOC.size//2)
                    proportion_i = biomass_initial/np.sum(biomass_initial)
                    #Maximum Biomass
                    B_max = np.max(B)             
                    #Initial Shannon
                    S_i = -np.sum(proportion_i*np.log(proportion_i))
                    #Maximum Shannon
                    S_max = np.max(S)
                    #Steady state Shannon
                    S_end = S[-1]
                    S_mid = np.take(S, DOC.size//2)
                    if t50_b1 != "NA":
                        #DOC at t50 of B1
                        DOC_t50_b1 = DOC[t50_b1]
                        #Biomass at t50 of B1
                        B_t50_b1 = Biomass[t50_b1]
                        #Biomass at t50
                        B_t50 = B_t50_b1
                        #Shannon at t50
                        S_t50 = S[t50]
                        #Shannon at t50 of B1
                        S_t50_b1 = S_t50
                    else:
                        DOC_t50_b1, B_t50_b1, B_t50, S_t50, S_t50_b1 = "NA", "NA", "NA", "NA", "NA"      
                    row.append([seed_sim,base_case, dom_n, bio_n, DOC_i, doc_input_i, DOC_end, t50, t50_b1, DOC_t50_b1, S_i, S_end, S_max, S_t50, S_t50_b1, bio_i, bio_end, B_max, B_t50, B_t50_b1, DOC_mid, S_mid, bio_mid, sim_status])
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
                    paras = sim_data['parameters']
                    os_i = paras["oxidation_state"]
                    vmax_base = paras['max_rate']
                    vmax_base_mean = np.sum(np.sum(vmax_base, axis = 1)*C[0,:])/DOC_i
                    vmax_base_sum = np.sum(vmax_base)
                    vmax_max_base = np.max(vmax_base)
                    vmax_max_base_os = os_i[np.where(vmax_base==np.max(vmax_base))[0][0]]
                    mean_os_initial_base = np.sum(os_i*C[0,:])/DOC_i
                    mean_os_end_base = np.sum(os_i*C[-1,:])/np.sum(C[-1,:])
                    c_b_row.append([seed_sim,base_case, dom_n, bio_n, DOC_i, DOC_end, vmax_base_mean, vmax_base_mean, vmax_base_mean/vmax_base_mean, vmax_base_sum, vmax_base_sum, vmax_max_base, vmax_max_base, vmax_max_base_os,vmax_max_base_os, mean_os_initial_base, mean_os_end_base, mean_os_end_base, S_i, C_B_r_min, C_B_r_max, C_B_r_median, C_B_r_mean, C_B_r_i, C_B_r_end, C_i_max_o_s, C_end_max_o_s, C_i_max_o_s_pc, C_end_max_o_s_pc])
                for baseline in ["b_2", "b_3", "b_4","b_5"]:
                    for label in ["a", "b", "c","d","e"]:
                        sim = baseline + "_" + label + "_"
                        if hr[sim][c_b][dom_init][seed_all]:
                            sim_data = hr[sim][c_b][dom_init][seed_all]
                            sim_seed_details = sim+"/"+c_b+"/"+dom_init+"/"+seed_all
                            init = sim_data['initial_conditions']
                            paras = sim_data['parameters']
                            dom_n, bio_n = sim_data['species_number']
                            x = sim_data['solution']
                            sim_status = seed_details[sim_seed_details]['sim_status']
                            C = np.asarray(x['dom'])
                            B = np.asarray(x['biomass'])
                            data = np.append(C,B,axis=1)
                            DOC = np.sum(C,axis=1)
                            proportion = B/B.sum(axis=1, keepdims = True)
                            S = -np.sum(proportion*np.log(proportion), axis = 1)
                            TOC = np.sum(C,axis=1) + np.sum(B, axis=1)
                            initial_data = sim_data['initial_conditions']
                            carbon_initial = np.asarray(initial_data['dom'])
                            #Carbon removal
                            #Carbon at end
                            DOC_end = DOC[-1]
                            DOC_mid = np.take(DOC, DOC.size//2)
                            #Initial carbon
                            DOC_i = np.sum(carbon_initial)
                            #Total input of DOC
                            doc_input_i = DOC_i + doc_input
                            #Time taken for DOC to be 50% of initial DOC
                            t50_arr = np.argwhere(np.round_(DOC/DOC_i, decimals = 2)==loss_crit)
                            if t50_arr.size > 0:
                                t50 = t50_arr[0][0]
                            else:
                                t50 = "NA" 
                            #Biomass and Shannon
                            Biomass = np.sum(B, axis = 1)
                            biomass_initial = np.asarray(initial_data['biomass'])
                            bio_i = np.sum(biomass_initial)
                            bio_end = Biomass[-1]
                            bio_mid = np.take(Biomass, DOC.size//2)
                            proportion_i = biomass_initial/np.sum(biomass_initial)
                            #Maximum Biomass
                            B_max = np.max(B)         
                            #Initial Shannon
                            S_i = -np.sum(proportion_i*np.log(proportion_i))
                            #Maximum Shannon
                            S_max = np.max(S)
                            #Steady state Shannon
                            S_end = S[-1]  
                            S_mid = np.take(S, DOC.size//2)
                            if (t50_b1 != "NA") and (DOC.size>=t50_b1):
                                #DOC at t50 of B1
                                DOC_t50_b1 = DOC[t50_b1]
                                #Biomass at t50 of B1
                                B_t50_b1 = Biomass[t50_b1]
                                #Shannon at t50 of B1
                                S_t50_b1 = S_t50
                            else:
                                DOC_t50_b1, B_t50_b1, S_t50_b1 = "NA", "NA", "NA"
                            if t50 != "NA":
                                #Biomass at t50
                                B_t50 = Biomass[t50]   
                                #Shannon at t50
                                S_t50 = S[t50]
                            else:
                                B_t50, S_t50 = "NA", "NA"
                            row.append([seed_sim,sim, dom_n, bio_n, DOC_i, doc_input_i, DOC_end, t50, t50_b1, DOC_t50_b1, S_i, S_end, S_max, S_t50, S_t50_b1, bio_i, bio_end, B_max, B_t50, B_t50_b1, DOC_mid, S_mid, bio_mid, sim_status])
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
                            paras = sim_data['parameters']
                            os_i = paras["oxidation_state"]
                            vmax = paras['max_rate']
                            vmax_mean = np.sum(np.sum(vmax, axis = 1)*C[0,:])/DOC_i
                            vmax_sum = np.sum(vmax)
                            vmax_max = np.max(vmax)
                            vmax_max_os = os_i[np.where(vmax==np.max(vmax))[0][0]]
                            mean_os_initial = np.sum(os_i*C[0,:])/DOC_i
                            mean_os_end = np.sum(os_i*C[-1,:])/np.sum(C[-1,:])
                            c_b_row.append([seed_sim,sim, dom_n, bio_n, DOC_i, DOC_end, vmax_base_mean, vmax_mean, vmax_mean/vmax_base_mean, vmax_base_sum, vmax_sum, vmax_max_base, vmax_max, vmax_max_base_os,vmax_max_os, mean_os_initial_base, mean_os_end_base, mean_os_end, S_i, C_B_r_min, C_B_r_max, C_B_r_median, C_B_r_mean, C_B_r_i, C_B_r_end, C_i_max_o_s, C_end_max_o_s, C_i_max_o_s_pc, C_end_max_o_s_pc])
        hr.close()

    diversity_data = pd.DataFrame.from_records(row, columns = ["Seed", "Sim_series", "carbon_species", "biomass_species", "DOC_initial", "DOC_input", "DOC_end", "T_50", "T_50_B1", "DOC_T50_B1","S_initial", "S_end", "S_max", "S_t_50", "S_t_50_b1", "Biomass_initial", "Biomass_end", "Biomass_max","Biomass_t_50", "Biomass_t_50_b1", "DOC_mid", "Biomass_mid", "S_mid", "Status"])
    c_b_dynamics_data = pd.DataFrame.from_records(c_b_row, columns = ["Seed", "Sim_series", "carbon_species", "biomass_species", "DOC_initial", "DOC_end", "vmax_base", "vmax_mean", "vmax_ratio", "vmax_sum_base", "vmax_sum", "vmax_max_base", "vmax_max", "vmax_max_os_base", "vmax_max_os", "mean_os_initial", "mean_os_end_base", "mean_os_end" , "S_initial", "C_B_min", "C_B_max", "C_B_median", "C_B_mean", "C_B_initial", "C_B_end", "DOC_initial_max_os", "DOC_end_max_os", "DOC_initial_max_os_pc", "DOC_end_max_os_pc"])

    print("The shape of the dataframe is ", diversity_data.shape)
    print("The dataframe contains the following data types ", diversity_data.dtypes)

    filename = os.path.join(results_dir, results_filename+result_fstring+"_diversity_data.pkl")
    diversity_data.to_pickle(filename)
    print ("Diversity with carbon data is saved here ", filename)

    print("The shape of the dataframe is ", c_b_dynamics_data.shape)
    print("The dataframe contains the following data types ", c_b_dynamics_data.dtypes)

    filename = os.path.join(results_dir, results_filename+result_fstring+"_c_b_dynamics_data.pkl")
    c_b_dynamics_data.to_pickle(filename)
    print ("Carbon dynamics are saved here ", filename)