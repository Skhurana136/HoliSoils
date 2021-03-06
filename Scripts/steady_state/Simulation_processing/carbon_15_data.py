## Import libraries
import os
import pandas as pd
import numpy as np

import h5py 

## LOAD RESULTS
#project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
project_dir = "D:\Projects\HoliSoils\data"
results_dir = os.path.join(project_dir, "results")

seed_sim_list = [420,13012022,13061989]

c_n = 15
loss_crit = 0.63
bio_n_series = [3,5,9,12,15,20,30]#,25,30]#[3,5,9,12,15,20,30]#originally: [3,5,9,12,15,20,25,30] bu 24 didn't work
ip = 1#0
input_factor = (1/ip)*5/365
init_dom_list = [1000,2000,5000,10000,15000]
row = []
failed_sim = []
complete_sim = []
results_filename = os.path.join(results_dir, 'carbon_' + str(c_n) + "_ip_"+str(ip)+"_diversity_data.pkl")

for seed_sim in seed_sim_list:
    details_subfolder = 'carbon_' + str(c_n) + '_'+str(seed_sim) + '_ip_' + str(ip)
    simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)

    hr = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'r')

    # Load all datasets and save their Shannon and diversity indices in a dataframe:
    sim_list = list(hr.keys())
    c_b_sim_list = list(hr[sim_list[0]].keys())
    dom_sim = list(hr[sim_list[0]][c_b_sim_list[0]].keys())
    seed_all = 'seed_'+str(seed_sim)

    filename = os.path.join(simulations_dir, "seeds_randoms.pkl")
    seed_details = pd.read_pickle(filename)

    hr = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'r')

    # Load all datasets and save their Shannon and diversity indices in a dataframe:
    sim_list = list(hr.keys())
    c_b_sim_list = list(hr[sim_list[0]].keys())
    dom_sim = list(hr[sim_list[0]][c_b_sim_list[0]].keys())
    #print("DOM initial removed")
    #seed_all = list(hr[sim_list[0]][c_b_sim_list[0]][dom_sim[0]].keys())[0]
    all_data_stored = list(hr[sim_list[0]][c_b_sim_list[0]][dom_sim[0]][seed_all].keys())

    simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
    filename = os.path.join(simulations_dir, "seeds_randoms.pkl")
    seed_details = pd.read_pickle(filename)
    failed_sim = []
    complete_sim = []
    for base in ["b_1"]:
        #print(base)
        for case in ["all"]:
            #print(case)
            for c_b_r in bio_n_series:
                #print(c_b_r)
                for dom_init in init_dom_list:
                    #print(dom_init)
                    sim_id = base + "_" + case + "_/bio_n_" + str(c_b_r) + "/dom_initial_" + str(dom_init) + "/"+seed_all
                    status = seed_details[sim_id]["sim_status"]
                    if status == "failed":
                        print(base, case, c_b_r, dom_init)
                        failed_sim.append([base, case, c_b_r, dom_init])
                    else:
                        complete_sim.append([base, case+"_", c_b_r, dom_init])
    for base in ["b_2", "b_3", "b_4", "b_5"]:
        #print(base)
        for case in ["a", "b", "c", "d", "e"]:
            #print(case)
            for c_b_r in bio_n_series:
                #print(c_b_r)
                for dom_init in init_dom_list:
                    #print(dom_init)
                    sim_id = base + "_" + case + "_/bio_n_" + str(c_b_r) + "/dom_initial_" + str(dom_init) + "/"+seed_all
                    status = seed_details[sim_id]["sim_status"]
                    if status == "failed":
                        print(base, case, c_b_r, dom_init)
                        failed_sim.append([base, case, c_b_r, dom_init])
                    else:
                        complete_sim.append([base, case+"_", c_b_r, dom_init])

    for b_n in bio_n_series:
        for t_dom_initial in init_dom_list:
            base_case = "b_1_all_"
            c_b = "bio_n_"+ str(b_n)
            dom_init = "dom_initial_" + str(t_dom_initial)
            seed_all = 'seed_'+str(seed_sim)
            doc_input = (t_dom_initial)*input_factor
            if hr[base_case][c_b][dom_init][seed_all]:
                sim_data = hr[base_case][c_b][dom_init][seed_all]
                init = sim_data['initial_conditions']
                paras = sim_data['parameters']
                dom_n, bio_n = sim_data['species_number']
                x = sim_data['solution']
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
                proportion_i = biomass_initial/np.sum(biomass_initial)
                #Maximum Biomass
                B_max = np.max(B)             
                #Initial Shannon
                S_i = -np.sum(proportion_i*np.log(proportion_i))
                #Maximum Shannon
                S_max = np.max(S)
                #Steady state Shannon
                S_end = S[-1]
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
                row.append([seed_sim, base_case, dom_n, bio_n, DOC_i, doc_input_i, DOC_end, t50, t50_b1, DOC_t50_b1, S_i, S_end, S_max, S_t50, S_t50_b1, bio_i, bio_end, B_max, B_t50, B_t50_b1])
            for baseline in ["b_2", "b_3", "b_4", "b_5"]:
                for label in ["a", "b", "c", "d", "e"]:
                    sim = baseline + "_" + label + "_"
                    if hr[sim][c_b][dom_init][seed_all]:
                        sim_data = hr[sim][c_b][dom_init][seed_all]
                        init = sim_data['initial_conditions']
                        paras = sim_data['parameters']
                        dom_n, bio_n = sim_data['species_number']
                        x = sim_data['solution']
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
                        proportion_i = biomass_initial/np.sum(biomass_initial)
                        #Maximum Biomass
                        B_max = np.max(B)         
                        #Initial Shannon
                        S_i = -np.sum(proportion_i*np.log(proportion_i))
                        #Maximum Shannon
                        S_max = np.max(S)
                        #Steady state Shannon
                        S_end = S[-1]  
                        if t50_b1 != "NA":
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
                        
                        row.append([seed_sim,sim, dom_n, bio_n, DOC_i, doc_input_i, DOC_end, t50, t50_b1, DOC_t50_b1, S_i, S_end, S_max, S_t50, S_t50_b1, bio_i, bio_end, B_max, B_t50, B_t50_b1])
    hr.close()

diversity_data = pd.DataFrame.from_records(row, columns = ["Seed","Sim_series", "carbon_species", "biomass_species", "DOC_initial", "DOC_input", "DOC_end", "T_50", "T_50_B1", "DOC_T50_B1","S_initial", "S_end", "S_max", "S_t_50", "S_t_50_b1", "Biomass_initial", "Biomass_end", "Biomass_max","Biomass_t_50", "Biomass_t_50_b1"])

print("The shape of the dataframe is ", diversity_data.shape)
print("The dataframe contains the following data types ", diversity_data.dtypes)

filename = os.path.join(results_dir, "diversity_data.pkl")
diversity_data.to_pickle(results_filename)
print ("Diversity with carbon data is saved here ", results_filename)