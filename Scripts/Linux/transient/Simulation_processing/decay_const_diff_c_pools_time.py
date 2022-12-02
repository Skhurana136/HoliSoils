## Import libraries
import os
import pandas as pd
import numpy as np
import sys
import h5py 

def derive_t_loss(sim_data, loss_criteria):
    paras = sim_data['parameters']
    os_i = np.asarray(paras["oxidation_state"])
    c_os_less_0 = np.where(os_i<0.)
    c_os_eq_0 = np.where(os_i==0.)
    c_os_gr_0 = np.where(os_i>0.)
    x = sim_data['solution']
    C = np.asarray(x['dom'])
    DOC = np.sum(C, axis = 1)
    C_less_0 = np.sum(np.sum(C[:,c_os_less_0], axis = 1),axis=1)
    C_eq_0 = np.sum(np.sum(C[:,c_os_eq_0], axis = 1),axis=1)
    C_gr_0 = np.sum(np.sum(C[:,c_os_gr_0], axis = 1),axis=1)
    carbon_initial = np.asarray(init['dom'])
    DOC_i = np.sum(carbon_initial)
    C_less_0_i = np.sum(carbon_initial[c_os_less_0])
    C_eq_0_i = np.sum(carbon_initial[c_os_eq_0])
    C_gr_0_i = np.sum(carbon_initial[c_os_gr_0])
    
    #DOC:
    t10_doc = np.argwhere(np.round_(DOC/DOC_i, decimals = 2)==loss_criteria[0])
    if t10_doc.size > 0:
        t10_doc_val = t10_doc[0][0]
        t20_doc_i = DOC[t10_doc_val]
        t20_doc = np.argwhere(np.round_(DOC[t10_doc_val:]/t20_doc_i, decimals = 2)==loss_criteria[0])
        if t20_doc.size>0:
            t20_doc_val = t20_doc[0][0]
            t30_doc_i = DOC[(t10_doc_val+t20_doc_val)]
            t30_doc = np.argwhere(np.round_(DOC[(t10_doc_val+t20_doc_val):]/t30_doc_i, decimals = 2)==loss_criteria[0])
            if t30_doc.size>0:
                t30_doc_val = t30_doc[0][0]
                t40_doc_i = DOC[(t10_doc_val+t20_doc_val+t30_doc_val)]
                t40_doc = np.argwhere(np.round_(DOC[(t10_doc_val+t20_doc_val+t30_doc_val):]/t40_doc_i, decimals = 2)==loss_criteria[0])
                if t40_doc.size>0:
                    t40_doc_val = t40_doc[0][0]
                    t50_doc_i = DOC[(t10_doc_val+t20_doc_val+t30_doc_val+t40_doc_val)]
                    t50_doc = np.argwhere(np.round_(DOC[(t10_doc_val+t20_doc_val+t30_doc_val+t40_doc_val):]/t50_doc_i, decimals = 2)==loss_criteria[0])
                    if t50_doc.size>0:
                        t50_doc_val = t50_doc[0][0]
                        t60_doc_i = DOC[(t10_doc_val+t20_doc_val+t30_doc_val+t40_doc_val+t50_doc_val)]
                        t60_doc = np.argwhere(np.round_(DOC[(t10_doc_val+t20_doc_val+t30_doc_val+t40_doc_val+t50_doc_val):]/t60_doc_i, decimals = 2)==loss_criteria[0])
                        if t60_doc.size>0:
                            t60_doc_val = t60_doc[0][0]
                        else:
                            t60_doc_val = "NA"
                    else:
                        t50_doc_val,t60_doc_val="NA","NA"
                else:
                    t40_doc_val,t50_doc_val,t60_doc_val = "NA","NA","NA"
            else:
                t30_doc_val,t40_doc_val,t50_doc_val,t60_doc_val = "NA","NA","NA","NA"
        else:
            t20_doc_val,t30_doc_val,t40_doc_val,t50_doc_val,t60_doc_val = "NA","NA","NA","NA","NA"
    else:
        t10_doc_val,t20_doc_val,t30_doc_val,t40_doc_val,t50_doc_val,t60_doc_val = "NA","NA","NA","NA","NA","NA"
    #reduced C:
    if C_less_0_i>0:
        t10_c_less_0 = np.argwhere(np.round_(C_less_0/C_less_0_i, decimals = 2)==loss_criteria[0])
        if t10_c_less_0.size > 0:
            t10_c_less_0_val = t10_c_less_0[0][0]
            t20_c_less_0_i = DOC[t10_c_less_0_val]
            t20_c_less_0 = np.argwhere(np.round_(DOC[t10_c_less_0_val:]/t20_c_less_0_i, decimals = 2)==loss_criteria[0])
            if t20_c_less_0.size>0:
                t20_c_less_0_val = t20_c_less_0[0][0]
                t30_c_less_0_i = DOC[(t10_c_less_0_val+t20_c_less_0_val)]
                t30_c_less_0 = np.argwhere(np.round_(DOC[(t10_c_less_0_val+t20_c_less_0_val):]/t30_c_less_0_i, decimals = 2)==loss_criteria[0])
                if t30_c_less_0.size>0:
                    t30_c_less_0_val = t30_c_less_0[0][0]
                    t40_c_less_0_i = DOC[(t10_c_less_0_val+t20_c_less_0_val+t30_c_less_0_val)]
                    t40_c_less_0 = np.argwhere(np.round_(DOC[(t10_c_less_0_val+t20_c_less_0_val+t30_c_less_0_val):]/t40_c_less_0_i, decimals = 2)==loss_criteria[0])
                    if t40_c_less_0.size>0:
                        t40_c_less_0_val = t40_c_less_0[0][0]
                        t50_c_less_0_i = DOC[(t10_c_less_0_val+t20_c_less_0_val+t30_c_less_0_val+t40_c_less_0_val)]
                        t50_c_less_0 = np.argwhere(np.round_(DOC[(t10_c_less_0_val+t20_c_less_0_val+t30_c_less_0_val+t40_c_less_0_val):]/t50_c_less_0_i, decimals = 2)==loss_criteria[0])
                        if t50_c_less_0.size>0:
                            t50_c_less_0_val = t50_c_less_0[0][0]
                            t60_c_less_0_i = DOC[(t10_c_less_0_val+t20_c_less_0_val+t30_c_less_0_val+t40_c_less_0_val+t50_c_less_0_val)]
                            t60_c_less_0 = np.argwhere(np.round_(DOC[(t10_c_less_0_val+t20_c_less_0_val+t30_c_less_0_val+t40_c_less_0_val+t50_c_less_0_val):]/t60_c_less_0_i, decimals = 2)==loss_criteria[0])
                            if t60_c_less_0.size>0:
                                t60_c_less_0_val = t60_c_less_0[0][0]
                            else:
                                t60_c_less_0_val = "NA"
                        else:
                            t50_c_less_0_val,t60_c_less_0_val="NA","NA"
                    else:
                        t40_c_less_0_val,t50_c_less_0_val,t60_c_less_0_val = "NA","NA","NA"
                else:
                    t30_c_less_0_val,t40_c_less_0_val,t50_c_less_0_val,t60_c_less_0_val = "NA","NA","NA","NA"
            else:
                t20_c_less_0_val,t30_c_less_0_val,t40_c_less_0_val,t50_c_less_0_val,t60_c_less_0_val = "NA","NA","NA","NA","NA"
        else:
            t10_c_less_0_val,t20_c_less_0_val,t30_c_less_0_val,t40_c_less_0_val,t50_c_less_0_val,t60_c_less_0_val = "NA","NA","NA","NA","NA","NA"
    else:
        t10_c_less_0_val,t20_c_less_0_val,t30_c_less_0_val,t40_c_less_0_val,t50_c_less_0_val,t60_c_less_0_val = "NA","NA","NA","NA","NA","NA"
    
    #necromass:
    #t10_c_eq_0 = np.argwhere(np.round_(C_eq_0/C_eq_0_i, decimals = 2)==loss_criteria[0])
    #if t10_c_eq_0.size > 0:
    #    t10_c_eq_0_val = t10_c_eq_0[0][0]
    #    t20_c_eq_0_i = DOC[t10_c_eq_0]
    #    t20_c_eq_0 = np.argwhere(np.round_(DOC[t10_c_eq_0_val:]/t20_c_eq_0_i, decimals = 2)==loss_criteria[0])
    #    if t20_c_eq_0.size>0:
    #        t20_c_eq_0_val = t20_c_eq_0[0][0]
    #        t30_c_eq_0_i = DOC[(t10_c_eq_0_val+t20_c_eq_0_val)]
    #        t30_c_eq_0 = np.argwhere(np.round_(DOC[(t10_c_eq_0_val+t20_c_eq_0_val):]/t30_c_eq_0_i, decimals = 2)==loss_criteria[0])
    #        if t30_c_eq_0.size>0:
    #            t30_c_eq_0_val = t30_c_eq_0[0][0]
    #        else:
    #            t30_c_eq_0_val="NA"
    #   else:
    #       t20_c_eq_0_val="NA"
    #        t30_c_eq_0_val="NA"            
    #else:
    t10_c_eq_0_val,t20_c_eq_0_val,t30_c_eq_0_val,t40_c_eq_0_val,t50_c_eq_0_val,t60_c_eq_0_val = "NA","NA","NA","NA","NA","NA"

    #oxidized C:
    if C_gr_0_i>0:
        t10_c_gr_0 = np.argwhere(np.round_(C_gr_0/C_gr_0_i, decimals = 2)==loss_criteria[0])
        if t10_c_gr_0.size > 0:
            t10_c_gr_0_val = t10_c_gr_0[0][0]
            t20_c_gr_0_i = DOC[t10_c_gr_0_val]
            t20_c_gr_0 = np.argwhere(np.round_(DOC[t10_c_gr_0_val:]/t20_c_gr_0_i, decimals = 2)==loss_criteria[0])
            if t20_c_gr_0.size>0:
                t20_c_gr_0_val = t20_c_gr_0[0][0]
                t30_c_gr_0_i = DOC[(t10_c_gr_0_val+t20_c_gr_0_val)]
                t30_c_gr_0 = np.argwhere(np.round_(DOC[(t10_c_gr_0_val+t20_c_gr_0_val):]/t30_c_gr_0_i, decimals = 2)==loss_criteria[0])
                if t30_c_gr_0.size>0:
                    t30_c_gr_0_val = t30_c_gr_0[0][0]
                    t40_c_gr_0_i = DOC[(t10_c_gr_0_val+t20_c_gr_0_val+t30_c_gr_0_val)]
                    t40_c_gr_0 = np.argwhere(np.round_(DOC[(t10_c_gr_0_val+t20_c_gr_0_val+t30_c_gr_0_val):]/t40_c_gr_0_i, decimals = 2)==loss_criteria[0])
                    if t40_c_gr_0.size>0:
                        t40_c_gr_0_val = t40_c_gr_0[0][0]
                        t50_c_gr_0_i = DOC[(t10_c_gr_0_val+t20_c_gr_0_val+t30_c_gr_0_val+t40_c_gr_0_val)]
                        t50_c_gr_0 = np.argwhere(np.round_(DOC[(t10_c_gr_0_val+t20_c_gr_0_val+t30_c_gr_0_val+t40_c_gr_0_val):]/t50_c_gr_0_i, decimals = 2)==loss_criteria[0])
                        if t50_c_gr_0.size>0:
                            t50_c_gr_0_val = t50_c_gr_0[0][0]
                            t60_c_gr_0_i = DOC[(t10_c_gr_0_val+t20_c_gr_0_val+t30_c_gr_0_val+t40_c_gr_0_val+t50_c_gr_0_val)]
                            t60_c_gr_0 = np.argwhere(np.round_(DOC[(t10_c_gr_0_val+t20_c_gr_0_val+t30_c_gr_0_val+t40_c_gr_0_val+t50_c_gr_0_val):]/t60_c_gr_0_i, decimals = 2)==loss_criteria[0])
                            if t60_c_gr_0.size>0:
                                t60_c_gr_0_val = t60_c_gr_0[0][0]
                            else:
                                t60_c_gr_0_val = "NA"
                        else:
                            t50_c_gr_0_val,t60_c_gr_0_val="NA","NA"
                    else:
                        t40_c_gr_0_val,t50_c_gr_0_val,t60_c_gr_0_val = "NA","NA","NA"
                else:
                    t30_c_gr_0_val,t40_c_gr_0_val,t50_c_gr_0_val,t60_c_gr_0_val = "NA","NA","NA","NA"
            else:
                t20_c_gr_0_val,t30_c_gr_0_val,t40_c_gr_0_val,t50_c_gr_0_val,t60_c_gr_0_val = "NA","NA","NA","NA","NA"
        else:
            t10_c_gr_0_val,t20_c_gr_0_val,t30_c_gr_0_val,t40_c_gr_0_val,t50_c_gr_0_val,t60_c_gr_0_val = "NA","NA","NA","NA","NA","NA"
    else:
        t10_c_gr_0_val,t20_c_gr_0_val,t30_c_gr_0_val,t40_c_gr_0_val,t50_c_gr_0_val,t60_c_gr_0_val = "NA","NA","NA","NA","NA","NA"
    
    results_list = [t10_doc_val,t10_c_less_0_val,t10_c_eq_0_val,t10_c_gr_0_val,t20_doc_val,t20_c_less_0_val,t20_c_eq_0_val,t20_c_gr_0_val,t30_doc_val,t30_c_less_0_val,t30_c_eq_0_val,t30_c_gr_0_val,t40_doc_val,t40_c_less_0_val,t40_c_eq_0_val,t40_c_gr_0_val,t50_doc_val,t50_c_less_0_val,t50_c_eq_0_val,t50_c_gr_0_val,t60_doc_val,t60_c_less_0_val,t60_c_eq_0_val,t60_c_gr_0_val]
    #print(results_list)
    return (results_list)

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
loss_crit_1 = 0.9
loss_crit_2 = 0.8
loss_crit_3 = 0.5
result_fstring = "_loss_temporal"

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
                    dom_n, bio_n = sim_data['species_number']
                    carbon_initial = np.asarray(init['dom'])
                    DOC_i = np.sum(carbon_initial)
                    c_loss_tim_points=derive_t_loss(sim_data, [loss_crit_1,loss_crit_2,loss_crit_3])
                    for c_pool, t10, t20, t30, t40, t50, t60 in zip(["DOC","reduced_C", "necromass", "oxidized_C"],c_loss_tim_points[:4],c_loss_tim_points[4:8],c_loss_tim_points[8:12],c_loss_tim_points[12:16],c_loss_tim_points[16:20],c_loss_tim_points[20:24]):
                        t10_base,t20_base, t30_base, t40_base,t50_base,t60_base = t10,t20,t30,t40,t50,t60
                        all_results_dictionary[dictionary_iter]={"Seed":seed_sim, "Sim_series":base_case, "carbon_species":dom_n, "biomass_species":bio_n, "DOC_initial":DOC_i,"C_pool": c_pool, "T10": t10, "T10_base":t10_base, "T20": t20, "T2base":t20_base,"T30": t30, "T30_base":t30_base,"T40": t40, "T40_base":t40_base,"T50": t50, "T50_base":t50_base,"T60": t60, "T60_base":t60_base}
                        dictionary_iter+=1
                for baseline in ["b_2", "b_3", "b_4","b_5"]:
                    for label in ["a", "b", "c","d","e"]:
                        sim = baseline + "_" + label + "_"
                        if hr[sim][c_b][dom_init][seed_all]:
                            sim_data = hr[sim][c_b][dom_init][seed_all]
                            sim_seed_details = sim+"/"+c_b+"/"+dom_init+"/"+seed_all
                            init = sim_data['initial_conditions']
                            carbon_initial = np.asarray(init['dom'])
                            DOC_i = np.sum(carbon_initial)
                            c_loss_tim_points=derive_t_loss(sim_data, [loss_crit_1,loss_crit_2,loss_crit_3])
                            for c_pool, t10, t20, t30, t40, t50, t60 in zip(["DOC","reduced_C", "necromass", "oxidized_C"],c_loss_tim_points[:4],c_loss_tim_points[4:8],c_loss_tim_points[8:12],c_loss_tim_points[12:16],c_loss_tim_points[16:20],c_loss_tim_points[20:24]):
                                all_results_dictionary[dictionary_iter]={"Seed":seed_sim, "Sim_series":sim, "carbon_species":dom_n, "biomass_species":bio_n, "DOC_initial":DOC_i,"C_pool": c_pool, "T10": t10, "T10_base":t10_base, "T20": t20, "T20_base":t20_base,"T30": t30, "T30_base":t30_base,"T40": t40, "T40_base":t40_base,"T50": t50, "T50_base":t50_base,"T60": t60, "T60_base":t60_base}
                                dictionary_iter+=1
        hr.close()
    print(dictionary_iter)
    decay_const_c_pools_data = pd.DataFrame.from_dict(all_results_dictionary,orient='index')#, columns = ["Seed", "Sim_series", "carbon_species", "biomass_species", "DOC_initial", "C_pool", "T10", "T10_base","T20", "T20_base","T50", "T50_base"])
    filename = os.path.join(results_dir, results_filename+result_fstring+"_temporal_decay_const_c_pools_data.pkl")
    decay_const_c_pools_data.to_pickle(filename)
    print ("Temporal decay constant for diff carbon pools are saved here ", filename)