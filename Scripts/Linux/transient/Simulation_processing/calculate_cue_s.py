## Import libraries
import os
import numpy as np
import sys
import h5py
import argparse

CLI=argparse.ArgumentParser()
CLI.add_argument(
  "--sim_label",  # name on the CLI - drop the `--` for positional/required parameters
  type=str,
  default="competition_adaptation",  # default if nothing is provided
)
CLI.add_argument(
  "--scenario",  # name on the CLI - drop the `--` for positional/required parameters
  type=str
)

user_args=CLI.parse_args()
scenario_folder = user_args.scenario
filestring = user_args.sim_label + '_carbon_'

def calculate_cue(sim_data):
    x = sim_data['solution']
    C = np.asarray(x['dom'])
    B = np.asarray(x['biomass'])
    cue = np.zeros((np.shape(B)[0],))
    for t in list(range(np.shape(B)[0])):
        delb = np.sum(B[t,:]) - np.sum(B[t-1,:])
        delc = np.sum(C[t,:]) - np.sum(C[t-1,:])
        cue = delb/delc
    return cue

def calculate_s_os_func_div (sim_data):
    #extract solution and parameters
    dom_n, bio_n = sim_data['species_number']
    x = sim_data['solution']
    C = np.asarray(x['dom'])
    B = np.asarray(x['biomass'])
    paras = sim_data['parameters']
    uptake_coeff = np.asarray(paras['carbon_uptake'])
    v_enz = np.tile(np.asarray(paras['exo_enzyme_rate']), (dom_n,1))
    v_params = np.asarray(paras['max_rate'])
    #Calculate Shannon diversity in time    
    proportion = B/B.sum(axis=1, keepdims = True)
    s = -np.sum(proportion*np.log(proportion), axis = 1)
    #Calculate mean oxidation state of carbon pool weighted by carbon concentration of each pool
    C_sum = np.sum(C, axis = 1)
    mean_wt_os = np.sum(np.array(C)*np.array(paras["oxidation_state"]).T, axis = 1)/C_sum
    #calculate parameters of carbon uptake and consumption as a bulk metric
    para_prod_sum = np.sum(v_enz*uptake_coeff*v_params,axis=0)
    # calculate variance of parameters and weight it by biomass in time
    mean_para_prod = np.mean(para_prod_sum)
    diff_para = para_prod_sum - mean_para_prod
    func_div = np.sum(B*(diff_para**2),axis=1)/np.sum(B, axis = 1)
    return s, mean_wt_os, func_div

## LOAD RESULTS

project_dir = os.path.join('/proj', 'hs_micro_div_072022', 'Project_data', 'transient', scenario_folder)
results_dir = os.path.join(project_dir, "results")
seed_sim_list = [610229235, 983307757, 643338060, 714504443, 277077803, 898393994, 420,13012022,13061989]

cn_list = [3,6,12,18]
bio_n_series = [4,8,12,16,20,24,28,32]
ip = 0
init_dom_list = [1000,2000,5000,10000,15000]
filestring = sys.argv[1] + '_carbon_' #null

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
                for base in ["b_1","b_2", "b_3", "b_4", "b_5"]:
                    if base == "b_1":
                        base_case = "b_1_all_"
                        sim_data = hrw[base_case][c_b][dom_init][seed_all]
                        cue_calc = calculate_cue(sim_data)
                        s_calc, os_calc, fdiv_calc = calculate_s_os_func_div(sim_data)
                        for data_label, new_data in zip(["cue", "shannon", "carbon_os", "func_div"], [cue_calc, s_calc, os_calc, fdiv_calc]):
                            dataset_name = base_case+"/"+c_b+"/"+dom_init+"/"+seed_all+"/solution"+"/"+data_label
                            if dataset_name in hrw.keys():
                                del hrw[base_case][c_b][dom_init][seed_all]['solution'][data_label]
                            hrw.create_dataset(dataset_name, data=new_data)
                    else:
                        for label in ["a", "b", "c","d","e"]:
                            sim = base + "_" + label + "_"
                            sim_data = hrw[sim][c_b][dom_init][seed_all]
                            cue_calc = calculate_cue(sim_data)
                            s_calc, os_calc, fdiv_calc = calculate_s_os_func_div(sim_data)
                            for data_label, new_data in zip(["cue", "shannon", "carbon_os", "func_div"], [cue_calc, s_calc, os_calc, fdiv_calc]):
                                dataset_name = sim+"/"+c_b+"/"+dom_init+"/"+seed_all+"/solution"+"/"+data_label
                                if dataset_name in hrw.keys():
                                    del hrw[sim][c_b][dom_init][seed_all]['solution'][data_label]
                                hrw.create_dataset(dataset_name, data=new_data)
        hrw.close()