# ## Import libraries
import numpy as np
import h5py
import os
import pickle
import csv
import math
import random
import sys
import argparse
from scipy.stats import skewnorm

from DS.solvers.diff_eqn_system import ReactionNetwork as rn
from DS.solvers.diff_eqn_system import generate_random_initial_conditions
from DS.solvers.diff_eqn_system import generate_random_boundary_conditions

project_dir = os.path.join('/proj', 'hs_micro_div_072022', 'Project_data', 'transient', 'gen_spec_lognorm_1_5x')

CLI=argparse.ArgumentParser()
CLI.add_argument(
  "--sim_label",  # name on the CLI - drop the `--` for positional/required parameters
  type=str,
  default="null",  # default if nothing is provided
)
CLI.add_argument(
  "--seeds_num",  # name on the CLI - drop the `--` for positional/required parameters
  nargs="*",  # 0 or more values expected => creates a list
  type=int,
  default=[610229235, 983307757, 643338060, 714504443, 277077803, 898393994, 420,13012022,13061989],  # default if nothing is provided
)
CLI.add_argument(
  "--carbon_num",  # name on the CLI - drop the `--` for positional/required parameters
  nargs="*",  # 0 or more values expected => creates a list
  type=int,
  default=[3,6,12,18],  # default if nothing is provided
)
CLI.add_argument(
  "--microbial_group_num",  # name on the CLI - drop the `--` for positional/required parameters
  nargs="*",  # 0 or more values expected => creates a list
  type=int,
  default=[4,8,12,16,20,24,28,32],  # default if nothing is provided
)
user_args=CLI.parse_args()
cn_list = user_args.carbon_num#[3,6,12,18]
bio_n_series = user_args.microbial_group_num#[4,8,12,16,20,24,28,32]
seed_sim_list = user_args.seeds_num
filestring = user_args.sim_label + '_carbon_'

ip = 0 #0 random scenarios
init_dom_list = [1000,2000,5000,10000,15000]

def run_sims (experiment, c_n, b_n, dom_initial, seed_sim, Switch_matrix, hw):
    np.random.seed(seed_sim)
    sim = experiment + "/bio_n_"+ str(b_n) + "/dom_initial_" + str(dom_initial) + "/seed_" + str(seed_sim)

    # declare a time vector (time window)
    t_span = [0,36500]
    t_step = 5
    t_span_list = np.arange(t_span[0], t_span [1],t_step)
    total_dom_initial = dom_initial
    dom_bio_ratio_initial = 10
    mean_dom_initial = 1000
    mean_bio_initial = 100

    dom_n = c_n
    bio_n = b_n
    # Initialize the same number of parameters and initial conditions:
    dom_initial, biomass_initial = generate_random_initial_conditions (dom_n, bio_n, mean_dom_initial, mean_bio_initial, total_dom_initial, dom_bio_ratio_initial)
    x0 = np.append(dom_initial, biomass_initial)
    mu_enz = 0.4
    mu_z = 0.2
    mu_v = 0.004
    mu_k = 1
    skewness_para = 1.5 #sigma/standard deviation
    # First order rate constant for the production of exoenzymes by each microbial group
    enzparams= abs(np.log(np.random.lognormal(mu_enz, sigma = skewness_para*mu_enz, size=bio_n)))
    #Fraction of depolymerized carbon pool that is used for microbial uptake (respiration + growth).
    zparams= abs(np.log(np.random.lognormal(mu_z, sigma = skewness_para*mu_z, size=bio_n*dom_n)))
    # M-M max rate constant for consumption of carbon compound by a particular microbial group
    vparams= abs(np.log(np.random.lognormal(mu_v, sigma = skewness_para*mu_v,size=bio_n*dom_n)))
    # M-M half saturation constant for consumption of carbon compound by a particular microbial group
    kparams = 600*abs(np.log(np.random.lognormal(mu_k, sigma = skewness_para*mu_k,size=dom_n*bio_n)))
    # Second order rate constant for quadratic mortality rate of microbial groups
    mparams = np.median(vparams.reshape(dom_n, bio_n),axis=0)/(5*np.sum(biomass_initial)) #earlier this factor was 5 in gen_spec_skew
    ##the next two lines are specifically to constrain the initialised dom and NOSC profile.
    ## In the unconstrained NOSC simulations, the next 2 lines are not included in the run scripts.
    np.random.seed(seed_sim)
    ox_state_1 = np.random.uniform(-0.5, 0.5, dom_n-1)
    ox_state = np.insert(ox_state_1, np.random.choice(len(ox_state_1), size=1), -0.2)
    ##specific lines end here
    
    rel_tol_arr = 10**-6
    abs_tol_arr = 10**-6

    carbon_input = generate_random_boundary_conditions(dom_n, 0, method_name = "user_defined")

    trial = rn(maximum_capacity=0.15*total_dom_initial,carbon_num = dom_n,bio_num = bio_n, carbon_input = carbon_input, necromass_distribution="notequal", competition = "True")
    
    trial.set_rate_constants(ox_state, enzparams, zparams, vparams, kparams, mparams)
    trial.rearrange_constants()
    trial.identify_components_natures(recalcitrance_criterion="oxidation_state")
    trial.reorder_constants_with_comp_nature()
    trial.adaptation(dom_initial)
    trial.microbe_carbon_switch(Switch_matrix)
    solution = trial.solve_network(x0, t_span, t_span_list, solv_method = 'LSODA', rel_tol = rel_tol_arr, abs_tol = abs_tol_arr, first_tim_step = 0.01, max_tim_step = 10)

    seed_dic = {sim : {'dom_number': dom_n, 'biomass_number': bio_n,
    'oxidation_state': ox_state,
    'enzyme_production_parameters': trial.v_enz,
    'uptake_parameters': trial.z,
    'max_rate_parameters': trial.v_params,
    'sat_const_parameters': trial.k_params,
    'efficiency_parameters': trial.y_params,
    'mortality_parameters': trial.m_params,
    'initial_conditions_dom': dom_initial,
    'initial_conditions_biomass': biomass_initial,
    'carbon_input_boundary': carbon_input,
    'sim_status':solution.status,
    'message': solution.message}}
        
    sim_array = solution.y.T

    dataset_category = sim
    
    dataset_name = dataset_category + "/species_number"
    hw.create_dataset(dataset_name, data=[dom_n, bio_n])
    
    dataset_name = dataset_category + "/parameters/oxidation_state"
    hw.create_dataset(dataset_name, data=ox_state)

    dataset_name = dataset_category + "/parameters/exo_enzyme_rate"
    hw.create_dataset(dataset_name, data=trial.v_enz)

    dataset_name = dataset_category + "/parameters/carbon_uptake"
    hw.create_dataset(dataset_name, data=trial.z)

    dataset_name = dataset_category + "/parameters/max_rate"
    hw.create_dataset(dataset_name, data=trial.v_params)
    
    dataset_name = dataset_category + "/parameters/half_saturation"
    hw.create_dataset(dataset_name, data=trial.k_params)
    
    dataset_name = dataset_category + "/parameters/yield_coefficient"
    hw.create_dataset(dataset_name, data=trial.y_params)

    dataset_name = dataset_category + "/parameters/mortality_const"
    hw.create_dataset(dataset_name, data=trial.m_params)
    
    dataset_name = dataset_category + "/initial_conditions/dom"
    hw.create_dataset(dataset_name, data=dom_initial)
    
    dataset_name = dataset_category + "/initial_conditions/biomass"
    hw.create_dataset(dataset_name, data=biomass_initial)

    dataset_name = dataset_category + "/solution/dom"
    hw.create_dataset(dataset_name, data=sim_array[:,:dom_n])
    
    dataset_name = dataset_category + "/solution/biomass"
    hw.create_dataset(dataset_name, data=sim_array[:,dom_n:])
    
    dataset_name = dataset_category + "/solution/status"
    hw.create_dataset(dataset_name, data = np.asarray(solution.status))

    return sim, seed_dic

def run_sim (random_seed_number, c_n):
    details_subfolder = filestring + str(c_n) + '_' + str(random_seed_number)+'_ip_' + str(ip)
    simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
    results_dir = os.path.join(project_dir, "results", details_subfolder)
    figures_dir = os.path.join(project_dir, "figures", details_subfolder)

    for sub_dir in [simulations_dir, results_dir, figures_dir]:
        if os.path.exists(sub_dir)=="True":
            print("Path exists already")
            break
        else:
            os.makedirs(sub_dir)

    hfile_to_write = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'w')

    empty_dic = {}
    failed_count = 0

    for N in bio_n_series:
        baseline = "b_1"
        label = "all"
        S_witch_b_1 = 1
        for t_dom_initial in init_dom_list:
                exp_details = baseline + "_" + label + "_"
                sim_id, seed_dictionary = run_sims (exp_details, c_n, N, t_dom_initial, random_seed_number, S_witch_b_1, hfile_to_write)
                if seed_dictionary[sim_id]["sim_status"] < 0:
                    failed_count +=1
                empty_dic.update(seed_dictionary)

    rng = np.random.default_rng()
    for N in bio_n_series:
        S_witches_4 = np.zeros((4,c_n,N))
        for n, baseline, activity_pc in zip([0,1,2,3],["b_2", "b_3", "b_4", "b_5"] ,[0.1, 0.25, 0.5, 0.75]):
            K = math.ceil(activity_pc*N*c_n)
            s = random.sample(range(0,N*c_n), K)
            s_r = list(int(np.floor(x/N)) for x in s)
            s_c = list( x - y*N for x,y in zip(s, s_r))
            S_witches_4[n,s_r,s_c] = 1
        for n, baseline, activity_pc in zip([0,1,2,3],["b_2", "b_3", "b_4", "b_5"] ,[0.1, 0.25, 0.5, 0.75]):
            for label, random_seed in zip(["a", "b", "c", "d", "e"], [1,2,3,4,5]):
                x_s_witch = rng.permutation(S_witches_4[n,:,:], axis = 1)
                for t_dom_initial in init_dom_list:
                    exp_details = baseline + "_" + label + "_"
                    sim_id, seed_dictionary = run_sims (exp_details, c_n, N, t_dom_initial, random_seed_number, x_s_witch, hfile_to_write)
                    if seed_dictionary[sim_id]["sim_status"] < 0:
                        failed_count +=1
                    empty_dic.update(seed_dictionary)

    # Save all results in HDF5 file
    hfile_to_write.close()

    print ("All results saved in file: ", hfile_to_write)
    print ("All simulations have ended. Failed simulations totaled: ", failed_count)
    pickle_file = os.path.join(simulations_dir,"seeds_randoms.pkl")
    # create a binary pickle file 
    f = open(pickle_file,"wb")
    # write the python object (dict) to pickle file
    pickle.dump(empty_dic,f)
    # close file
    f.close()
    print ("All seeds details saved in file: ", pickle_file)

    csv_file = os.path.join(simulations_dir,"seeds_randoms.csv")
    # open file for writing, "w" is writing
    w = csv.writer(open(csv_file, "w"))
    # loop over dictionary keys and values
    for key, val in empty_dic.items():
        # write every key and value to file
        w.writerow([key, val])
    print ("All seeds details saved in file: ", pickle_file)

    return None

for carbon_num in cn_list:
    for seed_sim in seed_sim_list:
        run_sim (seed_sim, carbon_num)
        print ("Completed simulations for seed ", seed_sim)