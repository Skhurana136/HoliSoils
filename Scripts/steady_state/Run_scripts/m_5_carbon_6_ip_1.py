# ## Import libraries
from random import seed
import numpy as np
import h5py
import os
import pickle
import csv
import math

from DS.solvers.diff_eqn_system import ReactionNetwork as rn
from DS.solvers.diff_eqn_system import generate_random_parameters
from DS.solvers.diff_eqn_system import generate_random_initial_conditions
from DS.solvers.diff_eqn_system import generate_random_boundary_conditions

#project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
project_dir = "D:\Projects\HoliSoils\data"

seed_sim_list = [420,13012022,13061989]
#np.random.seed(seed_sim)

c_n = 6
bio_n_series = [3,6,9,12]
ip = 1
init_dom_list = [1000,2000,5000,10000,15000]

def run_sims (experiment, c_n, b_n, dom_initial, seed_sim, Switch_matrix, hw):
    np.random.seed(seed_sim)
    sim = experiment + "/bio_n_"+ str(b_n) + "/dom_initial_" + str(dom_initial) + "/seed_" + str(seed_sim)

    # declare a time vector (time window)
    t_span = [0,50000]
    t_step = 0.1
    t_span_list = np.arange(t_span[0], t_span [1],t_step)
    total_dom_initial = dom_initial
    dom_bio_ratio_initial = 10
    mean_dom_initial = 1000
    mean_bio_initial = 100

    dom_n = c_n
    bio_n = b_n
    # Initialize the same number of parameters and initial conditions:
    dom_initial, biomass_initial = generate_random_initial_conditions (dom_n, bio_n, mean_dom_initial, mean_bio_initial, total_dom_initial, dom_bio_ratio_initial)
    
    ox_state, enzparams, zparams, vparams, kparams, mparams = generate_random_parameters(dom_n, bio_n,5*np.sum(biomass_initial))
    
    x0 = np.append(dom_initial, biomass_initial)
    
    carbon_input = generate_random_boundary_conditions(dom_n, total_dom_initial*t_step/(365*dom_n), method_name = "user_defined")

    trial = rn(maximum_capacity=5,carbon_num = dom_n,bio_num = bio_n, carbon_input = carbon_input, necromass_distribution="notequal")
    
    trial.set_rate_constants(ox_state, enzparams, zparams, vparams, kparams, mparams)
    trial.rearrange_constants()
    trial.identify_components_natures(recalcitrance_criterion="oxidation_state")
    trial.reorder_constants_with_comp_nature()
    trial.microbe_carbon_switch(Switch_matrix)
    solution = trial.solve_network(x0, t_span, t_span_list)

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
    'carbon_input_boundary': carbon_input}}

    tim = solution.t
    
    if tim.size<len(t_span_list):
        status = {'sim_status' : 'failed'}
        seed_dic[sim].update(status)
    else:
        status = {'sim_status' : 'complete'}
        seed_dic[sim].update(status)
        
    sim_array = solution.y.T.copy()[::int(5/t_step)]

    dataset_category = sim#dataset_category_1 + "/" + dataset_category_2
    
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
    
    return sim, seed_dic

def run_sim (random_seed_number):
    details_subfolder = 'carbon_'+str(c_n) + '_' + str(random_seed_number)+'_ip_' + str(ip)
    simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
    results_dir = os.path.join(project_dir, "results", details_subfolder)
    figures_dir = os.path.join(project_dir, "figures", details_subfolder)

    for sub_dir in [simulations_dir, results_dir, figures_dir]:
        if os.path.exists(sub_dir)=="True":
            print("Path exists already")
            break
        else:
            os.mkdir(sub_dir)

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
                if seed_dictionary[sim_id]["sim_status"] == "failed":
                    failed_count +=1
                empty_dic.update(seed_dictionary)

    rng = np.random.default_rng()
    for N in bio_n_series:
        S_witches_4 = np.zeros((4,c_n,N))
        for n, baseline, activity_pc in zip([0,1,2,3],["b_2", "b_3", "b_4", "b_5"] ,[0.1, 0.3, 0.5, 0.7]):
            K = math.ceil(activity_pc*N)
            s = np.array([1]*K + [0] * (N-K))
            S_witch_i = np.zeros((c_n, N))
            i=0
            while i<c_n:
                S_witch_i[i,:] = rng.permutation(s)
                i+=1
            S_witches_4[n,:,:] = S_witch_i
        for n, baseline, activity_pc in zip([0,1,2,3],["b_2", "b_3", "b_4", "b_5"] ,[0.1, 0.3, 0.5, 0.7]):
            for label, random_seed in zip(["a", "b", "c", "d", "e"], [1,2,3,4,5]):
                x_s_witch = rng.permutation(S_witches_4[n,:,:], axis = 1)
                for t_dom_initial in init_dom_list:
                    exp_details = baseline + "_" + label + "_"
                    sim_id, seed_dictionary = run_sims (exp_details, c_n, N, t_dom_initial, random_seed_number, x_s_witch, hfile_to_write)
                    if seed_dictionary[sim_id]["sim_status"] == "failed":
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

    # define a dictionary with key value pairs
    csv_file = os.path.join(simulations_dir,"seeds_randoms.csv")
    # open file for writing, "w" is writing
    w = csv.writer(open(csv_file, "w"))
    # loop over dictionary keys and values
    for key, val in empty_dic.items():
        # write every key and value to file
        w.writerow([key, val])
    print ("All seeds details saved in file: ", pickle_file)

    return None

for seed_sim in seed_sim_list:
    run_sim (seed_sim)
    print ("Completed simulations for seed ", seed_sim)