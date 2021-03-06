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
project_dir = "C:/Users/swkh9804/Documents/Projects/HoliSoils/data"

seed_sim = 13061989
#np.random.seed(seed_sim)

details_subfolder = 'carbon_3_ip_0_20k_days'
c_n = 3
bio_n_series = [2,3,5,6]

def run_sims (experiment, c_n, b_n, dom_initial, seed_sim, Switch_matrix, hw):
    np.random.seed(seed_sim)
    sim = experiment + "/bio_n_"+ str(N) + "/dom_initial_" + str(dom_initial) + "/seed_" + str(seed_sim)

    # declare a time vector (time window)
    t_span = [0,20000]
    t_span_list = np.arange(t_span[0], t_span [1],0.1)
    total_dom_initial = dom_initial
    dom_bio_ratio_initial = 10
    mean_dom_initial = 1000
    mean_bio_initial = 100

    dom_n = c_n
    bio_n = b_n
    # Initialize the same number of parameters and initial conditions:
    initial_conditions_file = os.path.join(project_dir, "simulations",details_subfolder[:-9], "simulations.h5")
    prev_sim_hfile = h5py.File(initial_conditions_file, mode = 'r')
    prev_sim_data = prev_sim_hfile[experiment]["bio_n_"+str(b_n)]['dom_initial_'+str(dom_initial)]['seed_'+str(seed_sim)]
    dom_initial = prev_sim_data['solution']['dom'][-1,:]
    biomass_initial = prev_sim_data['solution']['biomass'][-1,:]
    print(np.shape(dom_initial), np.shape(biomass_initial))
    params = prev_sim_data['parameters']
    ox_state = params['oxidation_state'][:]
    enzparams = params['exo_enzyme_rate'][:]
    zparams = params['carbon_uptake'][:]
    vparams = params['max_rate'][:]
    kparams = params['half_saturation'][:]
    mparams = params['mortality_const'][:]
    yparams = params['yield_coefficient'][:]

    x0 = np.append(dom_initial, biomass_initial)
    
    carbon_input = generate_random_boundary_conditions(dom_n, 0, method_name = "user_defined")#1*total_dom_initial/100

    trial = rn(maximum_capacity=5,carbon_num = dom_n,bio_num = bio_n, carbon_input = carbon_input, necromass_distribution="notequal")
    
    trial.set_rate_constants(yparams, ox_state, enzparams, zparams, vparams, kparams, mparams) 
    trial.constants_previously_defined()
    trial.identify_components_natures(recalcitrance_criterion='oxidation_state')
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
    
    return sim, seed_dic

def shuffle_switch_matrix (matrix_to_switch, shuffle_seed):
    np.random.seed(shuffle_seed)
    rng.permutation(matrix_to_switch, axis = 1)
    return matrix_to_switch

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
    for t_dom_initial in [1000,5000,10000,15000,20000]:
            exp_details = baseline + "_" + label + "_"
            sim_id, seed_dictionary = run_sims (exp_details, c_n, N, t_dom_initial, seed_sim, S_witch_b_1, hfile_to_write)
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
            for t_dom_initial in [1000,5000,10000,15000,20000]:
                exp_details = baseline + "_" + label + "_"
                sim_id, seed_dictionary = run_sims (exp_details, c_n, N, t_dom_initial, seed_sim, x_s_witch, hfile_to_write)
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