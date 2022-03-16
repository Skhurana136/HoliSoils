
# ## Import libraries
import numpy as np
import h5py
import os
import pickle
import csv

from DS.solvers.diff_eqn_system import ReactionNetwork as rn
from DS.solvers.diff_eqn_system import generate_random_parameters
from DS.solvers.diff_eqn_system import generate_random_initial_conditions
from DS.solvers.diff_eqn_system import generate_random_boundary_conditions

#project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
project_dir = "C:/Users/swkh9804/Documents/Projects/HoliSoils/data"

details_subfolder = 'n_1000_160322'
simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
results_dir = os.path.join(project_dir, "results", details_subfolder)
figures_dir = os.path.join(project_dir, "figures", details_subfolder)

for sub_dir in [simulations_dir, results_dir, figures_dir]:
    if os.path.exists(sub_dir)=="True":
        print("Path exists already")
        break
    else:
        os.mkdir(sub_dir)

hw = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'w')
# Run 100 random simulations
n=1000
# declare a time vector (time window)
t_span = [0,1000]
t = np.arange(t_span[0], t_span [1],0.01)

empty_dic = {}

for sim in list(range(n)):
    # Set seed
    np.random.seed(sim)

    # Number of DOM/Carbon species:

    dom_n = np.random.randint(5,50,1)[0]
    bio_n = np.random.randint(4,40,1)[0]

    # Initialize the same number of parameters and initial conditions:
    dom_initial, biomass_initial = generate_random_initial_conditions (dom_n, bio_n, 1000, 80)
    
    ox_state, enzparams, zparams, vparams, kparams, yparams, mparams = generate_random_parameters(dom_n, bio_n,5*biomass_initial)
    
    x0 = np.append(dom_initial, biomass_initial)
    
    carbon_input = generate_random_boundary_conditions(dom_n, 50, method_name = "user_defined")

    trial = rn(maximum_capacity=5,carbon_num = dom_n,bio_num = bio_n, carbon_input = carbon_input, necromass_distribution="notequal")
    
    trial.set_rate_constants(ox_state, enzparams, zparams, vparams, kparams, yparams,mparams)
    
    trial.identify_components_natures(recalcitrance_criterion="oxidation_state")

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
    empty_dic.update(seed_dic)

    solution = trial.solve_network(x0, t_span, t)

    tim = solution.t
    sim_array = solution.y.T

    dataset_category = str(sim)

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
    
    dataset_name = dataset_category + "/initial_conditions/dom"
    hw.create_dataset(dataset_name, data=dom_initial)
    
    dataset_name = dataset_category + "/initial_conditions/biomass"
    hw.create_dataset(dataset_name, data=biomass_initial)
    
    dataset_name = dataset_category + "/solution/dom"
    hw.create_dataset(dataset_name, data=sim_array[:,:dom_n])
    
    dataset_name = dataset_category + "/solution/biomass"
    hw.create_dataset(dataset_name, data=sim_array[:,dom_n:])


# Save all results in HDF5 file
hw.close()
print ("All results saved in file: ", hw)

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