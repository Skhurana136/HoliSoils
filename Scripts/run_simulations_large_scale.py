## Import libraries
import numpy as np
import h5py
import os

from DS.solvers.diff_eqn_system import ReactionNetwork as rn
from DS.solvers.diff_eqn_system import generate_random_parameters
from DS.solvers.diff_eqn_system import generate_random_initial_conditions
from DS.solvers.diff_eqn_system import generate_random_boundary_conditions

#project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
project_dir = "C:/Users/swkh9804/Documents/Projects/HoliSoils/data"
results_dir = os.path.join(project_dir, "simulations")

hw = h5py.File(os.path.join(results_dir,"simulations.h5"), mode = 'w')
# Run 1000 random simulations
n=1000
# declare a time vector (time window)
t_span = [0,1000]
t = np.arange(t_span[0], t_span [1],0.01)
for sim in list(range(n)):
    # Set seed
    np.random.seed(sim)

    # Number of DOM/Carbon species:

    dom_n = np.random.randint(5,50,1)[0]
    bio_n = np.random.randint(4,40,1)[0]

    # Initialize the same number of parameters and initial conditions:
    enzparams, vparams, kparams, yparams, mparams = generate_random_parameters(dom_n, bio_n,5)
    
    dom_initial, biomass_initial = generate_random_initial_conditions (dom_n, bio_n)
    
    x0 = np.append(dom_initial, biomass_initial)
    
    carbon_input = generate_random_boundary_conditions()

    seed_dic = {sim : {'dom_number': dom_n, 'biomass_number': bio_n,
    'enzyme_production_parameters': enzparams,
    'max_rate_parameters': vparams,
    'sat_const_parameters': kparams,
    'efficiency_parameters': yparams,
    'mortality_parameters': mparams,
    'initial_conditions_dom': dom_initial,
    'initial_conditions_biomass': biomass_initial,
    'carbon_input_boundary': carbon_input}}

    trial = rn(maximum_capacity=5,
        #carbon_mol_bio = 10,
        carbon_num = dom_n,
        bio_num = bio_n,
        carbon_input = carbon_input,
        sigmoid_coeff_stolpovsky = 0.01,
        enzyme_production_rate_constant = enzparams[0],
        efficiency_bio_uptake = enzparams[1],
        necromass_distribution="equal")  
    
    trial.set_rate_constants(vparams, kparams, yparams,mparams)
    
    trial.identify_components_natures()
    
    solution = trial.solve_network(x0, t_span, t)

    tim = solution.t
    sim_array = solution.y.T

    dataset_category = str(sim)

    dataset_name = dataset_category + "/species_number"
    hw.create_dataset(dataset_name, data=[dom_n, bio_n])
    
    dataset_name = dataset_category + "/parameters/max_rate"
    hw.create_dataset(dataset_name, data=vparams.reshape(bio_n, dom_n))
    
    dataset_name = dataset_category + "/parameters/half_saturation"
    hw.create_dataset(dataset_name, data=kparams.reshape(bio_n, dom_n))
    
    dataset_name = dataset_category + "/parameters/yield_coefficient"
    hw.create_dataset(dataset_name, data=yparams)
    
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