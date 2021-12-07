## Import libraries
import numpy as np
import h5py
import os

from DS.solvers.diff_eqn_system import ReactionNetwork as rn

project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
results_dir = os.path.join(project_dir, "simulations")

hw = h5py.File(os.path.join(results_dir,"simulations.h5"), mode = 'w')
# Run 1000 random simulations
n=1000
# declare a time vector (time window)
t = np.arange(0,100,0.01)
for sim in list(range(n)):
    # Set seed
    np.random.seed(sim)

    # Number of DOM/Carbon species:

    dom_n = np.random.randint(4,10,1)[0]
    bio_n = np.random.randint(2,10,1)[0]

    # Initialize the same number of parameters:

    vparams = np.random.uniform(0.2,10,bio_n*dom_n)
    kparams = np.random.randint(1000,8000,bio_n*dom_n)
    yparams = np.random.uniform(0.1,0.99,bio_n)
    mparams = np.mean(vparams.reshape(dom_n, bio_n),axis=0)/10

    # initial conditions
    dom_initial = np.random.randint(500,1000,dom_n)
    biomass_initial = np.random.randint(100,1000,bio_n)
    x0 = np.append(dom_initial,biomass_initial)

    trial = rn(maximum_capacity=5000,
            #carbon_mol_bio = 10,
            carbon_num = dom_n,
            bio_num = bio_n,
            carbon_input = np.max(mparams),
            sigmoid_coeff_stolpovsky = 0.1)

    trial.set_rate_constants(vparams, kparams, yparams,mparams)
    sim_array = trial.solve_network(x0, t)
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