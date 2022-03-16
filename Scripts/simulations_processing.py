# ## Import libraries
import os
import pandas as pd
import numpy as np
import h5py

from DS.analyses.steady_state import diversity_carbon, normalize_carbon

## LOAD RESULTS
#project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
project_dir = "C:/Users/swkh9804/Documents/Projects/HoliSoils/data"
details_subfolder = 'n_1000_160322'
simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
results_dir = os.path.join(project_dir, "results", details_subfolder)
figures_dir = os.path.join(project_dir, "figures", details_subfolder)


hr = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'r')

# Load all datasets and save their Shannon and diversity indices in a dataframe:

n = len(hr)

#emptydf = pd.DataFrame (columns = ["Sim", "carbon_species", "biomass_species", "Shannon", "DOC", "TOC"])
#emptydf_norm = pd.DataFrame (columns = ["Sim", "carbon_species", "biomass_species", "Shannon", "DOC", "TOC"])
row = []
for sim in list(range(n)):
    sim_data = hr[str(sim)]
    # Identify number of dom and microbial species:
    dom_n, bio_n = sim_data['species_number']
    x = sim_data['solution']
    C = np.asarray(x['dom'])
    B = np.asarray(x['biomass'])
    data = np.append(C,B,axis=1)
    neg_ind = np.where(data<0)
    if np.size(neg_ind)>0:
        pass        
    else:
        ## SHANNON DIVERSITY AND CARBON STOCK
        DOC = np.sum(C,axis=1)
        if len(DOC)<100000: #check for complete simulation
            pass
        else:
            proportion = B/B.sum(axis=1, keepdims = True)
            S = -np.sum(proportion*np.log(proportion), axis = 1)
            TOC = np.sum(C,axis=1) + np.sum(B, axis=1)
            DOC_end = DOC[-1]
            Biomass = np.sum(B, axis = 1)
            bio_max = np.max(Biomass)
            DOC_max = np.max(DOC)
            S_DOC_max = S[np.argmax(DOC)]
            S_max = np.max(S)
            DOC_S_max = DOC[np.argmax(S)]
            S_end = S[-1]
            initial_data = sim_data['initial_conditions']
            carbon_initial = np.asarray(initial_data['dom'])
            biomass_initial = np.asarray(initial_data['biomass'])
            proportion_i = biomass_initial/np.sum(biomass_initial)
            S_i = -np.sum(proportion_i*np.log(proportion_i))
            DOC_i = np.sum(carbon_initial)
            bio_i = np.sum(biomass_initial)
            bio_end = Biomass[-1]
            row.append([sim, dom_n, bio_n, DOC_i, S_i, S_end, DOC_max, S_DOC_max, S_max, DOC_S_max, DOC_end, bio_i, bio_end, bio_max])

diversity_data = pd.DataFrame.from_records(row, columns = ["Sim_series", "carbon_species", "biomass_species", "DOC_initial", "S_initial", "S_end","Max_DOC", "S_at_max_DOC", "S_max", "DOC_at_max_S", "DOC_end", "Biomass_initial", "Biomass_end", "Biomass_max"])

#print("The shape of the merged dataframe is ", emptydf.shape)
#print("The dataframe contains the following data types ", emptydf.dtypes)

filename = os.path.join(results_dir, "diversity_data.pkl")
diversity_data.to_pickle(filename)
print ("Diversity with carbon data is saved here ", filename)

#filename = os.path.join(results_dir, "diversity_normalized_data.pkl")
#emptydf_norm.to_pickle(filename)
#print ("Diversity with normalized carbon data is saved here ", filename)

hr.close()