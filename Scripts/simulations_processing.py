# ## Import libraries
import os
import pandas as pd
import numpy as np
import h5py

from DS.analyses.steady_state import diversity_carbon, normalize_carbon

## LOAD RESULTS
#project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
project_dir = "C:/Users/swkh9804/Documents/Projects/HoliSoils/data"
details_subfolder = 'from_spyder'
simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
results_dir = os.path.join(project_dir, "results", details_subfolder)
figures_dir = os.path.join(project_dir, "figures", details_subfolder)


hr = h5py.File(os.path.join(results_dir,"simulations.h5"), mode = 'r')

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
        proportion = B/B.sum(axis=1, keepdims = True)#np.sum(B,axis=0)
        S = -np.sum(proportion*np.log(proportion), axis = 1)
        TOC = np.sum(C,axis=1) + np.sum(B, axis=1)
        DOC = np.sum(C,axis=1)
        DOC_max = np.max(DOC)
        S_DOC_max = S[np.argmax(DOC)]
        S_max = np.max(S)
        DOC_S_max = DOC[np.argmax(S)]
        DOC_end = DOC[-1]
        #S, DOC, TOC = diversity_carbon(data[-2:,:], dom_n, bio_n)
        #S_series = pd.Series(S)
        #DOC_series = pd.Series(DOC)
        #TOC_series =pd.Series(TOC)
        #sim_series = pd.Series(np.zeros(len(S))+sim)
        #dom_series = pd.Series(np.zeros(len(S))+dom_n)
        #bio_series = pd.Series(np.zeros(len(S))+bio_n)
        #sim_df = pd.DataFrame ({"Sim":sim_series, "carbon_species":dom_series, "biomass_species":bio_series,
        #"Shannon":S_series, "DOC":DOC_series, "TOC":TOC_series})
        #emptydf = emptydf.append(sim_df, ignore_index=True)

        initial_data = sim_data['initial_conditions']
        carbon_initial = np.asarray(initial_data['dom'])
        biomass_initial = np.asarray(initial_data['biomass'])
        proportion_i = biomass_initial/np.sum(biomass_initial)
        
        S_i = -np.sum(proportion_i*np.log(proportion_i))
        DOC_i = np.sum(carbon_initial)

        #C_norm, B = normalize_carbon(data, dom_n, bio_n, carbon_initial, biomass_initial)
        #data_norm = np.append(C_norm[:,None],B,axis=1)
        #S, DOC, TOC = diversity_carbon(data_norm[-2:,:], dom_n, bio_n)
        #S_series = pd.Series(S)
        #DOC_series = pd.Series(DOC)
        #TOC_series =pd.Series(TOC)
        #sim_df = pd.DataFrame ({"Sim":sim_series, "carbon_species":dom_series, "biomass_species":bio_series,
        #"Shannon":S_series, "DOC":DOC_series, "TOC":TOC_series})
        #emptydf_norm = emptydf_norm.append(sim_df, ignore_index=True)
        row.append([sim, DOC_i, S_i, DOC_max, S_DOC_max, S_max, DOC_S_max, DOC_end])

diversity_data = pd.DataFrame.from_records(row, columns = ["Sim_series", "DOC_initial", "S_initial", "Max_DOC", "S_at_max_DOC", "S_max", "DOC_at_max_S", "DOC_end"])

#print("The shape of the merged dataframe is ", emptydf.shape)
#print("The dataframe contains the following data types ", emptydf.dtypes)

filename = os.path.join(results_dir, "diversity_data.pkl")
diversity_data.to_pickle(filename)
print ("Diversity with carbon data is saved here ", filename)

#filename = os.path.join(results_dir, "diversity_normalized_data.pkl")
#emptydf_norm.to_pickle(filename)
#print ("Diversity with normalized carbon data is saved here ", filename)

hr.close()