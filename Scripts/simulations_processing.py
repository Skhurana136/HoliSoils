## Import libraries
import os
import pandas as pd
import numpy as np
import h5py

from DS.analyses.steady_state import diversity_carbon, normalize_carbon

## LOAD RESULTS
project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
results_dir = os.path.join(project_dir, "simulations")
output_dir = os.path.join(project_dir, "results")

hr = h5py.File(os.path.join(results_dir,"simulations.h5"), mode = 'r')

# Load all datasets and save their Shannon and diversity indices in a dataframe:

n = len(hr)

emptydf = pd.DataFrame (columns = ["Sim", "carbon_species", "biomass_species", "Shannon", "DOC", "TOC"])
emptydf_norm = pd.DataFrame (columns = ["Sim", "carbon_species", "biomass_species", "Shannon", "DOC", "TOC"])

for sim in list(range(n)):
    sim_data = hr[str(sim)]
    # Identify number of dom and microbial species:
    dom_n, bio_n = sim_data['species_number']
    x = sim_data['solution']
    C = x['dom']
    B = x['biomass']
    data = np.append(C,B,axis=1)
    S, DOC, TOC = diversity_carbon(data[-2:,:], dom_n, bio_n)
    S_series = pd.Series(S)
    DOC_series = pd.Series(DOC)
    TOC_series =pd.Series(TOC)
    sim_series = pd.Series(np.zeros(len(S))+sim)
    dom_series = pd.Series(np.zeros(len(S))+dom_n)
    bio_series = pd.Series(np.zeros(len(S))+bio_n)
    sim_df = pd.DataFrame ({"Sim":sim_series, "carbon_species":dom_series, "biomass_species":bio_series,
    "Shannon":S_series, "DOC":DOC_series, "TOC":TOC_series})
    emptydf = emptydf.append(sim_df, ignore_index=True)

    initial_data = sim_data['initial_conditions']
    carbon_initial = initial_data['dom']
    biomass_initial = initial_data['biomass']

    C_norm, B = normalize_carbon(data, dom_n, bio_n, carbon_initial, biomass_initial)
    data_norm = np.append(C_norm[:,None],B,axis=1)
    S, DOC, TOC = diversity_carbon(data_norm[-2:,:], dom_n, bio_n)
    S_series = pd.Series(S)
    DOC_series = pd.Series(DOC)
    TOC_series =pd.Series(TOC)
    sim_df = pd.DataFrame ({"Sim":sim_series, "carbon_species":dom_series, "biomass_species":bio_series,
    "Shannon":S_series, "DOC":DOC_series, "TOC":TOC_series})
    emptydf_norm = emptydf_norm.append(sim_df, ignore_index=True)


print("The shape of the merged dataframe is ", emptydf.shape)
print("The dataframe contains the following data types ", emptydf.dtypes)

filename = os.path.join(output_dir, "diversity_data.pkl")
emptydf.to_pickle(filename)
print ("Diversity with carbon data is saved here ", filename)

filename = os.path.join(output_dir, "diversity_normalized_data.pkl")
emptydf_norm.to_pickle(filename)
print ("Diversity with normalized carbon data is saved here ", filename)