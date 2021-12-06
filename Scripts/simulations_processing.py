## Import libraries
import os
import pandas as pd
import numpy as np
import h5py

from DS.solvers.diff_eqn_system import diversity_carbon

## LOAD RESULTS
project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
results_dir = os.path.join(project_dir, "simulations")
hw = h5py.File(os.path.join(results_dir,"simulations.h5"), mode = 'r')

# Load all datasets and save their Shannone and diversity indices in a dataframe:

n = len(hw)

emptydf = pd.DataFrame (columns = ["Sim", "carbon_species", "biomass_species", "Shannon", "DOC", "TOC"])

for sim in list(range(n)):
    data = hw[str(sim)]
    # Identify number of dom and microbial species:
    dom_n, bio_n = data['species_number']
    x = data['solution']
    C = x['dom']
    B = x['biomass']
    data = np.append(C,B,axis=1)
    S, DOC, TOC = diversity_carbon(data, dom_n, bio_n)
    S_series = pd.Series(S)
    DOC_series = pd.Series(DOC)
    TOC_series =pd.Series(TOC)
    sim_series = pd.Series(np.zeros(len(S)+sim))
    dom_series = pd.Series(np.zeros(len(S)+dom_n))
    bio_series = pd.Series(np.zeros(len(S)+bio_n))
    sim_df = pd.DataFrame ([sim_series, dom_series, bio_series, S_series, DOC_series, TOC_series])
    emptydf = emptydf.append(sim_df, ignore_index=True)

print("The shape of the merged dataframe is ", emptydf.shape)
print("The dataframe contains the following data types ", emptydf.dtypes)