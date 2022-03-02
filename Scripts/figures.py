## Import libraries
import os
import pandas as pd
import numpy as np
import h5py

from DS.analyses.steady_state import diversity_carbon

import matplotlib.pyplot as plt
import matplotlib.lines as mlines

## LOAD RESULTS
#project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
project_dir = "C:/Users/swkh9804/Documents/Projects/HoliSoils/data"
results_dir = os.path.join(project_dir, "simulations")
output_dir = os.path.join(project_dir, "results")
figures_dir = os.path.join(project_dir, "figures")

hr = h5py.File(os.path.join(results_dir,"simulations.h5"), mode = 'r')

# Load all datasets and save their Shannon and diversity indices in a dataframe:

n = len(hr)

for sim in list(range(n)):
    sim_data = hr[str(sim)]
    # Identify number of dom and microbial species:
    dom_n, bio_n = sim_data['species_number']
    x = sim_data['solution']
    C = x['dom']
    B = x['biomass']
    data = np.append(C,B,axis=1)
    S, DOC, TOC = diversity_carbon(data, dom_n, bio_n)
    S_series = pd.Series(S)
    DOC_series = pd.Series(DOC)
    TOC_series =pd.Series(TOC)
    sim_series = pd.Series(np.zeros(len(S))+sim)
    dom_series = pd.Series(np.zeros(len(S))+dom_n)
    bio_series = pd.Series(np.zeros(len(S))+bio_n)

    # Figures
    # Time series of biomass and DOC concentration 
    plt.figure()
    plt.plot (B, linestyle = '--')#, label = "Bacteria")
    plt.plot (C, linestyle = '-')#, label = "Bacteria")
    plt.xlabel ("Time (days)")
    plt.ylabel ("Concentration")
    plt.legend()
    plt.text(s = "Seed: "+str(sim), x = 80, y = 800, fontsize = 12)
    solid_line = mlines.Line2D([], [], linestyle = '-', color='grey', linewidth=1.5, label='DOM')
    dashed_line = mlines.Line2D([], [], linestyle = '--', color='grey', linewidth=1.5, label='Biomass')
    first_legend = plt.legend(handles=[solid_line, dashed_line])
    plt.savefig(os.path.join(figures_dir, str(sim) + "_Time_series_concentration.png"), dpi = 300,
    layout = 'tight')

    # Time series of biomass and DOC concentration 
    fig, axes = plt.subplots(nrows=3, ncols=1, sharex = True, figsize = (6,8))
    axes.flat[0].plot(S)
    axes.flat[1].plot(DOC)
    axes.flat[2].plot(TOC)
    axes.flat[0].set_ylabel("Shannon diversity")
    axes.flat[1].set_ylabel("DOC")
    axes.flat[2].set_ylabel("TOC")
    axes.flat[2].set_xlabel("Time")
    axes.flat[0].text(s = "Seed: "+str(sim), x = 0.8, y = 0.8, transform = axes.flat[0].transAxes,
    fontsize = 12)
    plt.savefig(os.path.join(figures_dir, str(sim) + "_diversity_concentration.png"), dpi = 300,
    layout = 'tight')

hr.close()
                                        
filename = os.path.join(output_dir, "diversity_data.pkl")
diversity_data = pd.read_pickle(filename)

plt.figure()
plt.scatter(x = 'Shannon',y='DOC', data = diversity_data, s = 'carbon_species')
plt.xlabel ('Shannon diversity index')
plt.ylabel("DOC")