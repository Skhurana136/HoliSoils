## Import libraries
import os
import pandas as pd
import numpy as np
import h5py

from DS.analyses.steady_state import diversity_carbon

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns

## LOAD RESULTS
#project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
project_dir = "C:/Users/swkh9804/Documents/Projects/HoliSoils/data"
details_subfolder = 'paras_adjust_v2'
simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
results_dir = os.path.join(project_dir, "results", details_subfolder)
figures_dir = os.path.join(project_dir, "figures", details_subfolder)

hr = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'r')

# Load all datasets and save their Shannon and diversity indices in a dataframe:

sim_list = list(hr.keys())

for sim in sim_list:
    sim_data = hr[str(sim)]
    # Identify number of dom and microbial species:
    dom_n, bio_n = sim_data['species_number']
    x = sim_data['solution']
    C = np.asarray(x['dom'])
    B = np.asarray(x['biomass'])
    data = np.append(C,B,axis=1)

    proportion = B/B.sum(axis=1, keepdims = True)
    S = -np.sum(proportion*np.log(proportion), axis = 1)
    TOC = np.sum(C,axis=1) + np.sum(B, axis=1)
    DOC = np.sum(C,axis=1)
    DOC_max = np.max(DOC)
    S_DOC_max = S[np.argmax(DOC)]
    S_max = np.max(S)
    DOC_S_max = DOC[np.argmax(S)]
    DOC_end = DOC[-1]

    initial_data = sim_data['initial_conditions']
    carbon_initial = np.asarray(initial_data['dom'])
    biomass_initial = np.asarray(initial_data['biomass'])
    proportion_i = biomass_initial/np.sum(biomass_initial)
       
    S_i = -np.sum(proportion_i*np.log(proportion_i))
    DOC_i = np.sum(carbon_initial)

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
    plt.savefig(os.path.join(figures_dir, str(sim) + "_Time_series_concentration.png"), dpi = 300)

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
    plt.savefig(os.path.join(figures_dir, str(sim) + "_diversity_concentration.png"), dpi = 300)

hr.close()