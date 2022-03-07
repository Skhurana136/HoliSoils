#%%
## Import libraries
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from DS.solvers.diff_eqn_system import ReactionNetwork as rn
from DS.solvers.diff_eqn_system import generate_random_parameters
from DS.solvers.diff_eqn_system import generate_random_initial_conditions
from DS.solvers.diff_eqn_system import generate_random_boundary_conditions
from DS.analyses.steady_state import diversity_carbon, normalize_carbon

print ("All libraries loaded.")

#%%
# Assign project directory depending on where you are:
# Work computer:
project_dir = "C:/Users/swkh9804/Documents/Projects/HoliSoils/data"

#%%
# Personal computer:
project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
#%%
# Assign child directories:
details_subfolder = 'varied_c_07032022'
simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
results_dir = os.path.join(project_dir, "results", details_subfolder)
figures_dir = os.path.join(project_dir, "figures", details_subfolder)

for sub_dir in [simulations_dir, results_dir, figures_dir]:
    if os.path.exists(sub_dir)=="True":
        print("Path exists already")
        break
    else:
        os.mkdir(sub_dir)


print ("Project directories set up.")

#%%
# Initialize the entire system
# declare a time vector (time window)
t_span = [0,500]
t = np.arange(t_span[0], t_span [1],0.01)

#%%
# Number of DOM/Carbon species:
empty_dic = {}
for i in list(range(4)):
    np.random.seed(i)
    dom_n = np.random.randint(5,20,1)[0]
    bio_n = np.random.randint(4,20,1)[0]
    print(dom_n, bio_n)
    enzparams, vparams, kparams, yparams, mparams = generate_random_parameters(dom_n, bio_n,5)
    dom_initial, biomass_initial = generate_random_initial_conditions (dom_n, bio_n)
    x0 = np.append(dom_initial, biomass_initial)

    carbon_input = generate_random_boundary_conditions(dom_n, 50, method_name = "user_defined")

    seed_dic = {i : {'dom_number': dom_n, 'biomass_number': bio_n,
    'enzyme_production_parameters': enzparams,
    'max_rate_parameters': vparams,
    'sat_const_parameters': kparams,
    'efficiency_parameters': yparams,
    'mortality_parameters': mparams,
    'initial_conditions_dom': dom_initial,
    'initial_conditions_biomass': biomass_initial,
    'carbon_input_boundary': carbon_input}}
    
    empty_dic.update(seed_dic)

    trial = rn(maximum_capacity=5,
            #carbon_mol_bio = 10,
            carbon_num = dom_n,
            bio_num = bio_n,
            carbon_input = carbon_input,
            sigmoid_coeff_stolpovsky = 0.01,
            enzyme_production_rate_constant = enzparams[0],
            efficiency_bio_uptake = enzparams[1],
            necromass_distribution="notequal")
    trial.set_rate_constants(vparams, kparams, yparams,mparams)
    trial.identify_components_natures()
    solution = trial.solve_network(x0, t_span, t)

    tim = solution.t
    sim_array = solution.y.T

    # Shannon diversity
    S, DOC, TOC = diversity_carbon(sim_array, dom_n, bio_n)
    
    # Figures
    # Time series of biomass and DOC concentration 
    C1 = sim_array[:,trial.most_labile_c,]
    C2 = sim_array[:,trial.least_labile_c]
    C = sim_array[:,:dom_n]
    B = sim_array[:,dom_n:]

    if sim_array.any()<0:
        print("Negative values")
    else:
        proportion = B/B.sum(axis=1, keepdims = True)
        plog = np.log(proportion)
        puct = proportion*np.log(proportion)
        #Shannon = -np.sum(proportion*np.log(proportion), axis = 1)
        print(np.shape(puct))#, Shannon)
    #total_C_stock = np.sum(C,axis=1) + np.sum(B, axis=1)
    #C_stock = np.sum(C,axis=1)

    plt.figure()
    plt.plot (tim,B, linestyle = '--')#, label = "Bacteria")
    plt.plot (tim,C, linestyle = '-')#, label = "Bacteria")
    plt.xlabel ("Time (days)")
    plt.ylabel ("Concentration")
    plt.legend()
    plt.text(s = "Seed: "+str(i), x = 80, y = 800, fontsize = 12)
    solid_line = mlines.Line2D([], [], linestyle = '-', color='grey', linewidth=1.5, label='DOM')
    dashed_line = mlines.Line2D([], [], linestyle = '--', color='grey', linewidth=1.5, label='Biomass')
    first_legend = plt.legend(handles=[solid_line, dashed_line])
    plt.savefig(os.path.join(figures_dir, str(i) + "_Time_series_concentration.png"), dpi = 300)

    # Time series of biomass and DOC concentration 
    fig, axes = plt.subplots(nrows=3, ncols=1, sharex = True, figsize = (6,8))
    axes.flat[0].plot(tim, S)
    axes.flat[1].plot(tim, DOC)
    axes.flat[2].plot(tim, TOC)
    axes.flat[0].set_ylabel("Shannon diversity")
    axes.flat[1].set_ylabel("DOC")
    axes.flat[2].set_ylabel("TOC")
    axes.flat[2].set_xlabel("Time")
    axes.flat[0].text(s = "Seed: "+str(i), x = 0.8, y = 0.8, transform = axes.flat[0].transAxes,
    fontsize = 12)
    plt.savefig(os.path.join(figures_dir, str(i) + "_diversity_concentration.png"), dpi = 300)

#%%
index_to_plot = -8000
fig, axes = plt.subplots(nrows=2, ncols=1, sharex = True)
axes.flat[0].scatter(x = S[index_to_plot:], y = DOC[index_to_plot:])
axes.flat[1].scatter(x = S[index_to_plot:], y = TOC[index_to_plot:])
axes.flat[1].set_xlabel("Shannon diversity")
axes.flat[0].set_ylabel("DOC")
axes.flat[1].set_ylabel("TOC")
#%%
filename = os.path.join(output_dir, "diversity_data.pkl")
diversity_data = pd.read_pickle(filename)

#%%
plt.scatter(x = 'Shannon',y='DOC', data = diversity_data, s = 'carbon_species')