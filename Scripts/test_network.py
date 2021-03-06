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
details_subfolder = 'input_5000_test_0_1'
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
t_span = [0,5000]
t = np.arange(t_span[0], t_span [1],10)

#%%
empty_dic = {}
for i in list(range(107, 108)):
    np.random.seed(i)
    dom_n = 18#np.random.randint(5,20,1)[0]
    bio_n = 30#int(dom_n*1.5)#np.random.randint(4,20,1)[0]
    print(dom_n, bio_n)
    S_witch = np.random.randint(2, size=(dom_n, bio_n))
    #%%
    dom_initial, biomass_initial = generate_random_initial_conditions (dom_n, bio_n, 1000, 100, 5000, 10)

    print(np.sum(dom_initial), np.sum(biomass_initial))
    ox_state, enzparams, zparams, vparams, kparams, mparams = generate_random_parameters(dom_n, bio_n,50*np.sum(biomass_initial))
    enzparams = 0.4
    x0 = np.append(dom_initial, biomass_initial)

    carbon_input = generate_random_boundary_conditions(dom_n, 1*np.sum(dom_initial)/(100*dom_n), method_name = "user_defined")

    S_witch = 1#np.random.randint(2, size=(dom_n, bio_n))
    #np.random.shuffle(S_witch)
    #S_witch = np.ones((dom_n, bio_n))
    
    trial = rn(maximum_capacity=5,
            #carbon_mol_bio = 10,
            carbon_num = dom_n,
            bio_num = bio_n,
            carbon_input = carbon_input,
            sigmoid_coeff_stolpovsky = 0.01,
            necromass_distribution="oxidation_state")
    trial.set_rate_constants(ox_state, enzparams, zparams, vparams, kparams, mparams)
    trial.rearrange_constants()
    trial.identify_components_natures(recalcitrance_criterion="oxidation_state")
    trial.microbe_carbon_switch(S_witch)
    solution = trial.solve_network(x0, t_span, t)

    seed_dic = {i : {'dom_number': dom_n, 'biomass_number': bio_n,
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
    
    tim = solution.t
    sim_array = solution.y.T

    # Shannon diversity
    #S, DOC, TOC = diversity_carbon(sim_array, dom_n, bio_n)
    
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

    fig = plt.figure()
    plt.plot (tim,B, linestyle = '--')#, label = "Bacteria")
    plt.plot (tim,C, linestyle = '-')#, label = "Bacteria")
    plt.xlabel ("Time (days)")
    plt.ylabel ("Concentration")
    #plt.legend()
    plt.text(s = "Seed: "+str(i), x = 0.8, y = 0.8,
    fontsize = 12)
    solid_line = mlines.Line2D([], [], linestyle = '-', color='grey', linewidth=1.5, label='DOM')
    dashed_line = mlines.Line2D([], [], linestyle = '--', color='grey', linewidth=1.5, label='Biomass')
    first_legend = plt.legend(handles=[solid_line, dashed_line])
    plt.savefig(os.path.join(figures_dir, str(i) + "_Time_series_concentration.png"), dpi = 300)

#%%
dom_n = 18#np.random.randint(5,20,1)[0]
bio_n = dom_n*1.5#np.random.randint(4,20,1)[0]
print(dom_n, bio_n)
S_witch = np.random.randint(2, size=(dom_n, bio_n))
#%%
dom_initial, biomass_initial = sim_array[-1, :dom_n],sim_array[-1, dom_n:]
x0 = np.append(dom_initial, biomass_initial)
#%%
#trial = rn(maximum_capacity=5,
#        #carbon_mol_bio = 10,
#        carbon_num = dom_n,
#        bio_num = bio_n,
#        carbon_input = carbon_input,
#        sigmoid_coeff_stolpovsky = 0.01,
#        necromass_distribution="oxidation_state")
ox_state = trial.recalcitrance_state
enzparams = trial.v_enz
zparams = trial.z
vparams = trial.v_params
kparams = trial.k_params
mparams = trial.m_params
yparams = trial.y_params
#%%
trial.set_rate_constants(yparams, ox_state, enzparams, zparams, vparams, kparams, mparams)
trial.constants_previously_defined()
trial.identify_components_natures(recalcitrance_criterion="oxidation_state")
trial.microbe_carbon_switch(S_witch)
solution_10k = trial.solve_network(x0, t_span, t)

tim = solution_10k.t
sim_array = solution_10k.y.T

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

fig = plt.figure()
plt.plot (tim,B, linestyle = '--')#, label = "Bacteria")
plt.plot (tim,C, linestyle = '-')#, label = "Bacteria")
plt.xlabel ("Time (days)")
plt.ylabel ("Concentration")
#plt.legend()
plt.text(s = "Seed: "+str(i), x = 0.8, y = 0.8,
fontsize = 12)
solid_line = mlines.Line2D([], [], linestyle = '-', color='grey', linewidth=1.5, label='DOM')
dashed_line = mlines.Line2D([], [], linestyle = '--', color='grey', linewidth=1.5, label='Biomass')
first_legend = plt.legend(handles=[solid_line, dashed_line])
plt.savefig(os.path.join(figures_dir, str(i) + "_Time_series_concentration_next_ts.png"), dpi = 300)


#%%
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
# %%
c_n = 3
rng = np.random.default_rng()
for c_b_ratio in [3]:
    N = int(c_n * c_b_ratio)
    for baseline, activity_pc in zip(["b_2", "b_3", "b_4", "b_5"] ,[0.1, 0.3, 0.5, 0.7]):
        K = math.ceil(activity_pc*N)
        s = np.array([1]*K + [0] * (N-K))#.reshape(c_n, int(c_n*c_b_ratio))
        switch = np.zeros((c_n, N))#[s]*c_n
        #rng.permutation(switch, axis = 1)
        i=0
        while i<c_n:
           switch[i,:] = rng.permutation(s)
           i+=1
        for j in list(range(3)):
            x = rng.permutation(switch, axis = 1)
            print (c_b_ratio, N, baseline, K)
            #print(np.shape(switch))
            print(x)

        #for baseline, activity_pc in zip(["b_1", "b_2", "b_3", "b_4", "b_5"] ,[ 1, 0.1, 0.3, 0.5, 0.7]):
    