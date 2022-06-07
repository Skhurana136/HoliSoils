#%%
# ## Import libraries
import numpy as np
import os
import math

import matplotlib.pyplot as plt

from DS.solvers.diff_eqn_system import ReactionNetwork as rn
from DS.solvers.diff_eqn_system import generate_random_parameters
from DS.solvers.diff_eqn_system import generate_random_initial_conditions
from DS.solvers.diff_eqn_system import generate_random_boundary_conditions

#%%
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")

seed_sim_list = [643338060]#[610229235, 983307757, 643338060, 714504443, 277077803, 898393994, 420,13012022,13061989]

cn_list = [18]#[3,6,12,18]
bio_n_series = [4,8,16,32]
ip = 0 #0 random scenarios
init_dom_list = [1000,2000,5000,10000,15000]

c_n = 18
random_seed_number = 643338060
N = 16

#%%
baseline = "b_1"
label = "all"
S_witch_b_1 = 1
dom_initial = 1000

empty_dic = {}

experiment = baseline + "_" + label + "_"

np.random.seed(random_seed_number)
sim = experiment + "/bio_n_"+ str(N) + "/dom_initial_" + str(dom_initial) + "/seed_" + str(random_seed_number)

# declare a time vector (time window)
total_dom_initial = dom_initial
dom_bio_ratio_initial = 10
mean_dom_initial = 1000
mean_bio_initial = 100

dom_n = c_n
bio_n = N
# Initialize the same number of parameters and initial conditions:
dom_initial, biomass_initial = generate_random_initial_conditions (dom_n, bio_n, mean_dom_initial, mean_bio_initial, total_dom_initial, dom_bio_ratio_initial)

ox_state, enzparams, zparams, vparams, kparams, mparams = generate_random_parameters(dom_n, bio_n,5*np.sum(biomass_initial))

x0 = np.append(dom_initial, biomass_initial)

rel_tol_arr = 10**(math.floor(math.log10(np.min(x0)))-6)
abs_tol_arr = 10**-6

carbon_input = generate_random_boundary_conditions(dom_n, 0, method_name = "user_defined")

trial = rn(maximum_capacity=5,carbon_num = dom_n,bio_num = bio_n, carbon_input = carbon_input, necromass_distribution="notequal")

trial.set_rate_constants(ox_state, enzparams, zparams, vparams, kparams, mparams)
trial.rearrange_constants()
trial.identify_components_natures(recalcitrance_criterion="oxidation_state")
trial.reorder_constants_with_comp_nature()
trial.microbe_carbon_switch(S_witch_b_1)


#%%
t_end = 50000
t_step = 1
tim_sol = np.arange(0, t_end,5)
C_sol = np.empty((len(tim_sol),c_n))
B_sol = np.empty((len(tim_sol),N))
C_B_arr = x0
t = 0
while t <= t_end:
    sol_ivp = trial.solve_network(C_B_arr, [0,t_step], np.linspace(0, t_step, 6), solv_method = 'RK45', rel_tol = rel_tol_arr, abs_tol = abs_tol_arr, first_tim_step = 0.001, max_tim_step = 0.5)
    sol_array = sol_ivp.y
    C_B_arr = sol_array[:,-1]
    if t in tim_sol:
        print(random_seed_number, c_n, N, t)
        C_sol[np.argwhere(tim_sol==t),:] = C_B_arr[:c_n]
        B_sol[np.argwhere(tim_sol==t),:] = C_B_arr[c_n:]
    t += t_step
#%%
print(np.shape(C_sol), np.shape(B_sol))
solution = np.concatenate((C_sol, B_sol), axis = 1)
print(np.shape(solution))
#%%
plt.plot(C_sol, '-')
plt.plot(B_sol, '--')
plt.ylim((-1,100))
#%%
plt.plot(C_sol[:-10,5])
plt.plot(C_sol[:-10,7])
#%%
seed_dic = {sim : {'dom_number': dom_n, 'biomass_number': bio_n,
'oxidation_state': ox_state,
'enzyme_production_parameters': trial.v_enz,
'uptake_parameters': trial.z,
'max_rate_parameters': trial.v_params,
'sat_const_parameters': trial.k_params,
'efficiency_parameters': trial.y_params,
'mortality_parameters': trial.m_params,
'initial_conditions_dom': dom_initial,
'initial_conditions_biomass': biomass_initial,
'carbon_input_boundary': carbon_input,
'sim_status':sol_ivp.status,
'message': sol_ivp.message}}