## Import libraries
import os
import pandas as pd
import numpy as np
import sys
import h5py 
import itertools

def derive_t_loss(sim_data, loss_criteria):
    paras = sim_data['parameters']
    os_i = np.asarray(paras["oxidation_state"])
    c_os_less_0 = np.where(os_i<0.)
    c_os_gr_0 = np.where(os_i>0.)
    x = sim_data['solution']
    C = np.asarray(x['dom'])
    B = np.sum(np.asarray(x['biomass']), axis=1)
    DOC = np.sum(C, axis = 1)
    TOC = DOC+B
    C_less_0 = np.sum(np.sum(C[:,c_os_less_0], axis = 1),axis=1)
    C_gr_0 = np.sum(np.sum(C[:,c_os_gr_0], axis = 1),axis=1)
    #conc_series = np.array([DOC, C_less_0, C_gr_0]).T
    tim_points_num = len(loss_criteria)
    results_arr = np.empty((tim_points_num,5))
    for c_idx, c_tim_series in zip(list(range(4)), [TOC,DOC,C_less_0, C_gr_0]):
        #print(np.shape(c_tim_series))
        c_tim_series_i = c_tim_series[0]
        #print(c_tim_series_i)
        for tim in list(range(tim_points_num)):
            results_arr[tim, 0] = (loss_criteria[tim]+0.1)*100
            #t_c = np.argwhere(np.round_(c_tim_series/c_tim_series_i, decimals = 2)<=(loss_criteria[tim]))
            t_c = np.argwhere((c_tim_series/c_tim_series_i)<=(loss_criteria[tim]))
            #print(t_c.size)
            if t_c.size > 0:
                results_arr[tim, c_idx+1] = t_c[0][0]
                #print(c_idx, tim, results_arr[tim, c_idx+1])
            else:
                results_arr[tim, c_idx+1] = -1#np.nan
    #print(results_arr)
    subtract = results_arr#np.diff(results_arr, axis = 0)
    subtract[1:,1:] = np.subtract(results_arr[1:,1:],subtract[:-1,1:])
    #print(subtract)
    return subtract

def create_pd_dataset(data_val, c_val, b_val, seed_val, sim_val, doc_i_val):
    data_table = pd.DataFrame(data = data_val, columns = ["%C","TOC","DOC", "reduced_C", "oxidized_C"])
    dataset = pd.melt(data_table, id_vars=["%C"], value_vars=["TOC","DOC", "reduced_C", "oxidized_C"], var_name = "C_pool", value_name = "T_loss")
    carb_ser = pd.Series([c_val]*dataset.shape[0], copy=False,name = "carbon_species")
    bio_ser = pd.Series([b_val]*dataset.shape[0], copy=False,name = "biomass_species")
    seed_ser = pd.Series([seed_val]*dataset.shape[0], copy=False,name = "Seed")
    sim_ser = pd.Series([sim_val]*dataset.shape[0], copy = False, name = 'Sim_series')
    doc_ser = pd.Series([doc_i_val]*dataset.shape[0], copy = False, name = 'DOC_initial')
    data_c = dataset.join(carb_ser)
    data_cb = data_c.join(bio_ser)
    data_cbs = data_cb.join(seed_ser)
    data_cbss = data_cbs.join(sim_ser)
    data_cbssd = data_cbss.join(doc_ser)
    data_cbssd['T_loss'] = data_cbssd['T_loss'].clip(lower=0)
    #print(data_cbssd)
    return data_cbssd

## LOAD RESULTS
project_dir = os.path.join("D:\Projects", "HoliSoils","data","transient", sys.argv[1])
#project_dir = os.path.join('/proj', 'hs_micro_div_072022', 'Project_data', 'transient', sys.argv[1])
results_dir = os.path.join(project_dir, "results")
filestring =  'competition_adaptation_carbon_' #null
seed_sim_list = [610229235, 983307757, 643338060, 714504443, 277077803, 898393994, 420,13012022,13061989]

cn_list = [3,6,12,18]
bio_n_series = [4,8,12,16,20,24,28,32]
ip = 0
init_dom_list = [1000,2000,5000,10000,15000]
transient_switch = 0
input_factor = transient_switch*5/365
loss_criteria_series = [0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1, 0.]
result_fstring = "_loss_temporal"

files=[]
dictionary_iter = 0
for c_n, b_n, seed_sim, t_dom_initial in itertools.product(cn_list, bio_n_series,seed_sim_list, init_dom_list):
    # Load all datasets and save their Shannon and diversity indices in a dataframe
    seed_all = 'seed_'+str(seed_sim)        
    details_subfolder = filestring + str(c_n) + '_'+str(seed_sim) + '_ip_' + str(ip)
    simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
    hr = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'r')
    filename = os.path.join(simulations_dir, "seeds_randoms.pkl")
    seed_details = pd.read_pickle(filename)
    base_case = "b_1_all_"
    c_b = "bio_n_"+ str(b_n)
    dom_init = "dom_initial_" + str(t_dom_initial)
    doc_input = (t_dom_initial) * input_factor
    #print(c_n, seed_sim, b_n, t_dom_initial)
    if hr[base_case][c_b][dom_init][seed_all]:
        sim_data = hr[base_case][c_b][dom_init][seed_all]
        c_loss_tim_points=derive_t_loss(sim_data, loss_criteria_series)
        print(c_n, b_n, seed_sim, base_case, t_dom_initial)
        #print(c_loss_tim_points)
        pd_data = create_pd_dataset(c_loss_tim_points, c_n, b_n, seed_sim, base_case, t_dom_initial)
        files.append(pd_data)
        dictionary_iter+=1
    #for baseline in ["b_2", "b_3", "b_4","b_5"]:
    #    for label in ["a", "b", "c","d","e"]:
    #        sim = baseline + "_" + label + "_"
    #        if hr[sim][c_b][dom_init][seed_all]:
    #            sim_data = hr[sim][c_b][dom_init][seed_all]
    #            c_loss_tim_points=derive_t_loss(sim_data,loss_criteria_series)
    #            pd_data = create_pd_dataset(c_loss_tim_points, c_n, b_n, seed_sim, sim, t_dom_initial)
    #            files.append(pd_data)
    #            dictionary_iter+=1
    hr.close()
print("Results number processed: ", dictionary_iter)
decay_const_c_pools_data = pd.concat(files)
print("The shape of the dataframe is ", decay_const_c_pools_data.shape)
print("The dataframe contains the following data types ", decay_const_c_pools_data.dtypes)

filename = os.path.join(results_dir, filestring+"_decay_const_c_pools_data_initial_conditions.pkl")
decay_const_c_pools_data.to_pickle(filename)
print ("Temporal decay constant for diff carbon pools are saved here ", filename)

filename = os.path.join(results_dir, filestring+"_decay_const_c_pools_data_initial_conditions.csv")
decay_const_c_pools_data.to_csv(filename, index = False)
print ("Temporal decay constant for diff carbon pools are saved here ", filename)