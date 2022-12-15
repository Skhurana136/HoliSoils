## Import libraries
import os
import pandas as pd
import numpy as np
import sys
import h5py
import itertools

def calc_chars(data, timidx):
    x = data['solution']
    C = np.asarray(x['dom'])
    B = np.asarray(x['biomass'])
    NOSC = np.asarray(x['carbon_os'])
    CUE = np.asarray(x['cue'])
    FD = np.asarray(x['func_div'])
    S = np.asarray(x['shannon'])
    Biomass = np.sum(B, axis = 1)
    DOC = np.sum(C, axis = 1)
    TOC = DOC+Biomass
    search_idx=timidx
    nanargwhere = np.argwhere(search_idx==0.)
    search_idx[nanargwhere]=-1
    results_array = np.zeros((timidx.size,8))
    results_array = np.zeros((search_idx.size,12))
    results_array[:,0] = np.asarray([100,90,80,70,60,50,40,30,20,10])
    results_array[:,1] = S[search_idx]
    results_array[:,2] = DOC[search_idx]
    results_array[:,3] = Biomass[search_idx]
    results_array[:,4] = CUE[search_idx]
    results_array[:,5] = FD[search_idx]
    results_array[:,6] = NOSC[search_idx]
    results_array[:,7] = S[0]
    results_array[:,8] = Biomass[0]
    results_array[:,9] = CUE[0]
    results_array[:,10] = FD[0]
    results_array[:,11] = NOSC[0]
    #results_array[nanargwhere,1:] = np.nan
    return results_array

def create_pd_dataset(data_val, c_val, b_val, seed_val, sim_val, doc_i_val):
    dataset = pd.DataFrame(data = data_val, columns = ["%C", "S", "DOC", "Biomass", "CUE", "FD", "NOSC","S_initial", "Biomass_initial", "CUE_initial", "FD_initial", "NOSC_initial"])
    #dataset = pd.melt(data_table, id_vars=["%C"], value_vars=["S", "DOC", "Biomass", "CUE", "FD","NOSC"], var_name = "Characteristic", value_name = "Char_value")
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
    return data_cbssd

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient", sys.argv[1])
#project_dir = os.path.join('/proj', 'hs_micro_div_072022', 'Project_data', 'transient', sys.argv[1])
results_dir = os.path.join(project_dir, "results")
filestring = 'competition_adaptation_carbon_' #null
seed_sim_list = [610229235, 983307757, 643338060, 714504443, 277077803, 898393994, 420,13012022,13061989]

cn_list = [3,6,12,18]
bio_n_series = [4,8,12,16,20,24,28,32]
ip = 0
init_dom_list = [1000,2000,5000,10000,15000]
transient_switch = 0
input_factor = transient_switch*5/365
tim_file = os.path.join(results_dir, filestring + '_decay_const_c_pools_data_initial_conditions.pkl')
#print(tim_file)
tim_data = pd.read_pickle(tim_file)

files=[]
for c_n, b_n, seed_sim, t_dom_initial in itertools.product(cn_list, bio_n_series,seed_sim_list, init_dom_list):
    # Load all datasets and save their Shannon and diversity indices in a dataframe
    seed_all = 'seed_'+str(seed_sim)
    details_subfolder = filestring + str(c_n) + '_'+str(seed_sim) + '_ip_' + str(ip)
    simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
    hr = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'r')
    #filename = os.path.join(simulations_dir, "seeds_randoms.pkl")
    #seed_details = pd.read_pickle(filename)
    c_b = "bio_n_"+ str(b_n)
    dom_init = "dom_initial_" + str(t_dom_initial)
    doc_input = (t_dom_initial) * input_factor
    tim_subset = tim_data[(tim_data['C_pool']=='TOC')&(tim_data['carbon_species']==c_n)&(tim_data['Seed']==seed_sim)&(tim_data['biomass_species']==b_n)&(tim_data['DOC_initial'].round(decimals=0).astype(int)==int(t_dom_initial))].reset_index()
    base_case = "b_1_all_"
    char_tim_set = tim_subset[tim_subset.Sim_series==base_case]
    #print(c_n, seed_sim, b_n, t_dom_initial, sim, char_tim_set.shape)
    char_tim_arr = char_tim_set['T_loss'].values.astype(int)
    char_tim = np.cumsum(char_tim_arr)
    print(char_tim)
    sim_data = hr[base_case][c_b][dom_init][seed_all]
    rsdbcf = calc_chars(sim_data,char_tim)
    pd_data = create_pd_dataset(rsdbcf, c_n, b_n, seed_sim, base_case, t_dom_initial)
    files.append(pd_data)
    for baseline, label in itertools.product(["b_2", "b_3", "b_4","b_5"],["a", "b", "c","d","e"]):
        sim = baseline + "_" + label + "_"
        #print(sim)
        char_tim = tim_subset[tim_subset.Sim_series==sim]
        char_tim_arr = char_tim_set['T_loss'].values.astype(int)
        char_tim = np.cumsum(char_tim_arr)
        #print(char_tim)
        #print(c_n, seed_sim, b_n, t_dom_initial, sim, char_tim_set.shape)
        sim_data = hr[sim][c_b][dom_init][seed_all]
        rsdbcf = calc_chars(sim_data,char_tim)
        pd_data = create_pd_dataset(rsdbcf, c_n, b_n, seed_sim, sim, t_dom_initial)
        files.append(pd_data)
                
    hr.close()

cue_s_fd_data = pd.concat(files)
print(cue_s_fd_data.columns)
print("The shape of the dataframe is ", cue_s_fd_data.shape)
print("The dataframe contains the following data types ", cue_s_fd_data.dtypes)

filename = os.path.join(results_dir, filestring+"cue_wremaining_c_data.pkl")
cue_s_fd_data.to_pickle(filename)
print ("CUE data is saved here ", filename)

filename = os.path.join(results_dir, filestring+"cue_wremaining_c_data.csv")
cue_s_fd_data.to_csv(filename)
print ("CUE data is saved here ", filename)