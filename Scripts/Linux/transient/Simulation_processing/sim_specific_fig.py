#%%
## Import libraries
import os
import numpy as np
import pandas as pd
import h5py

import matplotlib.pyplot as plt
import seaborn as sns

# Functions to load
def lflatten(t):
    return [item for sublist in t for item in sublist]

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
simulations_dir = os.path.join(project_dir, "simulations")
results_dir = os.path.join(project_dir, "results")
figures_dir = os.path.join(project_dir, "figures")

ip = 0
seed_sim_list = [610229235, 983307757, 643338060, 714504443, 277077803, 898393994]#, 420,13012022,13061989]
styles = {"a":"darkgoldenrod", "b":"purple", "c":"indianred", "d":"steelblue", "e":"orange"}
carbon_num_list = [3, 6, 8, 11, 15, 18]
bio_n_series_master = [[2,3,5,6],[3,6,9,12],[4,6,8,12,16],[3,5,9,12,15,20,25,30],[3,6,9,12,15,18,24,30]]
init_dom_list = [1000,2000,5000,10000,15000]
time_span = np.linspace(0,50000, 10000)
#%%
for seed_sim in seed_sim_list:
    #for s, bio_n_series in zip(carbon_num_list, bio_n_series_master):
    s = carbon_num_list[-1]
    bio_n_series = bio_n_series_master[-1]
    details_subfolder = "carbon_"+ str(s) +"_" + str(seed_sim) + "_ip_" + str(ip)
    print(details_subfolder)
    simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
    figures_dir = os.path.join(project_dir, "figures", details_subfolder)
    hr = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'r')
    #row = []
    #for base, act_label in zip(["b_2", "b_3", "b_4", "b_5"], ["10%", "30%", "50%", "70%"]):
    base = "b_2"
    act_label = "10%"
        #print(base)
    for dom_init in init_dom_list[:1]:
        #print(case)
        for c_b_r in bio_n_series:
            np_list = []
            fig, ax1 = plt.subplots()
            plt.title (str(c_b_r) + " microbes with seed " + str(seed_sim))
            base_id = "b_1_all_/bio_n_" + str(c_b_r) + "/dom_initial_" + str(dom_init) + "/seed_" + str(seed_sim)
            base_data = hr[base_id]
            x = base_data['solution']
            C_sum_base = np.sum(np.asarray(x['dom']), axis = 1)
            B_sum_base = np.sum(np.asarray(x['biomass']), axis = 1)
            ax1.set_xlabel ("Time")
            ax1.set_ylabel ("Carbon")
            ax2 = ax1.twinx()
            ax2.set_ylabel ("Biomass")
            ax1.plot(C_sum_base, '-', color = "grey")
            ax2.plot(B_sum_base, '--', color = "grey")
            plt.savefig(os.path.join(figures_dir, str(c_b_r) + " microbes with seed " + str(seed_sim) + ".png"))
            #plt.close()
    hr.close
#%%
not_work = {}
work = {}
for seed_sim in seed_sim_list:
    #for s, bio_n_series in zip(carbon_num_list, bio_n_series_master):
    s = carbon_num_list[-1]
    bio_n_series = bio_n_series_master[-1]
    details_subfolder = "carbon_"+ str(s) +"_" + str(seed_sim) + "_ip_" + str(ip)
    print(details_subfolder)
    simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
    figures_dir = os.path.join(project_dir, "figures", details_subfolder)
    hr = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'r')
    #row = []
    #for base, act_label in zip(["b_2", "b_3", "b_4", "b_5"], ["10%", "30%", "50%", "70%"]):
    base = "b_2"
    act_label = "10%"
        #print(base)
    for dom_init in init_dom_list[:1]:
        #print(case)
        for c_b_r in bio_n_series:
            np_list = []
            base_id = "b_1_all_/bio_n_" + str(c_b_r) + "/dom_initial_" + str(dom_init) + "/seed_" + str(seed_sim)
            base_data = hr[base_id]
            x = base_data['solution']
            params = base_data['parameters']
            ox_state = params['oxidation_state'][:]
            enzparams = params['exo_enzyme_rate'][:]
            zparams = params['carbon_uptake'][:]
            vparams = params['max_rate'][:]
            kparams = params['half_saturation'][:]
            mparams = params['mortality_const'][:]
            yparams = params['yield_coefficient'][:]
            b_n = base_data['species_number'][1]
            para_dic = {'ox_state':ox_state.flatten(), 'enzparams':enzparams.flatten(), 'zparams': zparams.flatten(), 'vparams':vparams.flatten(), 'kparams':kparams.flatten(), 'mparams': mparams.flatten(), 'yparams': yparams.flatten()}
            if np.shape(x['dom'])[0]<time_span.size:
                not_work[str(seed_sim)+str(b_n)] = para_dic
            else:
                work[str(seed_sim)+str(b_n)]= para_dic
    hr.close

#%%
def plot_distribution (param_id):
    seed_not_work = list(not_work.keys())
    seed_work = list(work.keys())
    for s in seed_not_work:
        sns.kdeplot(not_work[s][param_id], shade = False)
    for s in seed_work:
        sns.kdeplot(work[s][param_id], label= s, shade = False, linestyle = "--")
    plt.title(param_id)
    
#%%
plot_distribution ('enzparams')
#%%
plot_distribution ('zparams')
#%%
plot_distribution ('kparams')
#%%
plot_distribution('mparams')
#%%
plot_distribution ('yparams')
#%%
# Load the carbon removal dataset
filename = os.path.join(results_dir, "combined_dataset.pkl")
diversity_data = pd.read_pickle(filename)
diversity_data['DOC_initial_int'] = round(diversity_data.DOC_initial, -3)
#%%
diversity_data = diversity_data.drop_duplicates()
base = diversity_data[diversity_data.Sim_series=='b_1_all_']
low_doc = base[base.DOC_initial_int==1000]
c18 = low_doc[low_doc.carbon_species==18]
sns.kdeplot(c18['DOC_removal'])
#%%
seed_not_work = list(not_work.keys())
seed_work = list(work.keys())
#%%
plt.figure()
s = seed_not_work[0][:-1]
b_nw = [6,15]
doc_seed = c18[c18.Seed == int(s)]
for b in b_nw:
    ye = doc_seed[doc_seed.biomass_species == b]['DOC_removal']
    xe = not_work[s+str(b)]['kparams']
    #for xe, ye in zip(xe, ye):
    plt.scatter(xe, [ye]*len(xe), label = b)
plt.legend(title = "Biomass_species")
plt.show()
#%%
def list_paras_to_plot(x_plot, y_plot):
    xnw, ynw, xw, yw = [], [], [], []
    for s in seed_not_work:
        ye = not_work[s][x_plot]
        xe = not_work[s][y_plot]
        xnw.append(xe.flatten())
        ynw.append(ye.flatten())
    for s in seed_work:
        ye = work[s][x_plot]
        xe = work[s][y_plot]
        xw.append(xe.flatten())
        yw.append(ye.flatten())
    xnw = lflatten(xnw)
    ynw = lflatten(ynw)
    xw = lflatten(xw)
    yw = lflatten(yw)
    
    return xnw, ynw, xw, yw    

def plot_paras (xnw, ynw, xw, yw, xl, yl):
    fig, axes = plt.subplots(1,2, figsize = (6,3), sharey = True, sharex= True)
    axes[0].hist2d(xnw, ynw, label = s, bins=(25, 25))#scatter(xe, ye, label = s, alpha = 0.7)
    axes[0].set_title("Not_work")
    axes[1].set_title("Work")
    axes[1].hist2d(xw, yw, label = s, bins=(25, 25))
    for a in axes[:]:
        a.set_xlabel(xl)
    axes[0].set_ylabel(yl)
    plt.show()

#%%
xvar = "vparams"
yvar = "kparams"
xpnw, ypnw, xpw, ypw = list_paras_to_plot(xvar, yvar)
plot_paras (xpnw, ypnw, xpw, ypw, xvar, yvar)