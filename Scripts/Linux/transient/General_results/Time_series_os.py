## Import libraries
import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import sys
## LOAD RESULTS
project_dir = os.path.join('/proj', 'hs_micro_div_072022', 'Project_data', 'transient', 'func_div')
filestring = sys.argv[1]

ip = 0
seed_sim_list = [610229235, 983307757, 643338060, 714504443, 277077803, 898393994, 420,13012022,13061989]
styles = {"a":"darkgoldenrod", "b":"purple", "c":"indianred", "d":"steelblue", "e":"orange"}
cn_list = [3,6,12,18]
bio_n_series = [4,8,16,32]
init_dom_list = [1000,2000,5000,10000,15000]
grey_patch = mpatches.Patch(color="black", label= 'Baseline')
gold_patch = mpatches.Patch(color="darkgoldenrod", alpha = 0.4, label= 'OS with \n variable microbial activity')#, alpha = 0.5)
linelist = [grey_patch, gold_patch]

time_span = np.linspace(0,36505, int(36500/5))
xticks_plot = time_span[::730].astype(int)
xticks_label = np.linspace(0,100,100)[::10].astype(int)

for seed_sim in seed_sim_list:
    for c_n in cn_list:
        details_subfolder = filestring + '_carbon_' + str(c_n) + '_'+str(seed_sim) + '_ip_0'
        simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
        figures_dir = os.path.join(project_dir, "figures", details_subfolder)
        hr = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'r')

        for base, act_label in zip(["b_2", "b_3", "b_4"], ["10%","25%", "50%", "75%"]):
            #print(base)
            for dom_init in init_dom_list:
                #print(case)
                for c_b_r in bio_n_series:
                    np_list = []
                    fig, ax1 = plt.subplots(figsize = (6,4))
                    plt.title (act_label + " microbial activity")# with initial DOM" + str(dom_init) )
                    base_id = "b_1_all_/bio_n_" + str(c_b_r) + "/dom_initial_" + str(dom_init) + "/seed_" + str(seed_sim)
                    base_data = hr[base_id]
                    x = base_data['solution']
                    paras_base = base_data["parameters"]
                    #C_sum_base = (1- (np.sum(np.asarray(x['dom']), axis = 1)/dom_init))*100
                    #B_sum_base = np.sum(np.asarray(x['biomass']), axis = 1)
                    C = x["dom"]
                    C_sum = np.sum(C, axis = 1)
                    mean_wt_os_base = np.sum(np.array(C)*np.array(paras_base["oxidation_state"]).T, axis = 1)/C_sum
                    ax1.plot(time_span, mean_wt_os_base, '-', color = "black", label = "baseline")
                    ax1.set_xticks(xticks_plot)
                    ax1.set_xticklabels(xticks_label)
                    ax1.set_xlabel ("Time [T]")
                    random_C = np.empty((time_span.size, 5))
                    random_B = np.empty((time_span.size, 5))
                    cases = ["a", "b", "c", "d", "e"]
                    for case in cases:
                        #print(dom_init)
                        sim_id = base + "_" + case + "_/bio_n_" + str(c_b_r) + "/dom_initial_" + str(dom_init) + "/seed_" + str(seed_sim)
                        if hr[sim_id]:
                            sim_data = hr[sim_id]
                            init = sim_data['initial_conditions']
                            paras = sim_data['parameters']
                            dom_n, bio_n = sim_data['species_number']
                            x = sim_data['solution']
                            C = np.asarray(x['dom'])
                            B = np.asarray(x['biomass'])                
                            C_sum = np.sum(C, axis = 1)
                            mean_wt_os = np.sum(np.array(C)*np.array(paras["oxidation_state"]).T, axis = 1)/C_sum
                            random_C[:,cases.index(case)] = mean_wt_os
                  
                    C_rand_plot_min = np.min(random_C, axis = 1)
                    C_rand_plot_max = np.max(random_C, axis = 1)
                    ax1.plot(time_span, C_rand_plot_min, '-', color = "darkgoldenrod")
                    ax1.plot(time_span, C_rand_plot_max, '-', color = "darkgoldenrod")
                    ax1.fill_between(time_span, y1 = C_rand_plot_max, y2 = C_rand_plot_min, color = "darkgoldenrod", alpha = 0.4)
                    ax1.set_xticks(xticks_plot)
                    ax1.set_xticklabels(xticks_label)
                    ax1.set_xlim(left = -10)
                    plt.legend(handles = linelist, ncol=2, bbox_to_anchor = (0.9,-0.18))
                    plt.tight_layout()
                    plt.savefig(os.path.join(figures_dir, "OS_Time_series_"+base+"all_cases_"+str(c_b_r)+"_"+str(dom_init)))
                    plt.close()
        hr.close
 