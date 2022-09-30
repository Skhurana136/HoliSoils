## Import libraries
import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import sys

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient","gen_spec_lognorm_1_5x")
filestring = sys.argv[1]
ip = 0
seed_sim_list = [610229235, 983307757, 643338060, 714504443, 277077803, 898393994, 420,13012022,13061989]
styles = {"a":"darkgoldenrod", "b":"purple", "c":"indianred", "d":"steelblue", "e":"orange"}
cn_list = [3,6,12,18]
bio_n_series = [4]#,8,16,32]
init_dom_list = [1000,2000,5000,10000,15000]
grey_carbon = mlines.Line2D([], [], linestyle = '-', color = "grey", marker = None, label='Carbon')
grey_biomass = mlines.Line2D([], [], linestyle = '--', color = "grey", marker = None, label='Biomass')
grey_patch = mpatches.Patch(color="black", label= 'Baseline')#, alpha = 0.5)
gold_patch = mpatches.Patch(color="darkgoldenrod", alpha = 0.4, label= 'Carbon with \n variable microbial activity')#, alpha = 0.5)
purple_patch = mpatches.Patch(color="purple", label= 'b')#, alpha = 0.5)
red_patch = mpatches.Patch(color="indianred", label= 'c')#, alpha = 0.5)
blue_patch = mpatches.Patch(color="steelblue", alpha = 0.4, label= 'Biomass with \nvariable microbial activity')#, alpha = 0.5)
orange_patch = mpatches.Patch(color = "orange", label =  'e')#, alpha = 0.5)
linelist = [grey_carbon, grey_biomass, grey_patch, gold_patch, blue_patch]#urple_patch, red_patch, ,orange_patch

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
                    C_sum_base = (1- (np.sum(np.asarray(x['dom']), axis = 1)/dom_init))*100
                    B_sum_base = np.sum(np.asarray(x['biomass']), axis = 1)
                    ax1.set_xlabel ("Time [T]")
                    ax1.set_ylabel ("Carbon consumed (%)")
                    ax2 = ax1.twinx()
                    ax2.set_ylabel ("Biomass [N $L^{-3}$]")
                    ax1.plot(time_span, C_sum_base, '-', color = "black")
                    ax2.plot(time_span, B_sum_base, '--', color = "black")
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

                            C_sum = (1-(np.sum(C, axis = 1)/dom_init))*100
                            B_sum = np.sum(B, axis = 1)
                            random_C[:,cases.index(case)] = C_sum 
                            random_B[:,cases.index(case)] = B_sum
                  
                    C_rand_plot_min = np.min(random_C, axis = 1)
                    C_rand_plot_max = np.max(random_C, axis = 1)
                    B_rand_plot_min = np.min(random_B, axis = 1)
                    B_rand_plot_max = np.max(random_B, axis = 1)
                    ax1.plot(time_span, C_rand_plot_min, '-', color = "darkgoldenrod")
                    ax1.plot(time_span, C_rand_plot_max, '-', color = "darkgoldenrod")
                    ax1.fill_between(time_span, y1 = C_rand_plot_max, y2 = C_rand_plot_min, color = "darkgoldenrod", alpha = 0.4)
                    ax2.plot(time_span, B_rand_plot_min, '--', color = "steelblue")    
                    ax2.plot(time_span, B_rand_plot_max, '--', color = "steelblue")
                    ax2.fill_between(time_span, y1 = B_rand_plot_max, y2 = B_rand_plot_min, color = "steelblue", alpha = 0.4)    
                    ax1.set_ylim(-10,100)
                    ax1.set_xticks(xticks_plot)
                    ax1.set_xticklabels(xticks_label)
                    ax2.set_xticks(xticks_plot)
                    ax2.set_xticklabels(xticks_label)
                    ax1.set_xlim(left = -10)
                    ax2.set_xlim(left = -10)
                    plt.legend(handles = linelist, ncol=2, bbox_to_anchor = (0.9,-0.18))
                    plt.tight_layout()
                    plt.savefig(os.path.join(figures_dir, "B_C_sum_Time_series_"+base+"all_cases_"+str(c_b_r)+"_"+str(dom_init)))
                    plt.close()
        hr.close
 