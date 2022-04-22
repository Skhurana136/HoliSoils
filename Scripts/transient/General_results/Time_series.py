## Import libraries
import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

## LOAD RESULTS
#project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
project_dir = "D:\Projects\HoliSoils\data"

ip = 2
seed_sim_list = [420,13012022,13061989]
styles = {"a":"darkgoldenrod", "b":"purple", "c":"indianred", "d":"steelblue", "e":"orange"}
carbon_num_list = [3, 6, 8, 15, 18]
#bio_n_series_master = [[2,3,5,6],[3,6,9,12],[4,6,8,16],[3,5,9,12,15,20,30],[3,9,12,15,18]] m_const:50, seed: 13061989
#bio_n_series_master = [[3,5,6],[3,6,9,12],[4,6,8,16],[3,5,9,12,15,20,30],[3,9,12,15,18]] #m_const:25, seed: 13061989
#bio_n_series_master = [[2,3,5,6],[3,6,9,12],[4,6,8,16],[3,5,9,12,15,20,30],[3,9,12,15,18]] #m_const:10, seed: 13061989
bio_n_series_master = [[2,3,5,6],[3,6,9,12],[4,6,8,12,16],[3,5,9,12,15,20,25,30],[3,6,9,12,15,18,24,30]] #m_const:10, seed:420
init_dom_list = [1000,2000,5000,10000,15000]
grey_carbon = mlines.Line2D([], [], linestyle = '-', color = "black", marker = None, label='Carbon')
grey_biomass = mlines.Line2D([], [], linestyle = '--', color = "black", marker = None, label='Biomass')
grey_patch = mpatches.Patch(color="grey", label= '100%')#, alpha = 0.5)
gold_patch = mpatches.Patch(color="darkgoldenrod", label= 'a')#, alpha = 0.5)
purple_patch = mpatches.Patch(color="purple", label= 'b')#, alpha = 0.5)
red_patch = mpatches.Patch(color="indianred", label= 'c')#, alpha = 0.5)
blue_patch = mpatches.Patch(color="steelblue", label= 'd')#, alpha = 0.5)
orange_patch = mpatches.Patch(color = "orange", label =  'e')#, alpha = 0.5)
linelist = [grey_carbon, grey_biomass, grey_patch, gold_patch, purple_patch, red_patch, blue_patch,orange_patch]

time_span = np.linspace(0,50000, 10000)

for seed_sim in seed_sim_list:
    for s, bio_n_series in zip(carbon_num_list[:2], bio_n_series_master[:2]):
        details_subfolder = "carbon_"+ str(s) +"_" + str(seed_sim) + "_ip_" + str(ip)
        print(details_subfolder)
        simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
        figures_dir = os.path.join(project_dir, "figures", details_subfolder)
        hr = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'r')
        row = []
        for base, act_label in zip(["b_2", "b_3", "b_4", "b_5"], ["10%", "30%", "50%", "70%"]):
            #print(base)
            for dom_init in init_dom_list:
                #print(case)
                for c_b_r in bio_n_series:
                    np_list = []
                    fig, ax1 = plt.subplots()
                    plt.title (act_label + "with initial DOM" + str(dom_init) )
                    base_id = "b_1_all_/bio_n_" + str(c_b_r) + "/dom_initial_" + str(dom_init) + "/seed_" + str(seed_sim)
                    base_data = hr[base_id]
                    x = base_data['solution']
                    C_sum_base = np.sum(np.asarray(x['dom']), axis = 1)
                    B_sum_base = np.sum(np.asarray(x['biomass']), axis = 1)
                    ax1.set_xlabel ("Time")
                    ax1.set_ylabel ("Carbon")
                    ax2 = ax1.twinx()
                    ax2.set_ylabel ("Biomass")
                    if np.shape(C_sum_base)[0] == time_span.size:
                        ax1.plot(time_span, C_sum_base, '-', color = "grey")
                        ax2.plot(time_span, B_sum_base, '--', color = "grey")
                    else:
                        np_list.append("base")
                    for case in ["a", "b", "c", "d", "e"]:
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
                            if np.shape(C)[0] == time_span.size:
                                C_sum = np.sum(C, axis = 1)
                                B_sum = np.sum(B, axis = 1)      
                                ax1.plot(time_span, C_sum, '-', color = styles[case])
                                ax2.plot(time_span, B_sum, '--', color = styles[case])
                            else:
                                np_list.append(case)
                    plt.legend(handles = linelist, ncol=2)
                    if len(np_list)>0:
                        print_string = ','.join(item for item in np_list)
                        plt.text(s = "NP:"+ print_string, x= 0, y = -0.5)
                    plt.savefig(os.path.join(figures_dir, "B_C_sum_Time_series_"+base+"all_cases_"+str(c_b_r)+"_"+str(dom_init)))
                    plt.close()
        hr.close
