## Import libraries
import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

## LOAD RESULTS
#project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
project_dir = "C:/Users/swkh9804/Documents/Projects/HoliSoils/data"

styles = {"a":"darkgoldenrod", "b":"purple", "c":"indianred", "d":"steelblue", "e":"orange"}
subfolders = ['carbon_3', 'carbon_6', 'carbon_8','carbon_15','carbon_18']
#bio_n_series_master = [[2,3,5,6],[3,6,9,12],[4,6,8,16],[3,5,9,12,15,20,30],[3,9,12,15,18]] m_const:50, seed: 13061989
#bio_n_series_master = [[3,5,6],[3,6,9,12],[4,6,8,16],[3,5,9,12,15,20,30],[3,9,12,15,18]] #m_const:25, seed: 13061989
#bio_n_series_master = [[2,3,5,6],[3,6,9,12],[4,6,8,16],[3,5,9,12,15,20,30],[3,9,12,15,18]] #m_const:10, seed: 13061989
bio_n_series_master = [[2,3,5,6],[3,6,9,12],[4,6,8,16],[3,5,9,12,15,20],[3,9,12,15,18]] #m_const:10, seed:420

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

for s, bio_n_series in zip(subfolders, bio_n_series_master):
    details_subfolder = s +"_seed_420_ip_0"
    print(details_subfolder)
    simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
    figures_dir = os.path.join(project_dir, "figures", details_subfolder)
    hr = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'r')
    row = []
    #for base in ["b_1"]:
    #    #print(base)
    #    for case in ["all"]:
    #        #print(case)
    #        for c_b_r in bio_n_series:
    #            #print(c_b_r)
    #            for dom_init in [1000,5000,10000,15000,20000]:
    #                #print(dom_init)
    #                sim_id = base + "_" + case + "_/bio_n_" + str(c_b_r) + "/dom_initial_" + str(dom_init) + "/seed_13061989"
    #                if hr[sim_id]:
    #                    sim_data = hr[sim_id]
    #                    init = sim_data['initial_conditions']
    #                    paras = sim_data['parameters']
    #                    dom_n, bio_n = sim_data['species_number']
    #                    x = sim_data['solution']
    #                    C = np.asarray(x['dom'])
    #                    B = np.asarray(x['biomass'])
    #                    fig, ax1 = plt.subplots()
    #                    ax1.set_xlabel ("Time")
    #                    ax1.set_ylabel ("Carbon")
    #                    ax1.plot(C, '-', label = "Carbon")
    #                    ax2 = ax1.twinx()
    #                    ax2.set_ylabel ("Biomass")
    #                    ax2.plot(B, '--', label = "Biomass")
    #                    fig.tight_layout()
    #                    plt.savefig(os.path.join(figures_dir, "Time_series_"+base+case+"_"+str(c_b_r)+"_"+str(dom_init)))
    #                    plt.close()
    #                    C_sum = np.sum(C, axis = 1)
    #                    B_sum = np.sum(B, axis = 1)
    #                    fig, ax1 = plt.subplots()
    #                    ax1.set_xlabel ("Time")
    #                    ax1.set_ylabel ("Carbon")
    #                    ax1.plot(C_sum, '-', label = "Carbon")
    #                    ax2 = ax1.twinx()
    #                    ax2.set_ylabel ("Biomass")
    #                    ax2.plot(B_sum, '--', label = "Biomass")
    #                    fig.tight_layout()
    #                    plt.savefig(os.path.join(figures_dir, "B_C_sum_Time_series_"+base+case+"_"+str(c_b_r)+"_"+str(dom_init)))
    #                    plt.close()
    for base, act_label in zip(["b_2", "b_3", "b_4", "b_5"], ["10%", "30%", "50%", "70%"]):
        #print(base)
        for dom_init in [1000,5000,10000,15000,20000]:
            #print(case)
            for c_b_r in bio_n_series:
                #print(c_b_r)
                fig, ax1 = plt.subplots()
                plt.title (act_label + "with initial DOM" + str(dom_init) )
                base_id = "b_1_all_/bio_n_" + str(c_b_r) + "/dom_initial_" + str(dom_init) + "/seed_420"
                base_data = hr[base_id]
                x = base_data['solution']
                C_sum = np.sum(np.asarray(x['dom']), axis = 1)
                B_sum = np.sum(np.asarray(x['biomass']), axis = 1)
                ax1.set_xlabel ("Time")
                ax1.set_ylabel ("Carbon")
                ax2 = ax1.twinx()
                ax2.set_ylabel ("Biomass")
                ax1.plot(time_span, C_sum, '-', color = "grey")
                ax2.plot(time_span, B_sum, '--', color = "grey")
                for case in ["a", "b", "c", "d", "e"]:
                    #print(dom_init)
                    sim_id = base + "_" + case + "_/bio_n_" + str(c_b_r) + "/dom_initial_" + str(dom_init) + "/seed_420"
                    if hr[sim_id]:
                        sim_data = hr[sim_id]
                        init = sim_data['initial_conditions']
                        paras = sim_data['parameters']
                        dom_n, bio_n = sim_data['species_number']
                        x = sim_data['solution']
                        C = np.asarray(x['dom'])
                        B = np.asarray(x['biomass'])                  
                        C_sum = np.sum(C, axis = 1)
                        B_sum = np.sum(B, axis = 1)      
                        ax1.plot(time_span, C_sum, '-', color = styles[case])
                        ax2.plot(time_span, B_sum, '--', color = styles[case])
                        #fig.tight_layout()
                plt.legend(handles = linelist, ncol=2)
                plt.savefig(os.path.join(figures_dir, "B_C_sum_Time_series_"+base+"all_cases_"+str(c_b_r)+"_"+str(dom_init)))
                plt.close()
    #for base in ["b_2", "b_3", "b_4", "b_5"]:
    #    #print(base)
    #    for case in ["a", "b", "c", "d", "e"]:
    #        #print(case)
    #        for c_b_r in bio_n_series:
    #            #print(c_b_r)
    #            for dom_init in [1000,5000,10000,15000,20000]:
    #                #print(dom_init)
    #                sim_id = base + "_" + case + "_/bio_n_" + str(c_b_r) + "/dom_initial_" + str(dom_init) + "/seed_13061989"
    #                if hr[sim_id]:
    #                    sim_data = hr[sim_id]
    #                    init = sim_data['initial_conditions']
    #                    paras = sim_data['parameters']
    #                    dom_n, bio_n = sim_data['species_number']
    #                    x = sim_data['solution']
    #                    C = np.asarray(x['dom'])
    #                    B = np.asarray(x['biomass'])
    #                    fig, ax1 = plt.subplots()
    #                    ax1.set_xlabel ("Time")
    #                    ax1.set_ylabel ("Carbon")
    #                    ax1.plot(C, '-', label = "Carbon")
    #                    ax2 = ax1.twinx()
    #                    ax2.set_ylabel ("Biomass")
    #                    ax2.plot(B, '--', label = "Biomass")
    #                    fig.tight_layout()
    #                    plt.savefig(os.path.join(figures_dir, "Time_series_"+base+case+"_"+str(c_b_r)+"_"+str(dom_init)))
    #                    plt.close()
    hr.close
