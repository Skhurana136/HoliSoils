## Import libraries
import os
import numpy as np
import h5py
import matplotlib.pyplot as plt

## LOAD RESULTS
#project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
project_dir = "C:/Users/swkh9804/Documents/Projects/HoliSoils/data"

subfolders = ['carbon_3', 'carbon_6', 'carbon_8','carbon_15','carbon_18']
bio_n_series_master = [[2,3,5,6],[3,6,9,12],[4,6,8,16],[3,5,9,12,15,20,30],[3,9,12,15,18]]

for details_subfolder, bio_n_series in zip(subfolders, bio_n_series_master):
    print(details_subfolder)
    simulations_dir = os.path.join(project_dir, "simulations", details_subfolder+"_ip_1pc")
    figures_dir = os.path.join(project_dir, "figures", details_subfolder+"_ip_1pc")
    hr = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'r')
    row = []
    for base in ["b_1"]:
        #print(base)
        for case in ["all"]:
            #print(case)
            for c_b_r in bio_n_series:
                #print(c_b_r)
                for dom_init in [1000,5000,10000,15000,20000]:
                    #print(dom_init)
                    sim_id = base + "_" + case + "_/bio_n_" + str(c_b_r) + "/dom_initial_" + str(dom_init) + "/seed_13061989"
                    if hr[sim_id]:
                        sim_data = hr[sim_id]
                        init = sim_data['initial_conditions']
                        paras = sim_data['parameters']
                        dom_n, bio_n = sim_data['species_number']
                        x = sim_data['solution']
                        C = np.asarray(x['dom'])
                        B = np.asarray(x['biomass'])
                        fig, ax1 = plt.subplots()
                        ax1.set_xlabel ("Time")
                        ax1.set_ylabel ("Carbon")
                        ax1.plot(C, '-', label = "Carbon")
                        ax2 = ax1.twinx()
                        ax2.set_ylabel ("Biomass")
                        ax2.plot(B, '--', label = "Biomass")
                        fig.tight_layout()
                        plt.savefig(os.path.join(figures_dir, "Time_series_"+base+case+"_"+str(c_b_r)+"_"+str(dom_init)))
                        plt.close()
                        C_sum = np.sum(C, axis = 1)
                        B_sum = np.sum(B, axis = 1)
                        fig, ax1 = plt.subplots()
                        ax1.set_xlabel ("Time")
                        ax1.set_ylabel ("Carbon")
                        ax1.plot(C_sum, '-', label = "Carbon")
                        ax2 = ax1.twinx()
                        ax2.set_ylabel ("Biomass")
                        ax2.plot(B_sum, '--', label = "Biomass")
                        fig.tight_layout()
                        plt.savefig(os.path.join(figures_dir, "B_C_sum_Time_series_"+base+case+"_"+str(c_b_r)+"_"+str(dom_init)))
                        plt.close()
    for base in ["b_2", "b_3", "b_4", "b_5"]:
        #print(base)
        for case in ["a", "b", "c", "d", "e"]:
            #print(case)
            for c_b_r in bio_n_series:
                #print(c_b_r)
                for dom_init in [1000,5000,10000,15000,20000]:
                    #print(dom_init)
                    sim_id = base + "_" + case + "_/bio_n_" + str(c_b_r) + "/dom_initial_" + str(dom_init) + "/seed_13061989"
                    if hr[sim_id]:
                        sim_data = hr[sim_id]
                        init = sim_data['initial_conditions']
                        paras = sim_data['parameters']
                        dom_n, bio_n = sim_data['species_number']
                        x = sim_data['solution']
                        C = np.asarray(x['dom'])
                        B = np.asarray(x['biomass'])
                        fig, ax1 = plt.subplots()
                        ax1.set_xlabel ("Time")
                        ax1.set_ylabel ("Carbon")
                        ax1.plot(C, '-', label = "Carbon")
                        ax2 = ax1.twinx()
                        ax2.set_ylabel ("Biomass")
                        ax2.plot(B, '--', label = "Biomass")
                        fig.tight_layout()
                        plt.savefig(os.path.join(figures_dir, "Time_series_"+base+case+"_"+str(c_b_r)+"_"+str(dom_init)))
                        plt.close()
                        C_sum = np.sum(C, axis = 1)
                        B_sum = np.sum(B, axis = 1)
                        fig, ax1 = plt.subplots()
                        ax1.set_xlabel ("Time")
                        ax1.set_ylabel ("Carbon")
                        ax1.plot(C_sum, '-', label = "Carbon")
                        ax2 = ax1.twinx()
                        ax2.set_ylabel ("Biomass")
                        ax2.plot(B_sum, '--', label = "Biomass")
                        fig.tight_layout()
                        plt.savefig(os.path.join(figures_dir, "B_C_sum_Time_series_"+base+case+"_"+str(c_b_r)+"_"+str(dom_init)))
                        plt.close()
    hr.close
