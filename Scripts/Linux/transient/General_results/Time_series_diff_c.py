## Import libraries
import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")

ip = 0
seed_sim_list = [610229235, 983307757, 643338060, 714504443, 277077803, 898393994, 420,13012022,13061989]
styles = {"a":"darkgoldenrod", "b":"purple", "c":"indianred", "d":"steelblue", "e":"orange"}
cn_list = [3,6,12,18]
bio_n_series = [4,8,16,32]
init_dom_list = [1000,2000,5000,10000,15000]
grey_carbon = mlines.Line2D([], [], linestyle = '-', color = "grey", marker = None, label='Carbon')
grey_biomass = mlines.Line2D([], [], linestyle = '--', color = "grey", marker = None, label='Biomass')

time_span = np.linspace(0,36505, int(36500/5))
xticks_plot = time_span[::730].astype(int)
xticks_label = np.linspace(0,100,100)[::10].astype(int)

for seed_sim in seed_sim_list:
    for c_n in cn_list:
        details_subfolder = 'v2_1c_adaptation_carbon_' + str(c_n) + '_'+str(seed_sim) + '_ip_0'
        simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
        figures_dir = os.path.join(project_dir, "figures", details_subfolder)
        hr = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'r')
        for base, act_label in zip(["b_2", "b_3", "b_4"], ["10%","25%", "50%", "75%"]):
            #print(base)
            for dom_init in init_dom_list:
                #print(case)
                for c_b_r in bio_n_series:
                    cases = ["a", "b", "c", "d", "e"]
                    for case in cases:
                        fig, ax1 = plt.subplots(figsize = (6,4))
                        plt.title (act_label + " microbial activity")
                        ax1.set_xlabel ("Time [T]")
                        ax1.set_ylabel ("Carbon [N $L^{-3}$]")
                        ax2 = ax1.twinx()
                        ax2.set_ylabel ("Biomass [N $L^{-3}$]")
                        sim_id = base + "_" + case + "_/bio_n_" + str(c_b_r) + "/dom_initial_" + str(dom_init) + "/seed_" + str(seed_sim)
                        if hr[sim_id]:
                            sim_data = hr[sim_id]
                            init = sim_data['initial_conditions']
                            paras = sim_data['parameters']
                            dom_n, bio_n = sim_data['species_number']
                            x = sim_data['solution']
                            C = np.asarray(x['dom'])
                            B = np.asarray(x['biomass'])                
                            ax1.plot(time_span, C, '-')
                            ax1.set_yscale("log")
                            ax1.set_xticks(xticks_plot)
                            ax1.set_xticklabels(xticks_label)
                            ax2.plot(time_span, B,'--')
                            ax2.set_xticks(xticks_plot)
                            ax2.set_xticklabels(xticks_label)
                        plt.tight_layout()
                        plt.savefig(os.path.join(figures_dir, "Time_series_"+base+"_case_"+case+"_bio_"+str(c_b_r)+"_"+str(dom_init)))
                        plt.close()
        hr.close
 