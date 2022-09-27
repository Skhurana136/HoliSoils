## Import libraries
import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import sys

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient" ,"gen_spec_skew")
filestring = sys.argv[1]

ip = 0
styles = {4:"darkgoldenrod", 8:"orange", 12:"indianred", 20:"purple", 28:"steelblue", 32: "darkblue"}
cmap_bio = plt.cm.get_cmap('Blues')

seed_sim_list = [610229235, 983307757, 643338060, 714504443, 277077803, 898393994, 420,13012022,13061989]
cn_list = [3,6,12,18]
bio_n_series = [4,8,12,16,20,24,28,32]
init_dom_list = [1000,2000,5000,10000,15000]
grey_carbon = mlines.Line2D([], [], linestyle = '-', color = "black", marker = None, label='Carbon')
grey_biomass = mlines.Line2D([], [], linestyle = '--', color = "black", marker = None, label='Biomass')
gold_patch = mpatches.Patch(color=cmap_bio(4*8), label= '4')
purple_patch = mpatches.Patch(color=cmap_bio(8*8), label= '8')
red_patch = mpatches.Patch(color=cmap_bio(12*8), label= '12')
blue_patch = mpatches.Patch(color=cmap_bio(20*8), label= '20')
blue_patch1 = mpatches.Patch(color=cmap_bio(28*8), label= '28')
blue_patch2 = mpatches.Patch(color=cmap_bio(32*8), label= '32')
linelist = [grey_carbon, grey_biomass]
biomass_patch = [gold_patch, purple_patch, red_patch, blue_patch, blue_patch1, blue_patch2]

time_span = np.linspace(0,36505, int(36500/5))
xticks_plot = time_span[::730].astype(int)
xticks_label = np.linspace(0,100,100)[::10].astype(int)

for seed_sim in seed_sim_list:
    for c_n in cn_list:
        details_subfolder = filestring + '_carbon_' + str(c_n) + '_'+str(seed_sim) + '_ip_0'
        simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
        figures_dir = os.path.join(project_dir, "figures", details_subfolder)
        hr = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'r')
        for dom_init in init_dom_list:
            fig1, ax1 = plt.subplots(figsize = (6,4))
            ax1.set_title(str(seed_sim) + "  " + str(c_n) + "  "  + str(dom_init))
            #fig2, ax2 = plt.subplots(figsize = (6,4))
            #ax2.set_title("Max_OS " + str(seed_sim) + "  " + str(c_n) + "  "  + str(dom_init))
            for c_b_r in bio_n_series:
                base_id = "b_1_all_/bio_n_" + str(c_b_r) + "/dom_initial_" + str(dom_init) + "/seed_" + str(seed_sim)
                base_data = hr[base_id]
                x = base_data['solution']
                paras = base_data["parameters"]
                C = x["dom"]
                C_sum = np.sum(x["dom"], axis = 1)
                max_idx = np.array(np.argmax(C, axis = 1))
                mean_wt_os = np.sum(np.array(C)*np.array(paras["oxidation_state"]).T, axis = 1)/C_sum
                ax1.plot(time_span, mean_wt_os, '-', color = cmap_bio(c_b_r*8), label = c_b_r)
                ax1.set_xticks(xticks_plot)
                ax1.set_xticklabels(xticks_label)
                #max_os = []
                #for i in max_idx:
                #    max_os.append(paras["oxidation_state"][i])
                #ax2.plot(time_span, max_os, '-', color = cmap_bio(c_b_r*8), label = c_b_r)
                #ax2.set_xticks(xticks_plot)
                #ax2.set_xticklabels(xticks_label)
            sm = plt.cm.ScalarMappable(cmap=cmap_bio)
            sm.set_array([])
            cbar = fig1.colorbar(sm)
            cbar.ax.get_yaxis().set_ticks([])
            for j, lab in enumerate(bio_n_series):
                cbar.ax.text(2.0, j/7, lab, ha='center', va='center')
            cbar.ax.get_yaxis().labelpad = 35
            cbar.ax.set_ylabel('biomass groups', rotation=270)
            ax1.set_ylabel ("Weighted mean oxidation state")
            ax1.set_xlabel ("Time [T]")
            fig1.tight_layout()
            fig1.savefig(os.path.join(figures_dir, "Baseline_Time_series_NOSC_"+str(dom_init)+".png"))
            #cbar = fig2.colorbar(sm)
            #cbar.ax.get_yaxis().set_ticks([])
            #for j, lab in enumerate(bio_n_series):
            #    cbar.ax.text(2.0, j/7, lab, ha='center', va='center')
            #cbar.ax.get_yaxis().labelpad = 35
            #cbar.ax.set_ylabel('biomass groups', rotation=270)
            #ax2.set_ylabel ("Max oxidation state")
            #ax2.set_xlabel ("Time [T]")
            #fig2.tight_layout()
            #fig2.savefig(os.path.join(figures_dir, "Baseline_Time_series_Max_OSC_"+str(dom_init)+".png"))
        hr.close
