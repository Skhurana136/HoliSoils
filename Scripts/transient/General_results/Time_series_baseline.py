## Import libraries
import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import sys

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient", "gen_spec_skew")
filestring = sys.argv[1]
ip = 0
styles = {4:"darkgoldenrod", 8:"orange", 12:"indianred", 20:"purple", 28:"steelblue", 32: "darkblue"}
cmap_bio = plt.cm.get_cmap('Blues')

seed_sim_list = [610229235, 983307757, 643338060, 714504443, 277077803, 898393994, 420,13012022,13061989]
cn_list = [3,6,12,18]
bio_n_series = [4,8,12,20,28,32]
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
            fig, ax1 = plt.subplots(figsize = (6,4))
            ax2 = ax1.twinx()
            for c_b_r in bio_n_series:
                base_id = "b_1_all_/bio_n_" + str(c_b_r) + "/dom_initial_" + str(dom_init) + "/seed_" + str(seed_sim)
                base_data = hr[base_id]
                x = base_data['solution']
                C_sum_base = (1- (np.sum(np.asarray(x['dom']), axis = 1)/dom_init))*100
                B_sum_base = np.sum(np.asarray(x['biomass']), axis = 1)
                ax1.set_xlabel ("Time [T]")
                ax1.set_ylabel ("Carbon consumed (%)")
                ax2.set_ylabel ("Biomass [N $L^{-3}$]")
                ax1.plot(time_span, C_sum_base, '-', color = cmap_bio(c_b_r*8), label = c_b_r)
                ax2.plot(time_span, B_sum_base, '--', color = cmap_bio(c_b_r*8))
                ax1.set_ylim(-10,100)
                ax1.set_xticks(xticks_plot)
                ax1.set_xticklabels(xticks_label)
                ax2.set_xticks(xticks_plot)
                ax2.set_xticklabels(xticks_label)
                ax1.set_xlim(left = -10)
                ax2.set_xlim(left = -10)
            legend1 = ax1.legend(handles = linelist, title = "Component", bbox_to_anchor = (0.75,-0.5), ncol = 2)
            #ax1.add_artist(legend1)
            sm = plt.cm.ScalarMappable(cmap=cmap_bio)
            sm.set_array([])
            cbar = plt.colorbar(sm, pad = 0.2, orientation = 'horizontal')
            cbar.ax.get_xaxis().set_ticks([])
            for j, lab in enumerate(bio_n_series):
                cbar.ax.text(j/5, -0.5, lab, ha='center', va='center')
            cbar.ax.get_xaxis().labelpad = -25
            cbar.ax.set_xlabel('biomass groups')#, rotation = 180)
            #legend2 = ax1.legend(handles = biomass_patch, ncol = 3, title = "Biomass groups", bbox_to_anchor = (0.5, -0.18))
            plt.tight_layout()
            plt.savefig(os.path.join(figures_dir, "Baseline_Time_series_categorical_"+str(dom_init)+".png"))
            plt.close()
        hr.close
 
#for seed_sim in seed_sim_list:
#    for c_n in cn_list:
#        details_subfolder = 'v2_1c_adaptation_carbon_' + str(c_n) + '_'+str(seed_sim) + '_ip_0'
#       simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
#        figures_dir = os.path.join(project_dir, "figures", details_subfolder)
#        hr = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'r')
#        for dom_init in init_dom_list:
#            for c_b_r in bio_n_series:
#                fig, ax1 = plt.subplots(figsize = (6,4))
#                ax2 = ax1.twinx()
#                base_id = "b_1_all_/bio_n_" + str(c_b_r) + "/dom_initial_" + str(dom_init) + "/seed_" + str(seed_sim)
#                base_data = hr[base_id]
#                x = base_data['solution']
#                C = x['dom']
#                B = x['biomass']
#                ax1.set_xlabel ("Time [T]")
#                ax1.set_ylabel ("Carbon [N $L^{-3}$]")
#                ax2.set_ylabel ("Biomass [N $L^{-3}$]")
#                ax1.plot(time_span, C, '-')
#                ax2.plot(time_span, B, '--')
#                ax1.set_xticks(xticks_plot)
#                ax1.set_xticklabels(xticks_label)
#                ax2.set_xticks(xticks_plot)
#                ax2.set_xticklabels(xticks_label)
#                ax1.set_xlim(left = -10)
#                ax2.set_xlim(left = -10)
#                plt.tight_layout()
#                plt.savefig(os.path.join(figures_dir, "Baseline_Time_series_diff_c_bio_"+str(c_b_r)+"_"+str(dom_init)+".png"))
#                plt.close()
#        hr.close
 