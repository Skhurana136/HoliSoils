#%%
## Import libraries
import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

## LOAD RESULTS
#project_dir = os.path.join("C:/", "Users", "swkh9804", "Documents", "Project_data", "HoliSoils", "transient", sys.argv[1])
project_dir = "C:/Users/swkh9804/Documents/Scripts/HoliSoils/Data"
figures_dir = os.path.join(project_dir, "Figures")
data_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
filestring = "competition_adaptation"
ip = 0
#styles = {4:"darkgoldenrod", 8:"orange", 12:"indianred", 20:"purple", 28:"steelblue", 32: "darkblue"}
cmap_bio = plt.cm.get_cmap('Oranges')
Seed_select = '898393994'
carbon_n = 6
bio_low=4
bio_high=20
doc_low = 2000
doc_high = 10000

simulation_subfolder = filestring + '_carbon_' + str(carbon_n) + '_'+str(Seed_select) + '_ip_0'
similar_dir = os.path.join(data_dir, 'gen_spec_lognorm_0_01', "simulations", simulation_subfolder)
dissimilar_dir = os.path.join(data_dir, 'gen_spec_lognorm', "simulations", simulation_subfolder)

similar_hr = h5py.File(os.path.join(similar_dir,"simulations.h5"), mode = 'r')
dissimilar_hr = h5py.File(os.path.join(dissimilar_dir,"simulations.h5"), mode = 'r')

bl_docl = "b_1_all_/bio_n_" + str(bio_low) + "/dom_initial_" + str(doc_low) + "/seed_" + str(Seed_select) + "/solution"
bl_doch = "b_1_all_/bio_n_" + str(bio_low) + "/dom_initial_" + str(doc_high) + "/seed_" + str(Seed_select) + "/solution"
bh_docl = "b_1_all_/bio_n_" + str(bio_high) + "/dom_initial_" + str(doc_low) + "/seed_" + str(Seed_select) + "/solution"
bh_doch = "b_1_all_/bio_n_" + str(bio_high) + "/dom_initial_" + str(doc_high) + "/seed_" + str(Seed_select) + "/solution"
sim_bl_docl = similar_hr[bl_docl]
sim_bl_doch = similar_hr[bl_doch]
dissim_bl_docl = dissimilar_hr[bl_docl]
dissim_bl_doch = dissimilar_hr[bl_doch]
sim_bh_docl = similar_hr[bh_docl]
sim_bh_doch = similar_hr[bh_doch]
dissim_bh_docl = dissimilar_hr[bh_docl]
dissim_bh_doch = dissimilar_hr[bh_doch]

#%%
lowvar_lowbn = mlines.Line2D([], [], linestyle = '-', color = "black", marker = '.', label='Low $V_{b}$ and low Nb')
highvar_lowbn = mlines.Line2D([], [], linestyle = 'dotted', color = "black", marker = '.', label='High $V_{b}$ and low Nb')
lowvar_highbn = mlines.Line2D([], [], linestyle = '-', color = "black", marker = 'x', label='Low $V_{b}$ and high Nb')
highvar_highbn = mlines.Line2D([], [], linestyle = 'dotted', color = "black", marker = 'x', label='High $V_{b}$ and high Nb')

linelist = [lowvar_lowbn, highvar_lowbn, lowvar_highbn, highvar_highbn]
var_styles= ['-', 'dotted']
bio_styles =['.','x']

orange_patch = mpatches.Patch(color=cmap_bio(doc_low/5000), label= '$C_{0}$ poor')
darkorange_patch = mpatches.Patch(color=cmap_bio(doc_high/5000), label= '$C_{0}$ rich')
patchlist = [orange_patch, darkorange_patch]

#%%
time_span = np.linspace(0,36505, int(36500/5))
xticks_plot = time_span[::730].astype(int)
xticks_label = np.linspace(0,100,100)[::10].astype(int)


fig,axes = plt.subplots(2,1, figsize = (5,6), sharex = True)
ax = axes.flatten()

for x,doc_lev,v in zip([sim_bl_docl, sim_bl_doch, dissim_bl_docl,dissim_bl_doch],[doc_low,doc_high, doc_low, doc_high], [0,0,1,1]):
    bn=bio_low
    C_sum = (1- (np.sum(np.asarray(x['dom']), axis = 1)/doc_lev))*100
    B_sum = np.sum(np.asarray(x['biomass']), axis = 1)
    ax[0].plot(time_span[::500], B_sum[::500], color = cmap_bio(doc_lev/5000), linestyle = var_styles[v], marker = bio_styles[0])
    ax[1].plot(time_span[::500], C_sum[::500], color = cmap_bio(doc_lev/5000), linestyle = var_styles[v],  marker = bio_styles[0])
    
for x,doc_lev,v in zip([sim_bh_docl,sim_bh_doch,dissim_bh_docl,dissim_bh_doch],[doc_low,doc_high, doc_low, doc_high], [0,0,1,1]):
    bn=bio_high
    C_sum = (1- (np.sum(np.asarray(x['dom']), axis = 1)/doc_lev))*100
    B_sum = np.sum(np.asarray(x['biomass']), axis = 1)
    ax[0].plot(time_span[::500], B_sum[::500], color = cmap_bio(doc_lev/5000), linestyle = var_styles[v], marker = bio_styles[1])
    ax[1].plot(time_span[::500], C_sum[::500], color = cmap_bio(doc_lev/5000), linestyle = var_styles[v],  marker = bio_styles[1])

ax[0].set_title('A', fontsize = 16, loc='left')
ax[1].set_title('B', fontsize = 16, loc='left')

ax[0].set_ylabel ("Biomass [N $L^{-3}$]")
ax[1].set_ylabel ("Carbon consumed (%)")
ax[1].set_xlabel ("Time [T]")
ax[0].set_ylim(bottom = 0)
ax[1].set_ylim(0,110)
ax[1].set_xticks(xticks_plot)
ax[1].set_xticklabels(xticks_label)
ax[1].set_xlim(left = 0)
legend1 = ax[0].legend(handles = linelist, bbox_to_anchor = (1.0,-1.45), ncol = 2, title = "Community characteristics", frameon=False)
legend2 = ax[1].legend(handles = patchlist, bbox_to_anchor = (0.75,-0.7), ncol = 2, title = "C availability", frameon = False)
plt.savefig(os.path.join(figures_dir, "Fig_XX_Time_series_concentration.png"), bbox_inches='tight', dpi = 300)