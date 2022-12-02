#%%
import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

import matplotlib.lines as mlines

from scipy.optimize import curve_fit

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
all_data = pd.read_csv(os.path.join(project_dir,"simulation_results_temporal_decay_const.csv"))
print(all_data.shape)
print(all_data.dtypes)
#%%
#compl = all_data.dropna(subset = ['ratio_t_50'])
act_full = all_data[all_data.Sim_series=="b_1_all_"]
#act_off = all_data[all_data.Sim_series!="b_1_all_"]
extr_full = act_full[(act_full['DOC_initial_int']==2000.)|(act_full['DOC_initial_int']==10000.)]
#extr_off = act_off[(act_off['DOC_initial_int']==2000.)|(act_off['DOC_initial_int']==10000.)]
#compl = compl.sort_values(by = 'active_H_c_connections')
#compl_act_off = compl[compl.activity<100]
#no_na = all_data.dropna()
#lognona = no_na.apply(lambda x: np.log10(x) if np.issubdtype(x.dtype, np.number) else x)
#%%
###--------------------------------------------------------
### COMPARISON OF DECAY CONSTANT OF DIFFERENT CARBON POOLS
###--------------------------------------------------------
#row_plots = ['Decay_constant_10','Decay_constant_20','Decay_constant_30']
row_plots = ['Decay_constant_10','Decay_constant_20','Decay_constant_30','Decay_constant_40','Decay_constant_50','Decay_constant_60']
col_plots = ['DOC', 'reduced_C', 'oxidized_C']#'necromass',
fig, axes = plt.subplots(len(row_plots),len(col_plots),sharex=True,sharey =True, figsize = (10,8))
ax = axes.flatten()
for i in list(range(len(col_plots))):
    subset = extr_full[extr_full.C_pool==col_plots[i]].reset_index()
    for j in list(range(len(row_plots))):
        axindx = j*len(col_plots) + i
        #print(axindx,subset[row_plots[j]].shape)
        g=sns.scatterplot(x=subset['FD_initial'],y=subset[row_plots[j]], hue = subset['DOC_initial_int'], ax=ax[axindx])
        #ax[axindx].scatter(x=subset['FD_initial'],y=subset[row_plots[j]])#, hue = subset['DOC_initial_int'], ax=ax[axindx])
        g.legend().remove()
        ax[axindx].set_xlabel("")
fig.supxlabel("Functional diversity (Variance)", fontsize = 14)
ax[0].set_title("DOC", fontsize = 14)
ax[1].set_title("reduced C", fontsize = 14)
ax[2].set_title("oxidized C", fontsize = 14)
ax[0].set_ylabel("1st\n10% loss", fontsize = 14)
ax[len(col_plots)*1].set_ylabel("2nd\n10% loss", fontsize = 14)
ax[len(col_plots)*2].set_ylabel("3rd\n10% loss", fontsize = 14)
ax[len(col_plots)*3].set_ylabel("4th\n10% loss", fontsize = 14)
ax[len(col_plots)*4].set_ylabel("5th\n10% loss", fontsize = 14)
ax[len(col_plots)*5].set_ylabel("6th\n10% loss", fontsize = 14)
plt.xscale("log")
#plt.yscale("log")
handles,labels=ax[axindx].get_legend_handles_labels()
plt.figlegend(handles,labels,title = 'C availability', fontsize = 12, title_fontsize = 12, bbox_to_anchor=(0.85,-0.1), ncol=5, loc = "lower right", borderpad=0.)
#%%
fig, axes = plt.subplots(6,1,sharex=True,figsize = (4,8))
ax = axes.flatten()
sns.kdeplot(data=act_full, x = 'Decay_constant_10', hue = 'C_pool', ax = ax[0])
sns.kdeplot(data=act_full, x = 'Decay_constant_20', hue = 'C_pool', ax = ax[1])
sns.kdeplot(data=act_full, x = 'Decay_constant_30', hue = 'C_pool', ax = ax[2])
sns.kdeplot(data=act_full, x = 'Decay_constant_40', hue = 'C_pool', ax = ax[3])
sns.kdeplot(data=act_full, x = 'Decay_constant_50', hue = 'C_pool', ax = ax[4])
sns.kdeplot(data=act_full, x = 'Decay_constant_60', hue = 'C_pool', ax = ax[5])
handles,labels=ax[-1].get_legend_handles_labels()
#plt.xscale("log")
plt.figlegend(handles,labels,title = 'C pool', fontsize = 12, title_fontsize = 12, bbox_to_anchor=(0.85,-0.1), ncol=5, loc = "lower right", borderpad=0.)
for a in ax[:]:
    a.legend().remove()
#%%
pair_df = extr_full[extr_full.C_pool=='DOC'][['Variance','DOC_initial_int','Decay_constant_10','Decay_constant_20','Decay_constant_30','Decay_constant_40','Decay_constant_50','Decay_constant_60']]
sns.pairplot(data=pair_df, hue = 'DOC_initial_int', corner =True, dropna=True)
#%%
pair_df = extr_full[extr_full.C_pool=='reduced_C'][['Variance','DOC_initial_int','Decay_constant_10','Decay_constant_20','Decay_constant_30','Decay_constant_40','Decay_constant_50','Decay_constant_60']]
sns.pairplot(data=pair_df, hue = 'DOC_initial_int', corner =True, dropna=True)
#%%
pair_df = extr_full[extr_full.C_pool=='oxidized_C'][['Variance','DOC_initial_int','Decay_constant_10','Decay_constant_20','Decay_constant_30','Decay_constant_40','Decay_constant_50','Decay_constant_60']]
sns.pairplot(data=pair_df, hue = 'DOC_initial_int', corner =True, dropna=True)

#%%
plt.figure(figsize=(8,4))
h = sns.scatterplot(x = "FD_initial", y = "Decay_constant_10", size = "carbon_species", hue = "biomass_species", data = act_full)
h.set_ylabel("Decay constant", fontsize = 16)
h.set_xlabel("Initial functional\ndiversity: Variance", fontsize = 16)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xscale("log")
leg_handles, leg_labels = h.get_legend_handles_labels()
leg_labels[0]="#Microbial groups"
leg_labels[7]="#Carbon groups"
plt.legend(leg_handles, leg_labels,bbox_to_anchor=(1,1.1), fontsize = 14)
plt.tight_layout()
#plt.savefig(os.path.join(project_dir, "imposed_func_div.png"), dpi = 300)
#%%
plt.figure(figsize=(8,4))
h = sns.scatterplot(x = "FD_initial", y = "Decay_constant_20", size = "carbon_species", hue = "biomass_species", data = act_full)
h.set_ylabel("Decay constant", fontsize = 16)
h.set_xlabel("Initial functional\ndiversity: Variance", fontsize = 16)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xscale("log")
leg_handles, leg_labels = h.get_legend_handles_labels()
leg_labels[0]="#Microbial groups"
leg_labels[7]="#Carbon groups"
plt.legend(leg_handles, leg_labels,bbox_to_anchor=(1,1.1), fontsize = 14)
plt.tight_layout()
#%%
plt.figure(figsize=(8,4))
h = sns.scatterplot(x = "FD_initial", y = "Decay_constant_30", size = "carbon_species", hue = "biomass_species", data = act_full)
h.set_ylabel("Decay constant", fontsize = 16)
h.set_xlabel("Initial functional\ndiversity: Variance", fontsize = 16)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xscale("log")
leg_handles, leg_labels = h.get_legend_handles_labels()
leg_labels[0]="#Microbial groups"
leg_labels[7]="#Carbon groups"
plt.legend(leg_handles, leg_labels,bbox_to_anchor=(1,1.1), fontsize = 14)
plt.tight_layout()
#%%
plt.figure(figsize=(8,4))
h = sns.scatterplot(x = "Decay_constant_10", y = "Decay_constant_20", size = "carbon_species", hue = "biomass_species", data = extr_full)
h.set_ylabel("Decay constant", fontsize = 16)
#h.set_xlabel("Initial functional\ndiversity: Variance", fontsize = 16)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
#plt.xscale("log")
leg_handles, leg_labels = h.get_legend_handles_labels()
leg_labels[0]="#Microbial groups"
leg_labels[7]="#Carbon groups"
plt.legend(leg_handles, leg_labels,bbox_to_anchor=(1,1.1), fontsize = 14)
plt.tight_layout()
#%%
plt.figure(figsize=(8,4))
h = sns.scatterplot(x = "Decay_constant_10", y = "Decay_constant_30", size = "carbon_species", hue = "biomass_species", data = extr_full)
h.set_ylabel("Decay constant", fontsize = 16)
#h.set_xlabel("Initial functional\ndiversity: Variance", fontsize = 16)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
#plt.xscale("log")
leg_handles, leg_labels = h.get_legend_handles_labels()
leg_labels[0]="#Microbial groups"
leg_labels[7]="#Carbon groups"
plt.legend(leg_handles, leg_labels,bbox_to_anchor=(1,1.1), fontsize = 14)
plt.tight_layout()