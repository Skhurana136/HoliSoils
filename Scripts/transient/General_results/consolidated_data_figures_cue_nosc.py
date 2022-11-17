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
all_data = pd.read_csv(os.path.join(project_dir,"simulation_results_diversity_parameters.csv"))
print(all_data.shape)
print(all_data.dtypes)
#%%
compl = all_data.dropna(subset = ['ratio_t_50'])
act_full = all_data[all_data.Sim_series=="b_1_all_"]
act_off = all_data[all_data.Sim_series!="b_1_all_"]
extr_full = act_full[(act_full['DOC_initial_int']==2000.)|(act_full['DOC_initial_int']==10000.)]
extr_off = act_off[(act_off['DOC_initial_int']==2000.)|(act_off['DOC_initial_int']==10000.)]
compl = compl.sort_values(by = 'active_H_c_connections')
compl_act_off = compl[compl.activity<100]

#%%
h = sns.violinplot(x = "carbon_species", y = "NOSC_initial", data = act_full)
plt.xlabel("Size of carbon pool")
plt.ylabel("Initial nature of carbon pool")
#%%
sns.jointplot(y = 'CUE_timscale', x = 'NOSC_timscale', hue = 'DOC_initial', data = extr_full)
plt.ylim((-1,1))
#%%
fig, axes = plt.subplots(1,3,sharex=True,sharey=True,figsize=(10,4))
sns.scatterplot(x = 'NOSC_initial', y = 'NOSC_timscale', hue = 'S_initial_int', style = 'DOC_initial_int', data = extr_full[extr_full.DOC_initial_int==10000], ax = axes[0])
sns.scatterplot(x = 'NOSC_initial', y = 'NOSC_maxbio', hue = 'S_initial_int', style = 'DOC_initial_int', data = extr_full[extr_full.DOC_initial_int==10000], ax=axes[1])
sns.scatterplot(x = 'NOSC_initial', y = 'NOSC_final', hue = 'S_initial_int', style = 'DOC_initial_int', data = extr_full[extr_full.DOC_initial_int==10000], ax = axes[2])
axes[0].set_title ("Characteristic\ntime scale")
axes[1].set_title("Peak biomass")
axes[2].set_title("100 years")
axes[0].set_ylabel("NOSC")
axes[0].set_ylim((-0.4,0.4))
axes[0].set_xlim((-0.4,0.4))
for a in axes:
    a.set_xlabel("Initial NOSC")
    if a!=axes[2]:
        a.get_legend().remove()
    else:
        a.legend(bbox_to_anchor=(0.25,-0.2), ncol = 2)
#%%
fig, axes = plt.subplots(1,3,sharex=True,sharey=True,figsize=(10,4))
sns.scatterplot(x = 'CUE_initial', y = 'CUE_timscale', size = 'S_initial_int', hue = 'DOC_initial_int', data = extr_full, ax = axes[0])
sns.scatterplot(x = 'CUE_initial', y = 'CUE_maxbio', size = 'S_initial_int', hue = 'DOC_initial_int', data = extr_full, ax=axes[1])
sns.scatterplot(x = 'CUE_initial', y = 'CUE_final', size = 'S_initial_int', hue = 'DOC_initial_int', data = extr_full, ax = axes[2])
axes[0].set_title ("Characteristic\ntime scale")
axes[1].set_title("Peak biomass")
axes[2].set_title("100 years")
axes[0].set_ylabel("CUE")
#axes[0].set_ylim((-1,1))
#axes[0].set_xlim((-1,1))
for a in axes:
    a.set_xlabel("Initial CUE")
    if a!=axes[2]:
        a.get_legend().remove()
    else:
        a.legend(bbox_to_anchor=(0.25,-0.2), ncol = 5, title = "H")
#%%
sns.jointplot(y = 'CUE_maxbio', x = 'NOSC_maxbio', hue = 'S_initial', data = extr_full)
plt.ylim(bottom=-0.05)
#%%
sns.scatterplot(y = 'CUE_timscale', x = 'S_initial', hue = 'DOC_initial_int',size = 'NOSC_initial', data = extr_full)
plt.legend(bbox_to_anchor=(1,1))
#%%
eqline = np.sort(all_data["NOSC_initial"])
sns.scatterplot(x = "NOSC_initial", y = "NOSC_maxbio", hue = 'DOC_initial', style ="Variance", data = all_data)
plt.plot(eqline, eqline, 'r-')
#%%
all_data['hue_var'] = all_data.biomass_species*all_data.activity/100#*all_data.carbon_species))
sns.scatterplot(x = "NOSC_initial", y = "NOSC_maxbio", hue = 'hue_var', size = 'DOC_initial', data = all_data)
plt.plot(eqline, eqline, 'r-')
#%%
sns.scatterplot('carbon_species', 'NOSC_initial', hue = "Seed",data = all_data)
#%%
subset = all_data[all_data.carbon_species==12]
eqline = np.sort(subset['NOSC_initial'])
plt.plot(eqline, eqline, 'r-')
sns.scatterplot(x = "NOSC_initial", y = "NOSC_maxbio", hue = 'FD_maxbio', data = subset)