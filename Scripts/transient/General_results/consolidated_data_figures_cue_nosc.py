#%%
import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

import matplotlib.lines as mlines

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
sim_suffixes = ["_0_01", "_0_1", "_0_5", "", "_1_5x"]
sim_suffixes_var = [0.01, 0.1, 0.5, 1, 1.5]

files = []
for p,a in zip(sim_suffixes, sim_suffixes_var):
    filename = os.path.join(project_dir, "gen_spec_lognorm" + p, "results", "competition_adaptation_carbon_cue_wremaining_c_data.pkl")
    data = pd.read_pickle(filename)
    var_arr = np.zeros((data.shape[0],))+a
    var_ser = pd.Series(var_arr, copy=False,name = "Variance")
    data_var = data.join(var_ser)
    files.append(data_var)
all_data = pd.concat(files)
print(all_data.columns)
all_data['DOC_initial_int'] = round(all_data.DOC_initial, -3)
all_data=all_data.dropna(subset=['%C'])
all_data = all_data.replace(np.nan, 0.)
#%%
#compl = all_data.dropna(subset = ['ratio_t_50'])
act_full = all_data[all_data.Sim_series=="b_1_all_"]
act_off = all_data[all_data.Sim_series!="b_1_all_"]
extr_full = act_full[(act_full['DOC_initial_int']==2000.)|(act_full['DOC_initial_int']==10000.)]
extr_off = act_off[(act_off['DOC_initial_int']==2000.)|(act_off['DOC_initial_int']==10000.)]
#compl = compl.sort_values(by = 'active_H_c_connections')
#compl_act_off = compl[compl.activity<100]
#%%
sns.displot(x='CUE', data=act_full, hue = 'Variance', row = 'DOC_initial', col='%C', kind='kde')
#%%
docs = [1000.,2000.,5000.,10000.,15000.]
cols = [3,6,12,18]#sim_suffixes_var
#select_columns = act_full.index[act_full['FD_initial_cut_n']==fd_bins[2]].tolist()[:5]
#select_columns = np.arange(769,770)#np.random.randint(0,1440,200)#[0,4,8,50,70,101,300,420,500,711,899,1000]#[10,40,80,1010,400,1400]#[0,4,8,101,420,1000]
#labelist=[90.0,81.0,72.9,65.6,59.0,53.1,47.8,43.0,38.7,34.9,31.4,28.2,25.4,22.9,20.6,18.5,16.7,15.0,13.5,12.2,10.9,9.8,8.9,8.0,7.2,6.5,5.8,5.2,4.7,4.2,3.8,3.4,3.1,2.8,2.5,2.3,2.0,1.8,1.6,1.4]
labelist=[90,80,70,60,50,40,30,20,10]
fig, axes = plt.subplots(nrows=5, ncols=len(cols), figsize=(10,10), sharex=True,sharey='row')
for iidx, i in enumerate(docs):
    for jidx, j in enumerate(cols):
        ax = axes[iidx][jidx]
        cue_data = act_full[(act_full['DOC_initial_int']==i)&(act_full['carbon_species']==j)].reset_index()
        ax.scatter(cue_data['%C'],cue_data['CUE'])
        if jidx==0:
            ax.set_ylabel(i)
        if iidx==0:
            ax.set_title(j)
labels = list(str(int(x)) for x in labelist[::8])
plt.ylim(bottom = 0)
#plt.xticks(np.arange(0, 9, step=1), labelist)
fig.supylabel("Normalized decay constant", fontsize = 16)
fig.supxlabel("% remaining C", fontsize = 16)
#%%
h = sns.violinplot(x = "carbon_species", y = "NOSC_initial", data = act_full)
plt.xlabel("Size of carbon pool")
plt.ylabel("Initial nature of carbon pool (NOSC)")
#%%
selected_data = extr_full
hue_var = 'Variance'#'carbon_species'
fig, axes = plt.subplots(2,4,sharey='row',sharex='row',figsize=(10,5))
sns.kdeplot('NOSC_initial', hue = hue_var, data = selected_data, ax = axes[0][0])
sns.kdeplot('NOSC_timscale', hue = hue_var, data = selected_data, ax = axes[0][1])
sns.kdeplot('NOSC_maxbio', hue = hue_var, data = selected_data, ax = axes[0][2])
sns.kdeplot('NOSC_final', hue = hue_var, data = selected_data, ax = axes[0][3])
sns.kdeplot('CUE_initial', hue = hue_var, data = selected_data, ax = axes[1][0])
sns.kdeplot('CUE_timscale', hue = hue_var, data = selected_data, ax = axes[1][1])
sns.kdeplot('CUE_maxbio', hue = hue_var, data = selected_data, ax = axes[1][2])
sns.kdeplot('CUE_final', hue = hue_var, data = selected_data, ax = axes[1][3])
axes[0,0].set_ylabel("NOSC")
axes[-1,0].set_ylabel("CUE")
axes[0,0].set_title("Initial\ntime step")
axes[0,1].set_title("Characteristic\ntime scale")
axes[0,2].set_title("At maximum\nbiomass")
axes[0,3].set_title("100 years")
for a in list(range(8)):
    if a != 7:
        axes.flatten()[a].legend().set_visible(False)
    axes.flatten()[a].set_xlabel("")

#%%
hue_var = 'Variance'#'DOC_initial'
sns.jointplot(y = 'CUE_timscale', x = 'NOSC_timscale', hue = hue_var, data = selected_data)
plt.ylim((-1,1))
#%%
hue_var = 'DOC_initial_int'
style_var = extr_full.carbon_species*extr_full.biomass_species
selected_data = extr_full#[extr_full.DOC_initial_int==10000]
x_var = 'NOSC_initial'
transparency = 0.5
fig, axes = plt.subplots(1,3,sharex=True,sharey=True,figsize=(10,4))
sns.scatterplot(x = x_var, y = 'NOSC_timscale', hue = hue_var, size = style_var, alpha = transparency, data = selected_data, ax = axes[0])
sns.scatterplot(x = x_var, y = 'NOSC_maxbio', hue = hue_var, size = style_var, alpha = transparency,data = selected_data, ax=axes[1])
sns.scatterplot(x = x_var, y = 'NOSC_final', hue = hue_var, size = style_var, alpha = transparency,data = selected_data, ax = axes[2])
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
        a.legend(bbox_to_anchor=(0.25,-0.2), ncol = 3)
#%%
hue_var = 'DOC_initial_int'
size_var = 'S_initial_int'
selected_data = extr_full#[extr_full.DOC_initial_int==10000]
x_var = 'CUE_initial'
x_var1 = 'NOSC_timscale'
x_var2 = 'NOSC_maxbio'
x_var3 = 'NOSC_final'
transparency = 0.7
fig, axes = plt.subplots(1,3,sharex=True,sharey=True,figsize=(10,4))
sns.scatterplot(x = x_var1, y = 'CUE_timscale', size = size_var, hue = hue_var, alpha = transparency, data = selected_data, ax = axes[0])
sns.scatterplot(x = x_var2, y = 'CUE_maxbio', size = size_var, hue = hue_var, alpha = transparency, data = selected_data, ax=axes[1])
sns.scatterplot(x = x_var3, y = 'CUE_final', size = size_var, hue = hue_var, alpha = transparency, data = selected_data, ax = axes[2])
axes[0].set_title ("Characteristic\ntime scale")
axes[1].set_title("Peak biomass")
axes[2].set_title("100 years")
axes[0].set_ylabel("CUE")
axes[0].set_xlim(-0.5,0.5)
for a in axes:
    a.set_xlabel("NOSC")#("Initial CUE")
    if a!=axes[2]:
        a.get_legend().remove()
    else:
        a.legend(bbox_to_anchor=(0.25,-0.2), ncol = 5, title = "H")
#%%
g = sns.FacetGrid(extr_full, row="DOC_initial_int", col = 'S_initial_int', hue='Variance')
g.map(sns.scatterplot,'NOSC_initial','CUE_initial')
g.add_legend()
#%%
g = sns.FacetGrid(extr_full, row="DOC_initial_int", col = 'S_initial_int', hue='Variance')
g.map(sns.scatterplot,'NOSC_initial','CUE_timscale')
g.add_legend()
#%%
g = sns.FacetGrid(extr_full, row="DOC_initial_int", col = 'S_initial_int', hue='Variance')
g.map(sns.scatterplot,'NOSC_initial','CUE_maxbio')
g.add_legend()
#%%
g = sns.FacetGrid(extr_full, row="DOC_initial_int", col = 'S_initial_int', hue = 'Variance')
g.map(sns.scatterplot,'NOSC_initial','CUE_final')
g.add_legend()
#%%
g = sns.pairplot(extr_full[['NOSC_initial', 'NOSC_timscale', 'NOSC_maxbio', 'NOSC_final','S_initial_int']],
hue="S_initial_int", diag_kind="kde")##, height=2.5)
#%%
g = sns.pairplot(extr_full[['CUE_initial', 'CUE_timscale', 'CUE_maxbio', 'CUE_final', 'S_initial_int']],
hue="S_initial_int", diag_kind="kde")##, height=2.5)
#%%
g = sns.pairplot(extr_full[['NOSC_initial', 'CUE_initial', 'CUE_timscale', 'CUE_maxbio', 'CUE_final', 'S_initial_int']],
hue="S_initial_int", diag_kind="kde")##, height=2.5)

#%%
sns.jointplot(y = 'CUE_maxbio', x = 'NOSC_maxbio', hue = 'S_initial', data = extr_full)
plt.ylim(bottom=-0.05)
#%%
sns.scatterplot(y = 'CUE_timscale', x = 'S_initial', hue = 'DOC_initial_int',size = 'NOSC_initial', data = extr_full)
plt.legend(bbox_to_anchor=(1,1))
#%%
eqline = np.sort(all_data["NOSC_initial"])
sns.scatterplot(x = "NOSC_initial", y = "NOSC_final", hue = 'DOC_initial', size = 'S_initial_int',style ="Variance", data = extr_full[extr_full.Variance>=0.5])#all_data)
plt.plot(eqline, eqline, 'r-')
plt.ylim((-0.5,0.5))
plt.legend(bbox_to_anchor=(1,1))
#%%
eqline = np.sort(all_data["NOSC_initial"])
sns.scatterplot(x = "NOSC_initial", y = "NOSC_final", hue = 'DOC_initial', size = 'S_initial_int',style ="Variance", data = extr_full[extr_full.FD_initial<0.0001])#all_data)
plt.plot(eqline, eqline, 'r-')
plt.ylim((-0.5,0.5))
plt.legend(bbox_to_anchor=(1,1))

#%%
eqline = np.sort(all_data["NOSC_initial"])
sns.scatterplot(x = "NOSC_initial", y = "NOSC_maxbio", hue = 'DOC_initial', style ="Variance", data = all_data)
plt.plot(eqline, eqline, 'r-')
plt.legend(bbox_to_anchor=(1,1))
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