#%%
import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
all_data = pd.read_csv(os.path.join(project_dir,"simulation_results_temporal_decay_const.csv"))
print(all_data.shape)
print(all_data.dtypes)
#%%
seeds= all_data.Seed.unique().tolist()
bios = all_data.biomass_species.unique().tolist()
carbs = all_data.carbon_species.unique().tolist()
carb_types = all_data.C_pool.unique().tolist()
docs = all_data.DOC_initial_int.unique().tolist()
sims = all_data.Sim_series.unique().tolist()
variances=all_data.Variance.unique().tolist()
T_cols = list(x for x in all_data.columns if 'T' in x)
old_t = 'T0'
for t_col in T_cols:
    if t_col == 'T0':
        all_data['dec0'] = 1/(all_data['T0']*5)
    else:
        #all_data['dec'+t_col[1:]] = 1/(all_data[t_col]*5) #execute when all calcs are for subsequent calcs
        all_data['T_diff'+t_col[1:]] = abs((all_data[t_col]-all_data[old_t])) #execute when all calcs are wrt initial conditions
        all_data['dec'+t_col[1:]] = 1/(all_data['T_diff'+t_col[1:]]*5) #execute when all calcs are wrt initial conditions
    all_data['dec_ratio'+t_col[1:]]=all_data['dec'+t_col[1:]]/all_data['dec0']#compl = all_data.dropna(subset = ['ratio_t_50'])
    old_t = t_col
all_data.replace([np.inf, -np.inf, np.nan], 0, inplace=True)
dec_ratio_columns = list(x for x in all_data.columns if 'dec_ratio' in x) 
#all_data[dec_ratio_columns].clip(lower=0, inplace=True)
act_full = all_data[all_data.Sim_series=="b_1_all_"]
extr_full = act_full[(act_full['DOC_initial_int']==2000.)|(act_full['DOC_initial_int']==10000.)]
#%%
###--------------------------------------------------------
### COMPARISON OF DECAY CONSTANT OF DIFFERENT CARBON POOLS
###--------------------------------------------------------
row_plots = ['dec3']
col_plots = ['DOC', 'reduced_C', 'oxidized_C']
fig, axes = plt.subplots(len(row_plots),len(col_plots),sharex=True, sharey=True, figsize = (10,4))
ax = axes.flatten()
for i in list(range(len(col_plots))):
    subset = extr_full[extr_full.C_pool==col_plots[i]].reset_index()
    for j in list(range(len(row_plots))):
        axindx = j*len(col_plots) + i
        #print(axindx,subset[row_plots[j]].shape)
        g=sns.scatterplot(x=subset['FD_initial'],y=subset[row_plots[j]]/subset['T0'], hue = subset['DOC_initial_int'], ax=ax[axindx])
        #ax[axindx].scatter(x=subset['FD_initial'],y=subset[row_plots[j]])#, hue = subset['DOC_initial_int'], ax=ax[axindx])
        g.legend().remove()
        ax[axindx].set_xlabel("")
fig.supxlabel("Functional diversity (Variance)", fontsize = 14)
ax[0].set_title("DOC", fontsize = 14)
ax[1].set_title("reduced C", fontsize = 14)
ax[2].set_title("Necromass", fontsize = 14)
ax[0].set_ylabel("Decay constant", fontsize = 14)
plt.xscale("log")
#plt.ylim(bottom=-1)
for a in ax[:]:
    a.axhline(y=0.0, c='maroon', linestyle='dashed')
handles,labels=ax[axindx].get_legend_handles_labels()
plt.figlegend(handles,labels,title = 'C availability', fontsize = 12, title_fontsize = 12, bbox_to_anchor=(0.85,-0.1), ncol=5, loc = "lower right", borderpad=0.)

#%%
###--------------------------------------------------------
### COMPARISON OF DECAY CONSTANT OF DIFFERENT CARBON POOLS TEMPORALLY
###--------------------------------------------------------
cols = ["DOC","reduced_C", "oxidized_C"]
select_columns = np.arange(769,770)#np.random.randint(0,1440,200)#[0,4,8,50,70,101,300,420,500,711,899,1000]#[10,40,80,1010,400,1400]#[0,4,8,101,420,1000]
#labelist=[90.0,81.0,72.9,65.6,59.0,53.1,47.8,43.0,38.7,34.9,31.4,28.2,25.4,22.9,20.6,18.5,16.7,15.0,13.5,12.2,10.9,9.8,8.9,8.0,7.2,6.5,5.8,5.2,4.7,4.2,3.8,3.4,3.1,2.8,2.5,2.3,2.0,1.8,1.6,1.4]
labelist=[90,80,70,60,50,40,30,20,10]
fig, axes = plt.subplots(nrows=5, ncols=len(cols), figsize=(10,10), sharex=True,sharey=True)
for iidx, i in enumerate(docs):
    for jidx, j in enumerate(cols):
        ax = axes[iidx][jidx]
        doc_data = act_full[(act_full['C_pool']==j)&(act_full['DOC_initial_int']==i)].reset_index()
        doc_t_data = doc_data[dec_ratio_columns].T
        ax.plot(doc_t_data.loc[:,select_columns])
        if jidx==0:
            ax.set_ylabel(i)
        if iidx==0:
            ax.set_title(j)
labels = list(str(int(x)) for x in labelist[::8])
plt.ylim(top=2)
plt.xticks(np.arange(0, 9, step=1), labelist)
fig.supylabel("Normalized decay constant", fontsize = 16)
fig.supxlabel("% remaining C", fontsize = 16)
#%%
###--------------------------------------------------------
### CONSOLIDATING RESULTS OF TEMPORAL VARIATION OF DECAY CONSTANT
###--------------------------------------------------------
group_all_data = all_data.groupby(['carbon_species','biomass_species','Seed','Sim_series','DOC_initial_int', 'C_pool'])
def f_argmax(x):
    #dec_ratio_vec = np.asarray(x[dec_ratio_columns].T)
    #dec_where_max = np.argmax(dec_ratio_vec)
    dec_remaining_c = labelist[np.argmax(np.asarray(x[dec_ratio_columns].T))]
    return dec_remaining_c

def f_argzero(x):
    #dec_ratio_vec = np.asarray(x[dec_ratio_columns].T)
    #dec_where_0 = np.argmin(dec_ratio_vec==0)[0][0]
    dec_remaining_c = labelist[np.argmin(np.asarray(x[dec_ratio_columns].T))]
    return dec_remaining_c

decay_max_series = group_all_data.apply(f_argmax)
decay_stop_series = group_all_data.apply(f_argzero)
#%%
def f_argmax(x):
    docwheremax = np.argmax(np.asarray(x.T))
    if docwheremax.size>0:
        dec_remaining_c = labelist[docwheremax]
    else:
        dec_remaining_c = 100
        docwheremax = -1
    return dec_remaining_c, docwheremax

def f_argzero(x):
    docwhere0 = np.argwhere(np.asarray(x.T)==0.)
    if docwhere0.size>0:
        docwhere0idx=docwhere0[0][0]
        dec_remaining_c = labelist[docwhere0idx]
    else:
        dec_remaining_c = 1
        docwhere0idx = -1
    return dec_remaining_c, docwhere0idx

data = act_full

for c in carbs:
    for b in bios:
        for s in seeds:
            for d in [2000.,10000.]:
                for sim in ['b_1_all_']:
                    for v in variances:
                        for ctype in carb_types:
                            subset = data[(data['C_pool']==ctype)&(data['carbon_species']==c)&(data['biomass_species']==b)&(data['Seed']==s)&(data['Sim_series']==sim)&(data['DOC_initial_int']==d)&(data['Variance']==v)][dec_ratio_columns]
                            decay_max, max_point = f_argmax(subset)
                            if max_point==-1:
                                decay_stop = 10
                            else:
                                decay_stop,stop_point = f_argzero(subset)
                            data.loc[(data['C_pool']==ctype)&(data['carbon_species']==c)&(data['biomass_species']==b)&(data['Seed']==s)&(data['Sim_series']==sim)&(data['DOC_initial_int']==d)&(data['Variance']==v),'decay_max']=decay_max
                            data.loc[(data['C_pool']==ctype)&(data['carbon_species']==c)&(data['biomass_species']==b)&(data['Seed']==s)&(data['Sim_series']==sim)&(data['DOC_initial_int']==d)&(data['Variance']==v),'decay_stop']=decay_stop

#%%
extr_full=data
extr_full['cbd']=extr_full.carbon_species*extr_full.biomass_species*extr_full.DOC_initial_int
extr_full.to_csv('C:/Users/swkh9804/Documents/random.csv', index=False)
#%%
extr_full = pd.read_csv('C:/Users/swkh9804/Documents/random.csv')
no_necro = extr_full[extr_full!='necromass']
#%%
no_necro['logdec0']=np.log10(extr_full.dec0)
g = sns.FacetGrid(col='C_pool', row = 'Variance', hue = 'logdec0', data = no_necro, palette="flare")
g.map(sns.scatterplot, 'decay_max', 'decay_stop')
g.set(ylim=(0, 100), xlim=(0,100), xlabel="%C for peak decomposition", ylabel="%C for C decomposition ending")
g.add_legend()
#%%
sns.scatterplot(x='decay_max', y='decay_stop',data = no_necro, style = 'Variance', hue = 'DOC_initial_int', alpha = 0.7)
plt.xlabel("%C at peak decomposition")
plt.ylabel("%C persistent")
plt.title("Fully active communities")
#%%
subset = extr_full[(extr_full['C_pool']=='oxidized_C')&(extr_full['Variance']==1.5)]
sns.jointplot(data=subset, x = 'decay_max', y = 'decay_stop', hue = 'carbon_species')
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
#Derive ratio of decay constants with respect to first decay constant
ratio_df = act_full
ratio_df['dec_20_ratio'] = ratio_df.Decay_constant_20/ratio_df.Decay_constant_10
ratio_df['dec_30_ratio'] = ratio_df.Decay_constant_30/ratio_df.Decay_constant_10
ratio_df['dec_40_ratio'] = ratio_df.Decay_constant_40/ratio_df.Decay_constant_10
ratio_df['dec_50_ratio'] = ratio_df.Decay_constant_50/ratio_df.Decay_constant_10
ratio_df['dec_60_ratio'] = ratio_df.Decay_constant_60/ratio_df.Decay_constant_10
#%%
#%%
pair_df = ratio_df[ratio_df.C_pool=='DOC'][['Variance','DOC_initial_int','Decay_constant_10','dec_20_ratio','dec_30_ratio','dec_40_ratio','dec_50_ratio','dec_60_ratio']]
sns.pairplot(data=pair_df, hue = 'DOC_initial_int', corner =True, dropna=True)
plt.xscale("log")
#%%
pair_df = ratio_df[ratio_df.C_pool=='reduced_C'][['Variance','DOC_initial_int','Decay_constant_10','dec_20_ratio','dec_30_ratio','dec_40_ratio','dec_50_ratio','dec_60_ratio']]
sns.pairplot(data=pair_df, hue = 'DOC_initial_int', corner =True, dropna=True)
#%%
pair_df = ratio_df[ratio_df.C_pool=='oxidized_C'][['Variance','DOC_initial_int','Decay_constant_10','dec_20_ratio','dec_30_ratio','dec_40_ratio','dec_50_ratio','dec_60_ratio']]
sns.pairplot(data=pair_df, hue = 'DOC_initial_int', corner =True, dropna=True)
#%%
#%%
fig, axes = plt.subplots(6,1,sharex=True,figsize = (4,8))
ax = axes.flatten()
#sns.histplot(data=ratio_df, binwidth = 0.001, x = 'Decay_constant_10', hue = 'C_pool', ax = ax[0])
sns.histplot(data=ratio_df, binwidth = 0.5, x = 'dec_20_ratio', hue = 'C_pool', ax = ax[1])
sns.histplot(data=ratio_df, binwidth = 0.5, x = 'dec_30_ratio', hue = 'C_pool', ax = ax[2])
sns.histplot(data=ratio_df, binwidth = 0.5, x = 'dec_40_ratio', hue = 'C_pool', ax = ax[3])
sns.histplot(data=ratio_df, binwidth = 0.5, x = 'dec_50_ratio', hue = 'C_pool', ax = ax[4])
sns.histplot(data=ratio_df,binwidth = 0.5, x = 'dec_60_ratio', hue = 'C_pool', ax = ax[5])
handles,labels=ax[-1].get_legend_handles_labels()
#plt.xscale("log")
plt.figlegend(handles,labels,title = 'C pool', fontsize = 12, title_fontsize = 12, bbox_to_anchor=(0.85,-0.1), ncol=5, loc = "lower right", borderpad=0.)
plt.xlim(right=10)
for a in ax[:]:
    a.legend().remove()
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