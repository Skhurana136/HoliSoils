#%%
import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
all_data = pd.read_csv(os.path.join(project_dir,"simulation_results_temporal_initial_conditions_decay_const.csv"))
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
all_data.replace([np.inf, -np.inf, np.nan], 0, inplace=True)
dec_ratio_columns = list(x for x in all_data.columns if 'dec_ratio' in x) 
#all_data[dec_ratio_columns].clip(lower=0, inplace=True)
#%%
###--------------------------------------------------------
### BINNING ACCORDING TO FUNCTIONAL DIVERSITY
###--------------------------------------------------------
plt.figure()
sns.histplot(np.log10(all_data['FD_initial']), bins = 5)
plt.title("Distribution of initialized functional diversity")
all_data['FD_initial_cut_n'] = pd.cut(np.log10(all_data['FD_initial']), bins=[-13,-7,-5,-1])
x = all_data['FD_initial_cut_n'].astype(str)
plt.figure()
plt.title("Binned functional diversity")
sns.histplot(x)
plt.xticks(rotation=45)
fd_bins=all_data.FD_initial_cut_n.unique().tolist()
#%%
###--------------------------------------------------------
### CONSOLIDATING RESULTS OF TEMPORAL VARIATION OF DECAY CONSTANT
###--------------------------------------------------------
def f_argmax(x):
    docwheremax = np.argmax(np.asarray(x[dec_ratio_columns].T))
    if docwheremax.size>0:
        dec_remaining_c = labelist[docwheremax]
    else:
        dec_remaining_c = 100
        docwheremax = -1
    return dec_remaining_c#, docwheremax

def f_argzero(x):
    docwhere0 = np.argwhere(np.asarray(x[dec_ratio_columns].T)==0.)
    if docwhere0.size>0:
        docwhere0idx=docwhere0[0][0]
        dec_remaining_c = labelist[docwhere0idx]
    else:
        dec_remaining_c = 1
        docwhere0idx = -1
    return dec_remaining_c#, docwhere0idx

all_data['decay_max'] = all_data.apply(f_argmax, axis=1)
all_data['decay_stop'] = all_data.apply(f_argzero, axis=1)

#%%
###--------------------------------------------------------
### SUBSET DATA TO WORK WITH SMALLER DATASETS TO AMKE SENS
###--------------------------------------------------------
act_full = all_data[all_data.Sim_series=="b_1_all_"]
extr_full = act_full[(act_full['DOC_initial_int']==2000.)|(act_full['DOC_initial_int']==10000.)]

#%%
###--------------------------------------------------------
### COMPARISON OF DECAY CONSTANT OF DIFFERENT CARBON POOLS
###--------------------------------------------------------
row_plots = ['dec0']
col_plots = ['DOC', 'reduced_C', 'oxidized_C']
data = extr_full[extr_full.FD_initial_cut_n==fd_bins[1]]
fig, axes = plt.subplots(len(row_plots),len(col_plots),sharex=True, sharey=True, figsize = (10,4))
ax = axes.flatten()
for i in list(range(len(col_plots))):
    subset = data[data.C_pool==col_plots[i]].reset_index()
    for j in list(range(len(row_plots))):
        axindx = j*len(col_plots) + i
        #print(axindx,subset[row_plots[j]].shape)
        g=sns.scatterplot(x=subset['FD_initial'].astype(str),y=subset[row_plots[j]]/subset['T0'], hue = subset['DOC_initial_int'], ax=ax[axindx])
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
select_columns = act_full.index[act_full['FD_initial_cut_n']==fd_bins[2]].tolist()[:5]
#select_columns = np.arange(769,770)#np.random.randint(0,1440,200)#[0,4,8,50,70,101,300,420,500,711,899,1000]#[10,40,80,1010,400,1400]#[0,4,8,101,420,1000]
#labelist=[90.0,81.0,72.9,65.6,59.0,53.1,47.8,43.0,38.7,34.9,31.4,28.2,25.4,22.9,20.6,18.5,16.7,15.0,13.5,12.2,10.9,9.8,8.9,8.0,7.2,6.5,5.8,5.2,4.7,4.2,3.8,3.4,3.1,2.8,2.5,2.3,2.0,1.8,1.6,1.4]
labelist=[90,80,70,60,50,40,30,20,10]
fig, axes = plt.subplots(nrows=5, ncols=len(cols), figsize=(10,10), sharex=True,sharey='row')
for iidx, i in enumerate(docs):
    for jidx, j in enumerate(cols):
        ax = axes[iidx][jidx]
        doc_data = act_full[(act_full['C_pool']==j)&(act_full['DOC_initial_int']==i)].reset_index()
        #select_columns = doc_data.index[doc_data['FD_initial_cut_n']==fd_bins[2]].tolist()
        non_0_cols = doc_data[doc_data.dec_ratio5!=0].index.tolist()
        eq_0_cols = doc_data[doc_data.dec_ratio5==0].index.tolist()
        doc_t_data = doc_data[dec_ratio_columns].T
        #ax.plot(doc_t_data.loc[:,select_columns])
        ax.plot(doc_t_data.loc[:,non_0_cols], c='lightgrey')
        if len(eq_0_cols)>0:
            ax.plot(doc_t_data.loc[:,eq_0_cols], c= "darkgoldenrod")
        if jidx==0:
            ax.set_ylabel(i)
        if iidx==0:
            ax.set_title(j)
labels = list(str(int(x)) for x in labelist[::8])
plt.ylim(bottom = 0)
plt.xticks(np.arange(0, 9, step=1), labelist)
fig.supylabel("Normalized decay constant", fontsize = 16)
fig.supxlabel("% remaining C", fontsize = 16)
#%%
act_full['logdec0']=np.log10(act_full.dec0)
g = sns.FacetGrid(col='C_pool', row = 'Variance', hue = 'DOC_initial_int', data = act_full, palette="flare")
g.map(sns.scatterplot, 'decay_max', 'decay_stop')
g.set(ylim=(0, 100), xlim=(0,100), xlabel="%C for peak decomposition", ylabel="%C for C decomposition ending")
g.add_legend()
#%%
act_full['carbonxbiomass_species']=act_full.carbon_species*act_full.biomass_species
sns.scatterplot(x='decay_max', y='decay_stop',data = act_full[act_full.C_pool=='DOC'], style = 'Variance', size = 'carbonxbiomass_species', hue = 'DOC_initial_int', hue_order = [15000.,10000.,5000.,2000.,1000.])
plt.xlim(0,100)
plt.ylim(-0,100)
plt.xlabel("%C at peak decomposition")
plt.ylabel("%C persistent")
plt.legend(bbox_to_anchor=(1.05,1))
plt.title("Fully active communities: DOC")
#%%
sns.scatterplot(x='decay_max', y='decay_stop',data = act_full[act_full.C_pool=='reduced_C'], style = 'Variance', size = 'carbonxbiomass_species', hue = 'DOC_initial_int', hue_order = [15000.,10000.,5000.,2000.,1000.])
plt.xlim(0,100)
plt.ylim(0,100)
plt.xlabel("%C at peak decomposition")
plt.ylabel("%C persistent")
plt.legend(bbox_to_anchor=(1.05,1))
plt.title("Fully active communities: reduced C")
#%%
sns.scatterplot(x='decay_max', y='decay_stop',data = act_full[act_full.C_pool=='oxidized_C'], style = 'Variance', size = 'carbonxbiomass_species', hue = 'DOC_initial_int', hue_order = [15000.,10000.,5000.,2000.,1000.])
plt.xlim(0,100)
plt.ylim(0,100)
plt.xlabel("%C at peak decomposition")
plt.ylabel("%C persistent")
plt.legend(bbox_to_anchor=(1.05,1))
plt.title("Fully active communities: oxidized C")
#%%
###--------------------------------------------------------
### WHAT KIND OF COMMUNITIES ALLOW FOR PREMATURE STOPPING OF C decomposition?
###--------------------------------------------------------
early_stop = act_full[act_full.dec_ratio5==0]
es_doc = early_stop[early_stop.C_pool=='DOC']
esd_cpoor = es_doc[es_doc.DOC_initial_int==1000.]
esd_crich = es_doc[es_doc.DOC_initial_int==15000.]
#%%
print(esd_cpoor.Variance.unique())
print(esd_crich.Variance.unique())
print(esd_cpoor.FD_initial_cut_n.unique())
print(esd_crich.FD_initial_cut_n.unique())
#%%
###--------------------------------------------------------
### FACTORIAL ANALYSIS:
# Compare time series of decay constants given different:
# initial FD, C availability,Variance
###--------------------------------------------------------

var_0_01 = act_full[act_full.Variance==0.01][dec_ratio_columns].values
var_0_1 = act_full[act_full.Variance==0.1][dec_ratio_columns].values
var_0_5 = act_full[act_full.Variance==0.5][dec_ratio_columns].values
var_1_0 = act_full[act_full.Variance==1.][dec_ratio_columns].values
var_1_5 = act_full[act_full.Variance==1.5][dec_ratio_columns].values
#Test for assumptions
#Test for normality
for dataset in [var_0_01, var_0_1, var_0_5, var_1_0, var_1_5]:
    print(scipy.stats.shapiro(dataset))
#all datasets are not normal
#Test for variance
print(scipy.stats.bartlett(var_0_01, var_0_1, var_0_5, var_1_0, var_1_5))
#%%
#Result: NA
#Assumptions to run factorial ANOVA are not met
kruskal = scipy.stats.kruskal(var_0_01, var_0_1, var_0_5, var_1_0, var_1_5, nan_policy = 'omit')
print(kruskal)
#pvalue is less than 0.05, so the decay constant is evolves differently for community type
#%%

doc_1k = act_full[act_full.DOC_initial_int==1000][dec_ratio_columns].values
doc_2k = act_full[act_full.DOC_initial_int==2000][dec_ratio_columns].values
doc_5k = act_full[act_full.DOC_initial_int==5000][dec_ratio_columns].values
doc_10k = act_full[act_full.DOC_initial_int==10000][dec_ratio_columns].values
doc_15k = act_full[act_full.DOC_initial_int==15000][dec_ratio_columns].values
#Test for assumptions
#Test for normality
for dataset in [doc_1k, doc_2k,doc_5k,doc_10k, doc_15k]:
    print(scipy.stats.shapiro(dataset))
#all datasets are not normal
#Test for variance
print(scipy.stats.bartlett(doc_1k, doc_2k,doc_5k,doc_10k, doc_15k))
#Result: NA
#%%
#Assumptions to run factorial ANOVA are not met
kruskal = scipy.stats.kruskal(doc_1k, doc_2k,doc_5k,doc_10k, doc_15k, nan_policy = 'omit')
print(kruskal)
#pvalue is less than 0.05, so the decay constant evolution is different for carbon availability
#%%
fd_1 = act_full[act_full.FD_initial_cut_n==fd_bins[0]][dec_ratio_columns].values
fd_2 = act_full[act_full.FD_initial_cut_n==fd_bins[1]][dec_ratio_columns].values
fd_3 = act_full[act_full.FD_initial_cut_n==fd_bins[2]][dec_ratio_columns].values
#Test for assumptions
#Test for normality
for dataset in [fd_1, fd_2,fd_3]:
    print(scipy.stats.shapiro(dataset))
#fd_1 and fd_2 are normal (p>0.05) and fd_3 is normal (p>0.01)
#%%
#ANOVA
F,p = scipy.stats.f_oneway(fd_1,fd_2,fd_3)
# p<0.01 for all time points, so they are not signficantly different from each other
#%%
#Assumptions to run factorial ANOVA are not met
kruskal = scipy.stats.kruskal(fd_1, fd_2,fd_3, nan_policy = 'omit')
print(kruskal)
#pvalue is less than 0.05, so the decay constant evolution is different for carbon availability
#%%
###--------------------------------------------------------
### FACTORIAL ANALYSIS:
# Compare time series of decay constants given different:
# initial FD, C availability,Variance
###--------------------------------------------------------
import statsmodels.api as sm
import statsmodels.formula.api as smf
#%%
max_vd = smf.ols('decay_max ~ Variance + DOC_initial_int', act_full)
max_vd = max_vd.fit()
print(max_vd.summary())
#%%
#%%
max_vdcb = smf.ols('decay_max ~ Variance + DOC_initial_int + carbon_species + biomass_species', act_full)
max_vdcb = max_vdcb.fit()
print(max_vdcb.summary())
#%%
max_fd = smf.ols('decay_max ~ FD_initial + DOC_initial_int', act_full)
max_fd = max_fd.fit()
print(max_fd.summary())
#%%
max_fdcb = smf.ols('decay_max ~ FD_initial + DOC_initial_int + carbon_species + biomass_species', act_full)
max_fdcb = max_fdcb.fit()
print(max_fdcb.summary())
#%%
stop_vd = smf.ols('decay_stop ~ Variance + DOC_initial_int', act_full)
stop_vd = stop_vd.fit()
print(stop_vd.summary())
#%%
#%%
stop_vdcb = smf.ols('decay_stop ~ Variance + DOC_initial_int + carbon_species + biomass_species', act_full)
stop_vdcb = stop_vdcb.fit()
print(stop_vdcb.summary())
#%%
stop_fd = smf.ols('decay_stop ~ FD_initial + DOC_initial_int', act_full)
stop_fd = stop_fd.fit()
print(stop_fd.summary())
#%%
stop_fdcb = smf.ols('decay_stop ~ FD_initial + DOC_initial_int + carbon_species + biomass_species', act_full)
stop_fdcb = stop_fdcb.fit()
print(stop_fdcb.summary())
#%%
act_off = all_data[all_data.Sim_series!="b_1_all_"]
act_off['cab']=act_off.carbon_species*act_off.biomass_species*act_off.activity/100
#%%
sns.scatterplot(x='decay_max', y='decay_stop',data = act_off[act_off.C_pool=='DOC'], style = 'Variance', size = 'cab', hue = 'DOC_initial_int', hue_order = [15000.,10000.,5000.,2000.,1000.])
plt.xlim(0,100)
plt.ylim(0,100)
plt.xlabel("%C at peak decomposition")
plt.ylabel("%C persistent")
plt.legend(bbox_to_anchor=(1.05,1))
plt.title("Partially active communities: DOC")
#%%
sns.scatterplot(x='decay_max', y='decay_stop',data = act_off[act_off.C_pool=='reduced_C'], style = 'Variance', size = 'cab', hue = 'DOC_initial_int', hue_order = [15000.,10000.,5000.,2000.,1000.])
plt.xlim(0,100)
plt.ylim(0,100)
plt.xlabel("%C at peak decomposition")
plt.ylabel("%C persistent")
plt.legend(bbox_to_anchor=(1.05,1))
plt.title("Partially active communities: reduced C")
#%%
sns.scatterplot(x='decay_max', y='decay_stop',data = act_off[act_off.C_pool=='oxidized_C'], style = 'Variance', size = 'cab', hue = 'DOC_initial_int', hue_order = [15000.,10000.,5000.,2000.,1000.])
plt.xlim(0,100)
plt.ylim(0,100)
plt.xlabel("%C at peak decomposition")
plt.ylabel("%C persistent")
plt.legend(bbox_to_anchor=(1.05,1))
plt.title("Partially active communities: oxidized C")
#%%
#sns.displot(data = act_off[act_off.C_pool=='DOC'], x = 'decay_max', hue = 'DOC_initial_int', kind = 'kde')
fig, axe = plt.subplots(3,2,sharex=True, sharey=True, figsize = (9,7))
axes = axe.flatten()
sns.scatterplot(data=act_full[act_full.C_pool=='DOC'], x = 'decay_max', y = 'decay_stop', hue ='DOC_initial_int', style = act_full[act_full.C_pool=='DOC']['FD_initial_cut_n'].astype(str), palette = "flare", ax = axes[0])
sns.scatterplot(data=act_full[act_full.C_pool=='reduced_C'], x = 'decay_max', y = 'decay_stop', hue ='DOC_initial_int', style = act_full[act_full.C_pool=='reduced_C']['FD_initial_cut_n'].astype(str), palette = "flare", ax = axes[2])
sns.scatterplot(data=act_full[act_full.C_pool=='oxidized_C'], x = 'decay_max', y = 'decay_stop', hue ='DOC_initial_int', style = act_full[act_full.C_pool=='oxidized_C']['FD_initial_cut_n'].astype(str), palette = "flare", ax = axes[4])
sns.scatterplot(data=act_off[act_off.C_pool=='DOC'], x = 'decay_max', y = 'decay_stop', hue ='DOC_initial_int', style = act_off[act_off.C_pool=='DOC']['FD_initial_cut_n'].astype(str), palette = "flare", ax = axes[1])
sns.scatterplot(data=act_off[act_off.C_pool=='reduced_C'], x = 'decay_max', y = 'decay_stop', hue ='DOC_initial_int', style = act_off[act_off.C_pool=='reduced_C']['FD_initial_cut_n'].astype(str), palette = "flare", ax = axes[3])
sns.scatterplot(data=act_off[act_off.C_pool=='oxidized_C'], x = 'decay_max', y = 'decay_stop', hue ='DOC_initial_int', style = act_off[act_off.C_pool=='oxidized_C']['FD_initial_cut_n'].astype(str), palette = "flare", ax = axes[5])
for a in axes:
    a.legend().remove()
    a.set_ylim(0,100)
    a.set_xlim(0,100)
axes[0].set_ylabel("DOC")
axes[2].set_ylabel("Reduced C")
axes[4].set_ylabel("Oxidized C")
axes[0].set_title("Fully active\ncommunities")
axes[1].set_title("Partially active\ncommunities")
handles,labels=axes[-1].get_legend_handles_labels()
labels[0]="C availability"
labels[-4]="f (log scale)"
plt.figlegend(handles,labels,fontsize = 12, bbox_to_anchor=(0.9,0.5), loc = "lower left", borderpad=0.)
fig.supxlabel("%carbon for peak decomposition")
fig.supylabel("Persistent carbon (%)")
plt.savefig(os.path.join(project_dir, "full_partial_decay_max_stop.png"), dpi = 300, bbox_inches='tight')
#%%
sns.jointplot(data=act_off[act_off.C_pool=='reduced_C'], x = 'decay_max', y = 'decay_stop', hue = 'cab')
#%%
sns.jointplot(data=act_off[act_off.C_pool=='oxidized_C'], x = 'decay_max', y = 'decay_stop', hue = 'DOC_initial_int')
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