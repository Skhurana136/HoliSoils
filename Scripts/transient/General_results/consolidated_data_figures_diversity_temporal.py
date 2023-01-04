#%%
import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

## LOAD RESULTS
project_dir = "C:/Users/swkh9804/Documents/Scripts/HoliSoils/Data"#os.path.join("D:/", "Projects", "HoliSoils","data","transient")
all_data = pd.read_csv(os.path.join(project_dir,"simulation_results_temporal_initial_conditions_decay_const.csv"))
print(all_data.shape)
print(all_data.dtypes)
seeds= all_data.Seed.unique().tolist()
bios = all_data.biomass_species.unique().tolist()
carbs = all_data.carbon_species.unique().tolist()
carb_types = all_data.C_pool.unique().tolist()
docs = all_data.DOC_initial_int.unique().tolist()
sims = all_data.Sim_series.unique().tolist()
variances=all_data.Variance.unique().tolist()
#all_data.replace([np.inf, -np.inf, np.nan], 0, inplace=True)
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
### SUBSET DATA TO WORK WITH SMALLER DATASETS TO MAKE SENS
###--------------------------------------------------------
act_full = all_data[all_data.Sim_series=="b_1_all_"]
extr_full = act_full[(act_full['DOC_initial_int']==2000.)|(act_full['DOC_initial_int']==10000.)]
#%%
###--------------------------------------------------------
### COMPARISON OF DECAY CONSTANT OF DIFFERENT CARBON POOLS
###--------------------------------------------------------
row_plots = [100.]
col_plots = ['TOC', 'DOC', 'reduced_C', 'oxidized_C']
data = extr_full#extr_full[extr_full.FD_initial_cut_n==fd_bins[1]]
fig, axes = plt.subplots(len(row_plots),len(col_plots),sharex=True, sharey=True, figsize = (10,4))
ax = axes.flatten()
for i in list(range(len(col_plots))):
    subset1 = data[data.C_pool==col_plots[i]].reset_index()
    for j in list(range(len(row_plots))):
        axindx = j*len(col_plots) + i
        subset = subset1[subset1['%C']==row_plots[j]]
        g=sns.scatterplot(x=subset['FD_initial'].astype(str),y=subset['decay_const_initial'], hue = subset['DOC_initial_int'], ax=ax[axindx], edgecolor='None', alpha = 0.7)
        g.legend().remove()
        ax[axindx].set_xlabel("")
ax[0].set_ylabel("Decay constant", fontsize = 14)
fig.supxlabel("Functional diversity (Variance)", fontsize = 14)
ax[0].set_title("DOC", fontsize = 14)
ax[1].set_title("reduced C", fontsize = 14)
ax[2].set_title("oxidized C", fontsize = 14)
plt.xscale("log")
for a in ax[:]:
    a.axhline(y=0.0, c='maroon', linestyle='dashed')
handles,labels=ax[axindx].get_legend_handles_labels()
plt.figlegend(handles,labels,title = 'C availability', fontsize = 12, title_fontsize = 12, bbox_to_anchor=(0.85,-0.1), ncol=5, loc = "lower right", borderpad=0.)
plt.show()
#%%
###--------------------------------------------------------
### COMPARISON OF DECAY CONSTANT OF DIFFERENT CARBON POOLS TEMPORALLY
###--------------------------------------------------------
row_plots = [90.,70.,50.,30.]
col_plots = ['TOC','DOC', 'reduced_C', 'oxidized_C']
data = extr_full#extr_full[extr_full.FD_initial_cut_n==fd_bins[1]]
fig, axes = plt.subplots(len(row_plots),len(col_plots),sharex=True, sharey=True, figsize = (10,10))
ax = axes.flatten()
for iidx, i in enumerate(row_plots):
    for jidx, j in enumerate(col_plots):
        axindx = jidx*len(col_plots) + iidx
        subset = data[(data['C_pool']==j)&(data['%C']==i)]
        g=sns.scatterplot(x=subset['FD_initial'].astype(str),y=subset['decay_ratio'], hue = subset['DOC_initial_int'], ax=ax[axindx], edgecolor='None', alpha = 0.7)
        g.legend().remove()
        ax[axindx].set_xlabel("")
        if i==0:
            ax[axindx].set_ylabel(str(row_plots[j])+"%", fontsize = 14)
        else:
            ax[axindx].set_ylabel("")
fig.supxlabel("Functional diversity (Variance)", fontsize = 14)
ax[0].set_title("DOC", fontsize = 14)
ax[1].set_title("reduced C", fontsize = 14)
ax[2].set_title("oxidized C", fontsize = 14)
plt.xscale("log")
plt.yscale("log")
for a in ax[:]:
    a.axhline(y=1.0, c='maroon', linestyle='dashed')
handles,labels=ax[axindx].get_legend_handles_labels()
plt.figlegend(handles,labels,title = 'C availability', fontsize = 12, title_fontsize = 12, bbox_to_anchor=(0.85,-0.1), ncol=5, loc = "lower right", borderpad=0.)
plt.show()
#%%
g= sns.FacetGrid (data = extr_full[extr_full['%C']<100.], col = 'C_pool', row = '%C', hue = 'DOC_initial_int', palette = sns.color_palette(['salmon','blue']), margin_titles=True,
row_order=[90.,80.,70.,60.,50.,40.,30.,20.,10.])
g.map(sns.scatterplot, 'decay_const_initial', 'decay_ratio')
g.set(xscale="log", ylim=(-5,25))
g.map(plt.axhline, y=-1, ls='--', c='black', label = 'End of decomposition')
g.map(plt.axhline, y=1, ls='--', c='maroon', label = 'Same speed')
g.add_legend()
#%%
g= sns.FacetGrid (data = extr_full[extr_full['%C']<100.], col = 'C_pool', row = '%C', hue = 'Variance', palette = 'magma_r', margin_titles=True,
row_order=[90.,80.,70.,60.,50.,40.,30.,20.,10.])
g.map(sns.scatterplot, 'decay_const_initial', 'decay_ratio')
g.set(xscale="log", ylim=(-5,25))
g.map(plt.axhline, y=-1, ls='--', c='black', label = 'End of decomposition')
g.map(plt.axhline, y=1, ls='--', c='maroon', label = 'Same speed')
g.add_legend()
#%%
###--------------------------------------------------------
### EXPLORE WHEN C DECOMPOSITION PEAKS AND WHEN IT STOPS
###--------------------------------------------------------
extr_full['size_var']=np.log10(extr_full.decay_const_initial)
all_data['color_var'] = all_data.active_H_c_connections
plt.figure()
g = sns.FacetGrid(data = all_data[all_data.C_pool=='TOC'], row = 'DOC_initial_int',col = 'Variance', hue = 'color_var', palette = 'magma_r', margin_titles=True)
g.map(sns.scatterplot, 'decay_max', 'decay_stop')
g.fig.subplots_adjust(top=0.9)
g.fig.suptitle('TOC')
plt.figure()
g = sns.FacetGrid(data = all_data[all_data.C_pool=='DOC'], row = 'DOC_initial_int',col = 'Variance', hue = 'color_var', palette = 'magma_r', margin_titles=True)
g.map(sns.scatterplot, 'decay_max', 'decay_stop')
g.fig.subplots_adjust(top=0.9)
g.fig.suptitle('DOC')
plt.figure()
g = sns.FacetGrid(data = all_data[all_data.C_pool=='reduced_C'], row = 'DOC_initial_int',col = 'Variance', hue = 'color_var', palette = 'magma_r', margin_titles=True)
g.map(sns.scatterplot, 'decay_max', 'decay_stop')
g.fig.subplots_adjust(top=0.9)
g.fig.suptitle('reduced_C')
plt.figure()
g = sns.FacetGrid(data = all_data[all_data.C_pool=='oxidized_C'], row = 'DOC_initial_int',col = 'Variance', hue = 'color_var', palette = 'magma_r', margin_titles=True)
g.map(sns.scatterplot, 'decay_max', 'decay_stop')
g.fig.subplots_adjust(top=0.9)
g.fig.suptitle('oxidized_C')
plt.show()
#%%
#act_full['logdec0']=np.log10(act_full.dec0)
g = sns.FacetGrid(col='C_pool', row = 'Variance', hue = 'DOC_initial_int', data = extr_full, palette="flare")
g.map(sns.scatterplot, 'decay_max', 'decay_stop')
g.set(ylim=(0, 110), xlim=(0,110), xlabel="%C for peak decomposition", ylabel="%C for C decomposition ending")
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
early_stop = act_full[(act_full['decay_ratio']==-1.)&(act_full['%C']==50.)]
es_doc = early_stop[early_stop.C_pool=='DOC']
esd_cpoor = es_doc[es_doc.DOC_initial_int<=2000.]
esd_crich = es_doc[es_doc.DOC_initial_int>2000.]
#%%
print(esd_cpoor.Variance.unique())
print(esd_crich.Variance.unique())
print(esd_cpoor.FD_initial_cut_n.unique())
print(esd_crich.FD_initial_cut_n.unique())
# Carbon poor conditions provide for storage as opposed to carbonrich conditions
#%%
###--------------------------------------------------------
### FACTORIAL ANALYSIS:
# Compare time series of decay constants given different:
# initial FD, C availability,Variance
###--------------------------------------------------------

var_0_01 = act_full[act_full.Variance==0.01]['decay_ratio']
var_0_1 = act_full[act_full.Variance==0.1]['decay_ratio']
var_0_5 = act_full[act_full.Variance==0.5]['decay_ratio']
var_1_0 = act_full[act_full.Variance==1.]['decay_ratio']
var_1_5 = act_full[act_full.Variance==1.5]['decay_ratio']
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

doc_1k = act_full[act_full.DOC_initial_int==1000]['decay_ratio']
doc_2k = act_full[act_full.DOC_initial_int==2000]['decay_ratio']
doc_5k = act_full[act_full.DOC_initial_int==5000]['decay_ratio']
doc_10k = act_full[act_full.DOC_initial_int==10000]['decay_ratio']
doc_15k = act_full[act_full.DOC_initial_int==15000]['decay_ratio']
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
fd_1 = act_full[act_full.FD_initial_cut_n==fd_bins[0]]['decay_ratio']
fd_2 = act_full[act_full.FD_initial_cut_n==fd_bins[1]]['decay_ratio']
fd_3 = act_full[act_full.FD_initial_cut_n==fd_bins[2]]['decay_ratio']
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
extr_off = act_off[(act_off['DOC_initial_int']==2000.)|(act_off['DOC_initial_int']==10000.)]
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
plt_full = extr_full
plt_off = extr_off
fig, axe = plt.subplots(4,2,sharex=True, sharey=True, figsize = (8,8))
axes = axe.flatten()
sns.scatterplot(data=plt_full[plt_full.C_pool=='TOC'], x = 'decay_max', y = 'decay_stop', hue ='DOC_initial_int', style = plt_full[plt_full.C_pool=='TOC']['FD_initial_cut_n'].astype(str), palette = "flare", markers = ['o', 'X', '^'], ax = axes[0])
sns.scatterplot(data=plt_full[plt_full.C_pool=='DOC'], x = 'decay_max', y = 'decay_stop', hue ='DOC_initial_int', style = plt_full[plt_full.C_pool=='DOC']['FD_initial_cut_n'].astype(str), palette = "flare",  markers = ['o', 'X', '^'], ax = axes[2])
sns.scatterplot(data=plt_full[plt_full.C_pool=='reduced_C'], x = 'decay_max', y = 'decay_stop', hue ='DOC_initial_int', style = plt_full[plt_full.C_pool=='reduced_C']['FD_initial_cut_n'].astype(str), palette = "flare",  markers = ['o', 'X', '^'], ax = axes[4])
sns.scatterplot(data=plt_full[plt_full.C_pool=='oxidized_C'], x = 'decay_max', y = 'decay_stop', hue ='DOC_initial_int', style = plt_full[plt_full.C_pool=='oxidized_C']['FD_initial_cut_n'].astype(str), palette = "flare",  markers = ['o', 'X', '^'], ax = axes[6])
sns.scatterplot(data=plt_off[plt_off.C_pool=='TOC'], x = 'decay_max', y = 'decay_stop', hue ='DOC_initial_int', style = plt_off[plt_off.C_pool=='TOC']['FD_initial_cut_n'].astype(str), palette = "flare",  markers = ['o', 'X', '^'], ax = axes[1])
sns.scatterplot(data=plt_off[plt_off.C_pool=='DOC'], x = 'decay_max', y = 'decay_stop', hue ='DOC_initial_int', style = plt_off[plt_off.C_pool=='DOC']['FD_initial_cut_n'].astype(str), palette = "flare",  markers = ['o', 'X', '^'], ax = axes[3])
sns.scatterplot(data=plt_off[plt_off.C_pool=='reduced_C'], x = 'decay_max', y = 'decay_stop', hue ='DOC_initial_int', style = plt_off[plt_off.C_pool=='reduced_C']['FD_initial_cut_n'].astype(str), palette = "flare",  markers = ['o', 'X', '^'], ax = axes[5])
sns.scatterplot(data=plt_off[plt_off.C_pool=='oxidized_C'], x = 'decay_max', y = 'decay_stop', hue ='DOC_initial_int', style = plt_off[plt_off.C_pool=='oxidized_C']['FD_initial_cut_n'].astype(str), palette = "flare",  markers = ['o', 'X', '^'], ax = axes[7])
for a in axes:
    a.legend().remove()
    a.set_ylim(0,110)
    a.set_xlim(0,110)
axes[0].set_ylabel("TOC")
axes[2].set_ylabel("DOC")
axes[4].set_ylabel("Reduced C")
axes[6].set_ylabel("Oxidized C")
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
g = sns.FacetGrid(data=extr_full[extr_full.C_pool=='TOC'], col = 'FD_initial_cut_n', row = 'DOC_initial_int', hue = 'active_H_c_connections', margin_titles=True, palette= 'flare')
g.map(sns.scatterplot, 'decay_max', 'decay_stop')
#g.add_legend()
#%%
extr_full['carbonxbiomass'] = extr_full.carbon_species*extr_full.biomass_species
#%%
cpooltoplot = 'reduced_C'
g = sns.FacetGrid(data=extr_full[extr_full.C_pool==cpooltoplot], col = 'FD_initial_cut_n', row = 'DOC_initial_int', hue = 'carbonxbiomass', margin_titles=True, palette= 'flare')
g.map(sns.scatterplot, 'decay_max', 'decay_stop')
g.set(xlim=(0,110), ylim=(0,110),xlabel="%C for peak decomposition", ylabel="%C persists")
g.set_titles(col_template="logf range= {col_name}", row_template="Initial C= {row_name}")#-13<logf<-7","-7<logf<-5","-5<logf<-1")#,"Initial C = 2000.","Initial C = 10000.")
g.fig.subplots_adjust(top=0.90)
g.fig.suptitle(cpooltoplot)
#%%
sns.jointplot(data=act_off[act_off.C_pool==cpooltoplot], x = 'decay_max', y = 'decay_stop', hue = 'cab')