#%%
import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

import matplotlib.lines as mlines

import statsmodels.api as sm
import statsmodels.formula.api as smf

## LOAD RESULTS
project_dir = "C:/Users/swkh9804/Documents/Scripts/HoliSoils/Data"
figures_dir = os.path.join(project_dir, "Figures")
all_data = pd.read_csv(os.path.join(project_dir,"simulation_results_temporal_initial_conditions_decay_const_with_paras.csv"))
print(all_data.shape)
print(all_data.dtypes)
#%%
compl = all_data#.dropna()#subset = ['decay_const_initial'])
act_full = compl[compl.Sim_series=="b_1_all_"]
maxid = act_full.groupby(["Seed","C_pool","DOC_initial_int","biomass_species", "carbon_species", "Variance"])['Biomass_ratio'].transform('idxmax').values
act_full['Biomass_max_ratio'] = act_full.loc[maxid,'Biomass_ratio'].values
act_full["Biomass_max_ratio"].mask(act_full["Biomass_max_ratio"] <1, 1.0, inplace=True)
act_full['fd_ratio_at_max_Biomass'] = act_full.loc[maxid,'FD_ratio'].values
act_full.loc[act_full["Biomass_max_ratio"]==1.0, 'fd_ratio_at_max_Biomass']=1.0
extr_full = act_full[(act_full['DOC_initial_int']==2000.)|(act_full['DOC_initial_int']==10000.)]
toc = extr_full[(extr_full['C_pool']=='TOC')&(extr_full['%C']==100.)]
### Log transform the data:
toc_all_data = act_full[act_full.C_pool=='TOC']
no_inf = toc_all_data.drop(toc_all_data[toc_all_data.decay_const_initial==np.inf].index)
log_data= no_inf.apply(lambda x: np.log10(x) if np.issubdtype(x.dtype, np.number) else x)
#%%
### Model 1: Y = f(f0, C0)
def bio_ratio_f_fd_doc (df):
    df_sub  = df[df.C_pool=='TOC']
    model = smf.ols('Biomass_max_ratio ~ FD_initial + DOC_initial_int', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results

def fd_ratio_f_fd_doc (df):
    df_sub  = df[df.C_pool=='TOC']
    model = smf.ols('fd_ratio_at_max_Biomass ~ FD_initial + DOC_initial_int', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results

def dec_cont_f_fd_doc (df, cpool_test):
    df_sub  = df[df.C_pool==cpool_test]
    model = smf.ols('decay_const_initial ~ FD_initial + DOC_initial_int', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results

def dec_stop_f_fd_doc (df, cpool_test):
    df_sub  = df[df.C_pool==cpool_test]
    model = smf.ols('decay_stop ~ FD_initial + DOC_initial_int', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results

### Generate models:
model_1_1_m, model_1_1_r = bio_ratio_f_fd_doc(log_data)
model_1_2_m, model_1_2_r = fd_ratio_f_fd_doc(log_data)
model_1_3_m, model_1_3_r = dec_cont_f_fd_doc(log_data, 'TOC')
model_1_4_m, model_1_4_r = dec_stop_f_fd_doc(log_data, 'TOC')
coefs_model_1 = pd.DataFrame()
coefs_model_1['Features'] = ['f0', 'C0']
coefs_model_1['Biomass_ratio'] = model_1_1_r.params[1:].values
coefs_model_1['FD_ratio'] = model_1_2_r.params[1:].values
coefs_model_1['decay_constant'] = model_1_3_r.params[1:].values
coefs_model_1['Ct'] = model_1_4_r.params[1:].values
coefs_model_1.set_index('Features', inplace=True)
#%%
### Model 2: Y = f(Nc, Nb, Vb, C0)
def bio_ratio_f_cbvdoc (df):
    df_sub  = df[df.C_pool=='TOC']
    model = smf.ols('Biomass_max_ratio ~ carbon_species + biomass_species + Variance + DOC_initial_int', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results

def fd_ratio_f_cbvdoc (df):
    df_sub  = df[df.C_pool=='TOC']
    model = smf.ols('fd_ratio_at_max_Biomass ~ carbon_species + biomass_species + Variance + DOC_initial_int', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results

def dec_cont_f_cbvdoc (df, cpool_test):
    df_sub  = df[df.C_pool==cpool_test]
    model = smf.ols('decay_const_initial ~ carbon_species + biomass_species + Variance + DOC_initial_int', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results

def dec_stop_f_cbvdoc (df, cpool_test):
    df_sub  = df[df.C_pool==cpool_test]
    model = smf.ols('decay_stop ~ carbon_species + biomass_species + Variance + DOC_initial_int', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results

### Generate models:
model_2_1_m, model_2_1_r = bio_ratio_f_cbvdoc(log_data)
model_2_2_m, model_2_2_r = fd_ratio_f_cbvdoc(log_data)
model_2_3_m, model_2_3_r = dec_cont_f_cbvdoc(log_data, 'TOC')
model_2_4_m, model_2_4_r = dec_stop_f_cbvdoc(log_data, 'TOC')
coefs_model_2 = pd.DataFrame()
coefs_model_2['Features'] = ['Nc', 'Nb', 'Vb', 'C0']
coefs_model_2['Biomass_ratio'] = model_2_1_r.params[1:].values
coefs_model_2['FD_ratio'] = model_2_2_r.params[1:].values
coefs_model_2['decay_constant'] = model_2_3_r.params[1:].values
coefs_model_2['Ct'] = model_2_4_r.params[1:].values
coefs_model_2.set_index('Features', inplace=True)

#%%
### Figure 3: Scatter plot sof ecosystem indicators against function diversity
datatoplot = toc[(toc['Variance']==0.01)|(toc['Variance']==1.0)]
x_model_1 = np.linspace(1E-10, 1E-3, num = 80)#np.array(datatoplot.FD_initial.unique().sort())
y_model_1_lowdoc = 10**(model_1_1_r.params[0])*(x_model_1**model_1_1_r.params[1])*(2000.**model_1_1_r.params[2])
y_model_1_highdoc = 10**(model_1_1_r.params[0])*(x_model_1**model_1_1_r.params[1])*(10000.**model_1_1_r.params[2])
y_model_2_lowdoc = 10**(model_1_2_r.params[0])*(x_model_1**model_1_2_r.params[1])*(2000.**model_1_2_r.params[2])
y_model_2_highdoc = 10**(model_1_2_r.params[0])*(x_model_1**model_1_2_r.params[1])*(10000.**model_1_2_r.params[2])
y_model_3_lowdoc = 10**(model_1_3_r.params[0])*(x_model_1**model_1_3_r.params[1])*(2000.**model_1_3_r.params[2])
y_model_3_highdoc = 10**(model_1_3_r.params[0])*(x_model_1**model_1_3_r.params[1])*(10000.**model_1_3_r.params[2])
y_model_4_lowdoc = 10**(model_1_4_r.params[0])*(x_model_1**model_1_4_r.params[1])*(2000.**model_1_4_r.params[2])
y_model_4_highdoc = 10**(model_1_4_r.params[0])*(x_model_1**model_1_4_r.params[1])*(10000.**model_1_4_r.params[2])
#%%
fig, axes = plt.subplots(4,1, figsize = (4,12), sharex = True,sharey=False)
axflat=axes.flatten()
sns.scatterplot(x = "FD_initial", y = "Biomass_max_ratio", hue = 'DOC_initial_int', style = "Variance", markers = ['o', 'v'], palette= ['goldenrod', 'olive'], facecolor = 'None', data = datatoplot, ax = axflat[0])
axflat[0].plot(x_model_1,y_model_1_lowdoc, c= "maroon", label = "low $C_{0}$")
axflat[0].plot(x_model_1,y_model_1_highdoc, c= "darkgreen", label = "high $C_{0}$")
axflat[0].axhline(y=1.0, c = 'grey', linestyle = 'dashed')
axflat[0].annotate('Dying\nsystems', xytext=(0.1,0.75),  xycoords='data',
            xy=(1E-10,1.0), textcoords='axes fraction',color='dimgrey',
            arrowprops=dict(facecolor='grey', width = 3, edgecolor = 'grey', shrink=0.05),
            horizontalalignment='left', verticalalignment='top',fontsize = 12)
axflat[0].set_ylabel("$B_{max}/B_{0}$", fontsize = 16, labelpad = 55)
axflat[0].set_title('A', fontsize = 16, loc='left')
axflat[0].text(1E-10, 2.4, '$R^{2}$ = 0.60', fontsize = 14)

axflat[1].plot(x_model_1, y_model_2_lowdoc, c= "maroon", label = "low $C_{0}$")
axflat[1].plot(x_model_1, y_model_2_highdoc, c= "darkgreen", label = "high $C_{0}$")
sns.scatterplot(x = "FD_initial", y = "fd_ratio_at_max_Biomass", hue = 'DOC_initial_int', style = "Variance", markers = ['o', 'v'], palette= ['goldenrod', 'olive'], data = datatoplot, ax = axflat[1])
axflat[1].set_title('B', fontsize = 16, loc='left')
axflat[1].axhline(y=1., c = 'grey', linestyle = 'dashed')
axflat[1].annotate('No change in\nfunctional\ndiversity', xytext=(0.1,0.75),  xycoords='data',
            xy=(1E-10,1.0), textcoords='axes fraction',color='dimgrey',
            arrowprops=dict(facecolor='grey', width = 3, edgecolor = 'grey', shrink=0.05),
            horizontalalignment='left', verticalalignment='top',fontsize = 12)
axflat[1].set_ylabel("$f_{B_{max}}/f_{0}$", fontsize = 16, labelpad = 45)
axflat[1].set_title('B', fontsize = 16, loc='left')
axflat[1].text(1E-10, 7.9, '$R^{2}$ = 0.50', fontsize = 14)

sns.scatterplot(x = "FD_initial", y = "decay_const_initial", hue = 'DOC_initial_int', style = "Variance", markers = ['o', 'v'], palette= ['goldenrod', 'olive'], data = datatoplot, ax = axflat[2])
axflat[2].plot(x_model_1, y_model_3_lowdoc, c= "maroon", label = "low $C_{0}$")
axflat[2].plot(x_model_1, y_model_3_highdoc, c= "darkgreen", label = "high $C_{0}$")
axflat[2].set_ylabel("$k_{10}(d^{-1})$", fontsize = 16, labelpad = 25)
axflat[2].set_title('C', fontsize = 16, loc='left')
axflat[2].set_xlabel("$f_{0}$", fontsize = 16)
axflat[2].text(1E-10, 0.0075, '$R^{2}$ = 0.75', fontsize = 14)

axflat[3].plot(x_model_1, y_model_4_lowdoc, c= "maroon", label = "low $C_{0}$")
axflat[3].plot(x_model_1, y_model_4_highdoc, c= "darkgreen", label = "high $C_{0}$")
sns.scatterplot(x = "FD_initial", y = "decay_stop", hue = 'DOC_initial_int', style = "Variance", markers = ['o', 'v'], palette= ['goldenrod', 'olive'], data = datatoplot, ax = axflat[3])
axflat[3].set_ylabel("$C_{T}/C_{0}$ (%)", fontsize = 16, labelpad = 55)
axflat[3].axhline(y=10., c = 'grey', linestyle = 'dashed')
axflat[3].annotate('Complete\ndecomposition', xytext=(0.9,0.82),  xycoords='data',
            xy=(1E-3,10), textcoords='axes fraction', color='dimgrey',
            arrowprops=dict(facecolor='grey', width = 3, edgecolor = 'grey', shrink=0.05),
            horizontalalignment='right', verticalalignment='top',fontsize = 12)
axflat[3].set_title('D', fontsize = 16, loc='left')
axflat[3].text(1E-10, 50, '$R^{2}$ = 0.39', fontsize = 14)
axflat[3].set_xlabel("$f_{0}$", fontsize = 16)
axflat[3].set_xscale("log")

for a in axflat:
    a.tick_params(axis='both', which='major', labelsize=14)
    a.legend().set_visible(False)

leg_handles, leg_labels = axflat[3].get_legend_handles_labels()
leg_labels[0] = "$C_{0}$ poor"
leg_labels[1] = "$C_{0}$ rich"
leg_labels[2]="Carbon availability"
leg_labels[3]="poor"
leg_labels[4]="rich"
l = axflat[1].legend(leg_handles[:2], leg_labels[:2], bbox_to_anchor=(0.1,-3.1), loc = 'lower left', ncol=3,title_fontsize = 16, fontsize = 16, handleheight = 1.2, columnspacing = 1., handletextpad = 0.3, frameon = False)
plt.text(3E-13, -25.5, 'Regression', fontsize = 16)
axflat[2].legend(leg_handles[2:5], leg_labels[2:5], bbox_to_anchor=(-0.5,-2.1), loc = 'lower left', ncol=3,title_fontsize = 16, fontsize = 16, handleheight = 1.2, columnspacing = 1., handletextpad = 0.3, frameon = False)
axflat[3].legend(leg_handles[5:], leg_labels[5:], bbox_to_anchor=(-0.5,-1.1), loc = 'lower left', ncol=3,title_fontsize = 16, fontsize = 16, handleheight = 1.2, columnspacing = 1., handletextpad = 0.3, frameon = False)
plt.tight_layout()
plt.savefig(os.path.join(figures_dir, "Fig_3_ecosystem_func_div.png"), bbox_inches='tight',pad_inches = 0, dpi = 300)
#%%
### Figure 4: Heatmaps of coefficients and R2 of Model 1.
sns.heatmap(coefs_model_1.transpose(), annot = True, fmt=".2f", center = 0, cmap = 'vlag', xticklabels=["$f_{0}$","$C_{0}$"], yticklabels=["$B_{max}/B_{0}$","$f_{B_{max}}/f_{0}$","$k_{10}$","$C_{T}/C_{0}$"])
plt.savefig(os.path.join(figures_dir, "Fig_4_model_1_coefs.png"), bbox_inches='tight',pad_inches = 0, dpi = 300)
#%%
### Figure 5: Heatmaps of coefficients and R2 of Model 2.
sns.heatmap(coefs_model_2.transpose(), annot = True, fmt=".2f", center = 0, cmap = 'vlag', xticklabels=["Nc","Nb", "$V_{b}$","$C_{0}$"], yticklabels=["$B_{max}/B_{0}$","$f_{B_{max}}/f_{0}$","$k_{10}$","$C_{T}/C_{0}$"])
plt.savefig(os.path.join(figures_dir, "Fig_5_model_2_coefs.png"), bbox_inches='tight',pad_inches = 0, dpi = 300)
#%%

###--------SUPPLEMENTARY FIGURES -------------###
#%%
fig, axes = plt.subplots(1,2, sharey = True, sharex = False, figsize = (10,4))
axe = axes.flatten()
#plt.figure(figsize=(8,4))
h = sns.scatterplot(x = "Variance", y = "FD_initial", hue = "carbon_species", size = "biomass_species", data = act_full, ax = axe[0])
h.set_xlabel("$V_{b}$", fontsize = 14)
h.set_ylabel("$f_{0}$", fontsize = 14)
leg_handles, leg_labels = h.get_legend_handles_labels()
leg_labels[0]="Nc"
leg_labels[5]="Nb"
axe[0].legend(leg_handles, leg_labels,bbox_to_anchor=(1.,1.05), fontsize = 14, frameon = False)
axe[0].set_title("A", fontsize = 14, loc = "left")

g = sns.scatterplot(x = "S_initial", y = "FD_initial", hue = 'carbon_species', style = "Variance", data = act_full, ax = axe[1])
g.set_xlabel("$H_{0}$", fontsize = 14)
leg_handles, leg_labels = g.get_legend_handles_labels()
leg_labels[0]="Nc"
leg_labels[5]="$V_{b}$"
axe[1].legend(leg_handles, leg_labels,bbox_to_anchor=(1.,1.05), fontsize = 14, frameon = False)
axe[1].set_title("B", fontsize = 14, loc = "left")

for a in axe:
    a.tick_params(axis='both', labelsize=12)
    a.set_yscale("log")
plt.tight_layout()
plt.savefig(os.path.join(figures_dir, "Fig_S1_imposed_func_div.png"), bbox_inches = 'tight', dpi = 300)
#%%
toc_all_data['FD_initial_cut_n'] = pd.cut(np.log10(toc_all_data['FD_initial']), bins=[-13,-7,-5,-1])
x = toc_all_data['FD_initial_cut_n'].astype(str)

#%%
### Line plot of slowing down decay constant with decreasing carbon
#plotdata_grp = toc_all_data.groupby(["Seed", "C_pool", "Sim_series", "Variance","DOC_initial_int", "carbon_species", "biomass_species"])
x_lin = np.linspace(100,10,10)
row_plots = [2000.0, 10000.]#[0.01, 0.1, 0.5, 1.0, 1.5]#x.unique().tolist()
col_plots = ['TOC','DOC', 'reduced_C', 'oxidized_C']
data = extr_full#[extr_full.FD_initial_cut_n==fd_bins[1]]
fig, axes = plt.subplots(len(row_plots),len(col_plots),sharex=True, sharey='row', figsize = (8,4))
ax = axes.flatten()
for iidx, i in enumerate(row_plots):
    for jidx, j in enumerate(col_plots):
        axindx = iidx*len(col_plots) + jidx
        subset = data[(data['C_pool']==j)&(data['DOC_initial_int']==i)]#['FD_initial_cut_n'].astype(str)==i)]
        actual = subset.groupby(["Seed", "Variance", "carbon_species", "biomass_species"])
        for (key) in actual.groups.keys():
            y = actual.get_group(key)['decay_ratio'].values
            if y[-4]==-1.:
                pass
            else:    
                ax[axindx].plot(x_lin, y, c = 'grey', alpha = 0.2)
        for (key) in actual.groups.keys():
            y = actual.get_group(key)['decay_ratio'].values
            if y[-4]==-1.:
                ax[axindx].plot(x_lin, y, c = 'darkgoldenrod')
            else:    
                pass
        if iidx==1:
            ax[axindx].set_xlabel("%$C_{0}$")
        if jidx==0:
            ax[axindx].set_ylabel(i, fontsize = 14)
        else:
            ax[axindx].set_ylabel("")
fig.supylabel("$C_{0}$", fontsize = 14)
ax[0].set_title("TOC", fontsize = 14)
ax[1].set_title("DOC", fontsize = 14)
ax[2].set_title("reduced C", fontsize = 14)
ax[3].set_title("oxidized C", fontsize = 14)
plt.ylim(bottom=-1.)
plt.xlim((100.,10.))
for a in ax[:]:
    a.axhline(y=1.0, c='maroon', linestyle='dashed')
plt.savefig(os.path.join(figures_dir, "Fig_S3_decay_const_temporal.png"), bbox_inches = 'tight', dpi = 300)
#%%
#### Microbial diversity, chemical diveristy in scenarios of incomplete carbon decomposition
group = extr_full[extr_full.C_pool=="TOC"].groupby(["Seed","Variance", "DOC_initial_int","carbon_species", "biomass_species"])
key_dic = []
for (key) in group.groups.keys():
    y = group.get_group(key)['decay_ratio'].values
    if y[-4]==-1.:
        key_dic.append(key)
#%%
vmax_dic = []
kmean_dic = []
nosc_dic = []
for k in key_dic:
    vmax_dic.append(group.get_group(k)['vmax_mean'])
    kmean_dic.append(group.get_group(k)['k_mean'])
    nosc_dic.append(group.get_group(k)['NOSC_initial'])
#%%
### PLOT TIME SERIES PLOTS OF THESE 7 SCENARIOS
### LOOK AT CONCENTRATIONS, and Carbon pool
import h5py

## LOAD RESULTS
data_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
filestring = "competition_adaptation"
ip = 0

sc1 = key_dic[4]
Seed_select = sc1[0]
var_select = sc1[1]
doc_low = sc1[2]
carbon_n = sc1[3]
bio_low= sc1[4]

simulation_subfolder = filestring + '_carbon_' + str(carbon_n) + '_'+str(Seed_select) + '_ip_0'
similar_dir = os.path.join(data_dir, 'gen_spec_lognorm_0_'+str(var_select)[-1], "simulations", simulation_subfolder)
#similar_dir = os.path.join(data_dir, 'gen_spec_lognorm', "simulations", simulation_subfolder)
#similar_dir = os.path.join(data_dir, 'gen_spec_lognorm_1_5x', "simulations", simulation_subfolder)
similar_hr = h5py.File(os.path.join(similar_dir,"simulations.h5"), mode = 'r')

bl_docl = "b_1_all_/bio_n_" + str(bio_low) + "/dom_initial_" + str(doc_low) + "/seed_" + str(Seed_select)
sim_bl_docl = similar_hr[bl_docl+'/solution']
paras = similar_hr[bl_docl+'/parameters']
nosc = paras['oxidation_state'][:]
vprod = np.sum(np.asarray(paras['max_rate'])*np.asarray(paras['carbon_uptake'])*np.asarray(paras['exo_enzyme_rate']), axis=0)
x = sim_bl_docl
C = np.asarray(x['dom'])
B = np.asarray(x['biomass'])
C_sum = (1- (np.sum(np.asarray(x['dom']), axis = 1)/doc_low))*100
B_sum = np.sum(np.asarray(x['biomass']), axis = 1)

C_pc = (1- (np.asarray(x['dom'])/np.asarray(x['dom'])[0,:]))*100

nosc_argidx = np.argsort(nosc)
vm_argidx = np.argsort(vprod)

time_span = np.linspace(0,36505, int(36500/5))
xticks_plot = time_span[::730].astype(int)
xticks_label = np.linspace(0,100,100)[::10].astype(int)
cmap_bio = plt.cm.get_cmap('plasma_r')
#cmap_bio = plt.cm.viridis_r#((nosc-nosc.min())/(nosc.max()-nosc.min()))
#color_idx = np.arange(0, 100, step = 10, dtype = int)
fig,axes = plt.subplots(2,1, figsize = (5,6), sharex = True)
ax = axes.flatten()

for cn in list(range(carbon_n)):
    print(cn, nosc_argidx[cn])
    ax[0].plot(time_span, C_pc[:, nosc_argidx[cn]], c = cmap_bio(10+(50*cn)), label = nosc[nosc_argidx[cn]])
ax[0].legend(bbox_to_anchor=(1,1.))
ax[0].set_ylabel ("%$C_{0}$")
for bn in list(range(bio_low)):
    ax[1].plot(time_span, B[:,vm_argidx[bn]],c = cmap_bio(10+(50*bn)), label = vprod[vm_argidx[bn]])
ax[1].legend(bbox_to_anchor = (1.65,1.))
ax[1].set_ylabel ("B [N $L_{-3}$]")
ax[1].set_xlabel("Time (days)")
plt.suptitle(sc1)
#%%
print(np.shape(np.asarray(paras['max_rate'])))
#%%
print(paras['max_rate'][:])#[[2,3,4],:])#*np.asarray(paras['exo_enzyme_rate'])*np.asarray(paras['carbon_uptake']
#%%
for x,doc_lev,v in zip([sim_bl_docl],[doc_low], [0,0,1,1]):
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
#%%
sns.scatterplot(x = "Biomass_max_ratio", y = "fd_ratio_at_max_Biomass", hue = 'DOC_initial_int', style = "Variance", size = "carbon_species",data = toc[toc.biomass_species==4])
plt.axhline(y=1, c = 'red', label = "Same community")
plt.axvline(x=1, c = 'black', label = 'Dying community')
plt.xlim((0.5,5))
plt.ylim(top = 12)
plt.yscale("log")
plt.xlabel("Ratio: Biomass")
plt.ylabel("Ratio: Functional diversity")
plt.legend(bbox_to_anchor=(1, 1))
#%%
