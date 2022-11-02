#%%
import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

import matplotlib.lines as mlines

from scipy import stats
from scipy.optimize import curve_fit

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
sim_suffixes = ["_0_01", "_0_1", "_0_5", "", "_1_5x"]
sim_suffixes_var = [0.01, 0.1, 0.5, 1, 1.5]

files = []
for p,a in zip(sim_suffixes, sim_suffixes_var):
    filename = os.path.join(project_dir, "gen_spec_lognorm" + p, "results", "competition_adaptation_carbon__loss_0.9_cue_combined_dataset.pkl")
    data = pd.read_pickle(filename)
    var_arr = np.zeros((data.shape[0],))+a
    var_ser = pd.Series(var_arr, copy=False,name = "Variance")
    data_var = data.join(var_ser)
    files.append(data_var)
fd_data = pd.concat(files)
print(fd_data.columns)
fd_data['DOC_initial_int'] = round(fd_data.DOC_initial, -3)
#%%
files = []
for p,a in zip(sim_suffixes, sim_suffixes_var):
    filename = os.path.join(project_dir, "gen_spec_lognorm" + p, "results", "gen_spec_lognorm" + p + "_parameters.csv")
    data = pd.read_csv(filename)
    var_arr = np.zeros((data.shape[0],))+a
    var_ser = pd.Series(var_arr, copy=False,name = "Variance")
    data_var = data.join(var_ser)
    files.append(data_var)
para_data = pd.concat(files)
print(para_data.columns)
cols_to_merge = para_data.columns.to_list()[2:]
#%%
fd_para_data = pd.merge(fd_data, para_data[cols_to_merge], on = ["Seed", "Variance", "biomass_species", "carbon_species", "Sim_series"])
#%%
### CHARACTERISTIC REACTION TIME SCALE
files = []
for p,a in zip(sim_suffixes, sim_suffixes_var):
    filename = os.path.join(project_dir, "gen_spec_lognorm" + p, "results", "competition_adaptation_carbon__loss_0.9_combined_dataset.pkl")
    data = pd.read_pickle(filename)
    var_arr = np.zeros((data.shape[0],))+a
    var_ser = pd.Series(var_arr, copy=False,name = "Variance")
    data_var = data.join(var_ser)
    files.append(data_var)
tim_data = pd.concat(files)
print(tim_data.columns)
tim_data['DOC_initial_int'] = round(tim_data.DOC_initial, -3)
cols_to_merge = tim_data.columns.to_list()
#%%
#%%
fd_para_tim_data = pd.merge(fd_data, tim_data[cols_to_merge], on = ["Seed", "Variance", "DOC_initial_int","biomass_species", "carbon_species", "Sim_series"])
#%%

tim_data['DOC_initial_int'] = round(tim_data.DOC_initial, -3)
tim_data['S_initial_int'] = round(tim_data.S_initial, 1)
tim_data['ratio_t_50'] = tim_data.T_50/tim_data.T_50_B1
compl = tim_data.dropna(subset = ['ratio_t_50'])
init_doc_list = np.sort(list(compl.DOC_initial_int.unique()))
activity_list = np.sort(list(compl.activity.unique()))
s_initial_list = np.sort(list(compl.S_initial_int.unique()))
gen_act_tim = tim_data[tim_data.Sim_series=="b_1_all_"]
act_off_tim = tim_data[tim_data.Sim_series!="b_1_all_"]
#%%
for v in [0.01, 0.1, 0.5, 1., 1.5]:
    for s in seed_list:
        for c in [3,6,12,18]:
            for b in [4,8,16,32]:
                sub_scb = para_data[(para_data["Variance"]==v)&(para_data["Seed"].astype(int)==s)&(para_data["carbon_species"].astype(int)==c)&(para_data["biomass_species"].astype(int)==b)]
                #v_mean_base = sub_scb[sub_scb.Sim_series.isin(['b_4_a_', 'b_4_b_','b_4_c_','b_4_d_','b_4_e_'])]['vmax_mean'].mean()
                v_mean_base = sub_scb[sub_scb.Sim_series.isin(['b_1_all_'])]['vmax_mean'].mean()
                para_data.loc[(para_data["Variance"]==v)&(para_data["Seed"].astype(int)==s)&(para_data["carbon_species"].astype(int)==c)&(para_data["biomass_species"].astype(int)==b), "vmax_mean_base"]=v_mean_base
para_data['vmax_ratio'] = para_data.vmax_mean/para_data.vmax_mean_base


for v in [0.01, 0.1, 0.5, 1., 1.5]:
    for s in seed_list:
        for c in [3,6,12,18]:
            for b in [4,8,16,32]:
                for d in init_doc_list:
                    sub_scb = all_data[(all_data["Variance"]==v)&(all_data["Seed"].astype(int)==s)&(all_data["carbon_species"].astype(int)==c)&(all_data["biomass_species"].astype(int)==b)&(all_data["DOC_initial_int"].astype(int)==d)]
                    fd_base = sub_scb[sub_scb.Sim_series=="b_1_all_"]['FD_initial'].values[0]
                    all_data.loc[(all_data["Variance"]==v)&(all_data["Seed"].astype(int)==s)&(all_data["carbon_species"].astype(int)==c)&(all_data["biomass_species"].astype(int)==b)&(all_data["DOC_initial_int"].astype(int)==d), "FD_initial_base"]=fd_base

all_data['DOC_initial_int'] = round(all_data.DOC_initial, -3)
all_data['S_initial_int'] = round(all_data.S_initial, 1)
all_data['NOSC_reldel_maxbio'] = (all_data.NOSC_maxbio/all_data.NOSC_initial - 1) * 100
all_data['FD_reldel_maxbio'] = (all_data.FD_maxbio/all_data.FD_initial - 1) * 100
all_data['Biomass_reldel_maxbio'] = (all_data.Biomass_maxbio/all_data.Biomass_initial - 1) * 100
init_doc_list = np.sort(list(all_data.DOC_initial_int.unique()))
activity_list = np.sort(list(all_data.activity.unique()))
s_initial_list = np.sort(list(all_data.S_initial_int.unique()))
seed_list = np.sort(list(all_data.Seed.unique()))
B2_paras = pd.merge(B2_df, para_data[cols_to_merge], on = ["Seed", "Variance", "biomass_species", "carbon_species", "Sim_series"])
B2_paras.drop_duplicates()
B2_paras_fd = pd.merge(B2_paras, all_data, on = ["Seed", "Variance", "biomass_species", "carbon_species", "Sim_series"])
B2_paras_fd.drop_duplicates()
B2_paras_fd["norm_dec_const"] = 1/B2_paras_fd['ratio_t_50'].to_numpy(dtype = float)
B2_paras_fd["decay_constant"] = 1/B2_paras_fd['t_50_days'].to_numpy(dtype = float)
gen_act = all_data[all_data.Sim_series=="b_1_all_"]
act_off = all_data[all_data.Sim_series!="b_1_all_"]
#%%
#h = sns.jointplot(x = "Variance", y = "FD_initial", hue = "carbon_species", size = "biomass_species", data = gen_act)#hue = 'S_initial_int',
h = sns.scatterplot(x = "Variance", y = "FD_initial", size = "carbon_species", hue = "biomass_species", data = gen_act)#hue = 'S_initial_int',
#h.set_axis_labels("Log normal variance (% of mean)","Initial functional diversity")
h.set_xlabel("Log normal variance (x mean)")
h.set_ylabel("Initial functional diversity")
plt.yscale("log")
plt.legend(bbox_to_anchor=(1,1))
#%%
## OR
sns.boxplot('S_initial_int', 'FD_initial', data = all_data)
plt.yscale("log")
#%%
h = sns.violinplot(x = "carbon_species", y = "NOSC_initial", data = gen_act)#hue = 'S_initial_int',
plt.xlabel("Size of carbon pool")
plt.ylabel("Initial nature of carbon pool")
#%%
sns.jointplot(y = 'CUE_timscale', x = 'NOSC_timscale', hue = 'DOC_initial', data = gen_act)
plt.ylim((-1,1))
#%%
#%%
fig, axes = plt.subplots(1,3,sharex=True,sharey=True,figsize=(10,4))
sns.scatterplot(x = 'NOSC_initial', y = 'NOSC_timscale', hue = 'S_initial_int', data = gen_act[gen_act.DOC_initial_int>=10000], ax = axes[0])
sns.scatterplot(x = 'NOSC_initial', y = 'NOSC_maxbio', hue = 'S_initial_int', data = gen_act[gen_act.DOC_initial_int>=10000], ax=axes[1])
sns.scatterplot(x = 'NOSC_initial', y = 'NOSC_final', hue = 'S_initial_int', data = gen_act[gen_act.DOC_initial_int>=10000], ax = axes[2])
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
        a.legend(bbox_to_anchor=(0.25,-0.2), ncol = 5, title = "H")
#%%
#%%
fig, axes = plt.subplots(1,3,sharex=True,sharey=True,figsize=(10,4))
sns.scatterplot(x = 'CUE_initial', y = 'CUE_timscale', hue = 'S_initial_int', data = gen_act[gen_act.DOC_initial_int<10000], ax = axes[0])
sns.scatterplot(x = 'CUE_initial', y = 'CUE_maxbio', hue = 'S_initial_int', data = gen_act[gen_act.DOC_initial_int<10000], ax=axes[1])
sns.scatterplot(x = 'CUE_initial', y = 'CUE_final', hue = 'S_initial_int', data = gen_act[gen_act.DOC_initial_int<10000], ax = axes[2])
axes[0].set_title ("Characteristic\ntime scale")
axes[1].set_title("Peak biomass")
axes[2].set_title("100 years")
axes[0].set_ylabel("CUE")
axes[0].set_ylim((-1,1))
#axes[0].set_xlim((-1,1))
for a in axes:
    a.set_xlabel("Initial CUE")
    if a!=axes[2]:
        a.get_legend().remove()
    else:
        a.legend(bbox_to_anchor=(0.25,-0.2), ncol = 5, title = "H")
#%%
sns.jointplot(y = 'CUE_maxbio', x = 'NOSC_maxbio', hue = 'DOC_initial', data = gen_act)
plt.ylim(bottom=-0.05)
#%%
gen_act['FD_ratio'] = gen_act.FD_maxbio/gen_act.FD_initial
gen_act['FDxcarbon'] = gen_act.FD_initial*gen_act.DOC_initial
sns.scatterplot(x = "FD_initial", y = "FD_ratio", hue = 'DOC_initial', style = "Variance", data = gen_act)
plt.xscale("log")
plt.xlabel("Initial functional diversity x Carbon availability")
plt.ylabel("Functional diversity: Peak/Initial")
plt.legend(bbox_to_anchor=(1, 1))
#%%
sns.kdeplot(x="FDxcarbon", y="FD_ratio", shade = 'True', thresh = 0.005,data = gen_act)
plt.xlabel("Initial functional diversity x Carbon availability")
plt.ylabel("Gain in functional diversity: Peak/Initial")
plt.legend(bbox_to_anchor=(1, 1))
#%%
gen_act['Biomass_ratio'] = gen_act.Biomass_maxbio/gen_act.Biomass_initial
sns.scatterplot(x = "FDxcarbon", y = "Biomass_ratio", hue = 'DOC_initial', style = "Variance", data = gen_act)
plt.xscale("log")
plt.xlabel("Initial functional diversity x Carbon availability")
plt.ylabel("Gain in biomass: Peak/Initial")
plt.legend(bbox_to_anchor=(1, 1))
#%%
sns.kdeplot(x="FDxcarbon", y="Biomass_ratio", shade = 'True', thresh = 0.005,data = gen_act)
plt.xlabel("Initial functional diversity x Carbon availability")
plt.ylabel("Gain in biomass: Peak/Initial")
plt.legend(bbox_to_anchor=(1, 1))
#%%
sns.scatterplot(x = "S_initial", y = "FD_initial", hue = 'carbon_species', style = "Variance", data = gen_act)
plt.yscale("log")
plt.xlabel("Initial Shannon diversity index")
plt.ylabel("Initial functional diversity")
plt.legend(bbox_to_anchor=(1, 1))
#%%
sns.scatterplot(x = "Biomass_ratio", y = "FD_ratio", hue = 'DOC_initial_int', style = "Variance", size = "carbon_species",data = gen_act)
plt.axhline(y=1, c = 'red', label = "Same community")
plt.axvline(x=1, c = 'black', label = 'Dying community')
plt.yscale("log")
plt.xlabel("Ratio: Biomass")
plt.ylabel("Ratio: Functional diversity")
plt.legend(bbox_to_anchor=(1, 1))
#%%
#%%
sns.kdeplot(x = "Biomass_ratio", y = "FD_ratio", shade = True, thresh = 0.005,data = gen_act)
plt.axhline(y=1, c = 'red', label = "Same community")
plt.axvline(x=1, c = 'black', label = 'Dying community')
#plt.yscale("log")
plt.xlabel("Ratio: Biomass")
plt.ylabel("Ratio: Functional diversity")
plt.legend(bbox_to_anchor=(1, 1))
### LOSS IN ACTIVITY
#%%
act_off['FD_ratio'] = act_off.FD_maxbio/act_off.FD_initial
act_off['FDxcarbon'] = act_off.FD_initial*act_off.DOC_initial*100/(act_off.activity*act_off.biomass_species*act_off.carbon_species)
sns.scatterplot(x = "FDxcarbon", y = "FD_ratio", hue = 'DOC_initial_int', style = "Variance", data = act_off)
plt.xscale("log")
plt.xlabel("Initial functional diversity x carbon availability")
plt.ylabel("Functional diversity: Peak/Initial")
plt.legend(bbox_to_anchor=(1, 1))
#%%
sns.kdeplot(x = "FD_initial", y = "FD_ratio", shade = True,data = act_off)
plt.xlabel("Initial functional diversity x carbon availability")
plt.ylabel("Functional diversity: Peak/Initial")
plt.legend(bbox_to_anchor=(1, 1))
#%%
act_off['Biomass_ratio'] = act_off.Biomass_maxbio/act_off.Biomass_initial
sns.scatterplot(x = "FD_initial", y = "Biomass_ratio", hue = 'DOC_initial', style = "Variance", data = act_off)
plt.xscale("log")
plt.xlabel("Initial functional diversity")
plt.ylabel("Gain in biomass: Peak/Initial")
plt.legend(bbox_to_anchor=(1, 1))
#%%
sns.kdeplot(x = "FD_initial", y = "Biomass_ratio", shade = True,data = act_off)
plt.xlabel("Initial functional diversity")
plt.ylabel("Gain in biomass: Peak/Initial")
plt.legend(bbox_to_anchor=(1, 1))
#%%
g = sns.FacetGrid(act_off, col="DOC_initial_int", row="activity",hue = 'carbon_species', palette = "Blues")
g.map(sns.scatterplot,'FD_initial','Biomass_ratio')
plt.xscale("log")
#%%
sns.scatterplot(x = "Biomass_ratio", y = "FD_ratio", hue = 'activity', style = "Variance", size = "carbon_species",data = act_off)
plt.axhline(y=1, c = 'red', label = "Same community")
plt.axvline(x=1, c = 'black', label = 'Dying community')
plt.yscale("log")
plt.xlabel("Ratio: Biomass")
plt.ylabel("Ratio: Functional diversity")
plt.legend(bbox_to_anchor=(1, 1))
#%%

#%%
## Predicting decay constant
compl["x_var"] = compl['S_initial']*compl['carbon_species']*compl["activity"]/100
compl["y_var"] = 1/compl['t_50_days'].to_numpy(dtype = float)    
g = sns.FacetGrid(compl, row="Variance", col="activity",hue = 'DOC_initial_int', palette = "rocket_r")
g.map(sns.scatterplot,'x_var','y_var')
#%%
## Decay constant as a function of initial diversity?
## General communities first
gen_fdiv_tim = pd.merge(gen_act_tim, gen_act, on = ['Seed', 'Variance', 'Sim_series','carbon_species', 'biomass_species', 'DOC_initial_int'])
gen_fdiv_tim["x_var"] = gen_fdiv_tim['S_initial_x']*gen_fdiv_tim['carbon_species']*gen_fdiv_tim["activity_x"]/100
gen_fdiv_tim["y_var"] = 1/gen_fdiv_tim['t_50_days'].to_numpy(dtype = float)    
#%%
sns.scatterplot('FD_initial', 'y_var', hue = 'DOC_initial_int', style = 'Variance', data = gen_fdiv_tim)
plt.xscale("log")
plt.ylabel("Decay constant (/day)")
plt.xlabel("Initial functional diversity")
plt.legend(bbox_to_anchor=(1,1))
#%%
actoff_fdiv_tim = pd.merge(act_off_tim, act_off, on = ['Seed', 'Variance', 'Sim_series','carbon_species', 'biomass_species', 'DOC_initial_int'])
actoff_fdiv_tim["x_var"] = actoff_fdiv_tim['S_initial_x']*actoff_fdiv_tim['carbon_species']*actoff_fdiv_tim["activity_x"]/100
actoff_fdiv_tim["y_var"] = 1/actoff_fdiv_tim['t_50_days'].to_numpy(dtype = float)    
sns.scatterplot('FD_initial', 'y_var', hue = 'DOC_initial_int', style = 'Variance', data = actoff_fdiv_tim)
plt.xscale("log")
plt.ylabel("Decay constant (/day)")
plt.xlabel("Initial functional diversity")
plt.legend(bbox_to_anchor=(1,1))
#%%
actoff_fdiv_tim["FD_initial_ratio"] = actoff_fdiv_tim.FD_initial/actoff_fdiv_tim.FD_initial_base
actoff_fdiv_tim["y_var"] = 1/actoff_fdiv_tim['ratio_t_50'].to_numpy(dtype = float)
sns.scatterplot('FD_initial','y_var', data = actoff_fdiv_tim, style = "Variance", alpha = 0.5, hue = "DOC_initial_int")
plt.xscale("log")
plt.yscale("log")
plt.ylabel("Normalized decay constant (-)")
plt.xlabel("Initial functional diversity")
plt.legend(bbox_to_anchor = (1,1))
#%%
actoff_fdiv_tim["FD_initial_ratio"] = actoff_fdiv_tim.FD_initial/actoff_fdiv_tim.FD_initial_base
sns.scatterplot('FD_initial_ratio','y_var', data = actoff_fdiv_tim, style = "Variance", alpha = 0.5, hue = "DOC_initial_int")
plt.axhline(y=1, c = 'red', label = 'No change in\ndecay constant')
plt.axvline(x=1, c = 'black', label = 'No change in\nfunctional diversity')
plt.xscale("log")
plt.yscale("log")
plt.ylabel("Normalized decay constant (-)")
plt.xlabel("Normalized functional diversity")
plt.legend(bbox_to_anchor = (1,1))
#%%
g = sns.FacetGrid(actoff_fdiv_tim, col="Variance",col_wrap=2,hue = 'DOC_initial_int', palette = "Blues")
g.map(sns.scatterplot,'FD_initial_ratio','y_var')
g.set(yscale="log")
g.set(xscale="log")
g.add_legend()
#%%
actoff_fdiv_tim["FD_cov"]= actoff_fdiv_tim.FD_initial/actoff_fdiv_tim.vmax_mean
#%%

sns.scatterplot('FD_ratio','y_var', data = actoff_fdiv_tim, style = "Variance", alpha = 0.5, hue = "DOC_initial_int")
#plt.ylim(bottom=1)
#plt.xscale("log")
#plt.yscale("log")
plt.ylabel("Normalized decay constant (-)")
plt.xlabel("Peak functional diversity")
plt.legend(bbox_to_anchor = (1,1))

#%%
compl["x_var"] = (compl['S_initial']*compl['carbon_species']*compl["activity"]/100)
compl["y_var"] = 1/compl['ratio_t_50'].to_numpy(dtype = float)
sns.scatterplot('x_var','y_var', data = compl, style = "Variance", alpha = 0.5, hue = "DOC_initial_int")
#plt.ylim(bottom=1)
plt.yscale("log")
plt.legend(bbox_to_anchor = (1,1))
## In some scenarios, decay constant is higher than
## the base case. Drill into what makes this happen.

## PREDICT FUNCTION PARAMETERS OF DECAY CONSTANT ##
#%%
def powlaw(x, a, b):
    return a * (x**b)

def sig(x, k, c):
    return c + 1/(1+((x)**k))**(1/k)

def glf (x, a, k, c, q, b, v):
    return a + (k-a)/(c + q* np.exp(-b*x))**(1/v)

def glf_simp (x, a):
    return (1 + np.exp(-x))**(-a)

def sqf (x, a, c):
    return a * np.sqrt(x) + c

def sig_exp(x,k,c):
    return c + 1/(1+np.exp(-k*x))
#%%
gen_act_tim["x_var"] = (gen_act_tim['S_initial']*gen_act_tim['carbon_species']*gen_act_tim["activity"]/100)
gen_act_tim["y_var"] = 1/gen_act_tim['t_50_days'].to_numpy(dtype = float)
gen_act_tim = gen_act_tim.sort_values(by = 'x_var')
#%%
sns.scatterplot('x_var','y_var', data = gen_act_tim, style = "Variance", alpha = 0.5, hue = "DOC_initial_int")
#plt.ylim(bottom=1)
#plt.yscale("log")
plt.legend(bbox_to_anchor = (1,1))
for doc in [1000,2000,5000,10000,15000]:
    X = gen_act_tim[gen_act_tim.DOC_initial_int==doc]["x_var"]
    y = gen_act_tim[gen_act_tim.DOC_initial_int==doc]['y_var']
    meany = np.mean(y)
    popt_func, pcov_e = curve_fit(powlaw, xdata=X, ydata=y, p0=[0.4,0.7])
    ypred = powlaw(X, popt_func[0], popt_func[1])
    yerr = 1 - (np.sum((ypred-y)**2)/np.sum((y - meany)**2))
    plt.plot(X, ypred, 'r-', alpha = 0.6, label = doc)
    print(yerr)
plt.legend()

#%%
## PREDICT FUNCTION PARAMETERS FOR IMPACT ON DECAY CONSTANT ##
#%%
non_100 = compl[compl.activity<100]
non_100["exp_x_var"] = non_100['S_initial']*non_100["activity"]/100
non_100 = non_100.sort_values(by=["exp_x_var"])
X = non_100['exp_x_var']
y = 1/non_100['ratio_t_50'].to_numpy(dtype = float)
meany = np.mean(y)
plt.scatter(X,y, marker = '.', label = "data", alpha = 0.1)
plt.yscale("log")
plt.ylabel("Normalized decay constant")
plt.xlabel("Active Shannon diversity")
popt_func, pcov_e = curve_fit(powlaw, xdata=X, ydata=y, p0=[0.4,0.7])
ypred = powlaw(X, popt_func[0], popt_func[1])
yerr = 1 - (np.sum((ypred-y)**2)/np.sum((y - meany)**2))
plt.plot(X, ypred, 'r-', alpha = 0.6, label = 'exp_func')
plt.text(1.5, 0.04, "p1: "+str(round(popt_func[0],2)))
plt.text(1.5, 0.02, "p2: "+str(round(popt_func[1],2)))
#plt.text(1.5, 0.18, "p3: "+str(round(popt_func[2],2)))
plt.text(1.5, 0.01, "R2: "+str(round(yerr,2)))
plt.savefig(os.path.join(project_dir, "exponential_predict_impact_decay_constant_compet_adapt.png"), dpi = 300)
## ALTERNATIVE BASELINE?
#%%
seed_list = list(tim_data.Seed.unique())
c_sp_list = list(tim_data.carbon_species.unique())
bio_sp_list = list(tim_data.biomass_species.unique())
row = []
for v in [0.01, 0.1, 0.5, 1., 1.5]:
    sub0 = tim_data[tim_data.Variance==v]
    for seed in seed_list:
        sub1 = sub0[sub0.Seed==seed]
        for doci in init_doc_list:
            sub2 = sub1[sub1.DOC_initial_int==doci]
            for c_sp in c_sp_list:
                sub3 = sub2[sub2.carbon_species==c_sp]
                for b_sp in bio_sp_list:
                    sub4 = sub3[sub3.biomass_species==b_sp]
                    #sub = all_data[(all_data['Seed']==seed)&]
                    sub50 = sub4[sub4.activity==50]
                    median_t_50_B2 = np.median(sub50["T_50"])
                    row.append([seed, v, doci, c_sp, b_sp, median_t_50_B2])
T_50_B2 = pd.DataFrame.from_records(row, columns = ["Seed", "Variance", "DOC_initial_int", "carbon_species", "biomass_species", "T_50_B2"])

#%%
B2_df = pd.merge(tim_data, T_50_B2, on = ["Seed", "Variance", "DOC_initial_int", "carbon_species", "biomass_species"]).reset_index()
B2_df = B2_df.drop(['index'], axis = 1)
B2_df["ratio_t_50_b2"] = B2_df["T_50"]/B2_df["T_50_B2"]
B2_df["dec_const_b2_norm"] = 1/B2_df.ratio_t_50_b2
B2_df = B2_df.dropna()
#%%
B2_df["x_b2"] = B2_df.S_initial*(B2_df.activity/100 - 50/100)
B2_df["x_b2_v2"] = B2_df.S_initial
B2_df["x_b2_v3"] = (B2_df.activity/100 - 50/100)
plt.scatter(data = B2_df, x="x_b2", y="dec_const_b2_norm",marker = '.', alpha = 0.1)
plt.yscale("log")
plt.ylabel("Normalized decay constant")
plt.xlabel("Change in active Shannon diversity")
#plt.savefig(os.path.join(figures_dir, "impact_decay_constant_compet_adapt_S_active_B2.png"), dpi = 300, bbox_inches='tight')

### PREDICT IMPACT OF CHANGING DIVERSITY on DECAY CONSTANT B2
#%%
def logistic_func(x, L, b):
    return L/(1+np.exp(-b*(x)))
B2_df["exp_x_var"] = B2_df.S_initial*(B2_df.activity/100 - 50/100)#/B2_df.DOC_initial#
B2_df = B2_df.sort_values(by=["exp_x_var"])
#%%
fig, axes = plt.subplots(3,2, figsize=(8,10), sharex = True)
ax = axes.flatten()
for v, a in zip([0.01, 0.1,0.5,1.,1.5],[0,1,2,3,4]):
    X = B2_df[B2_df.Variance==v]['exp_x_var']
    y = B2_df[B2_df.Variance==v]["dec_const_b2_norm"]
    mean = np.mean(y)
    ax[a].scatter(X,y, marker = '.', label = "Data", alpha = 0.1)
    #ax[a].ylabel("Normalized decay constant")
    #ax[a].xlabel("Change in activity of microbial community given\nShannon diversity from 50% active microbes baseline")
    popt_func, pcov_e = curve_fit(logistic_func, xdata=X, ydata=y, p0=[2, 0.5])
    ypred = logistic_func(X, popt_func[0], popt_func[1])
    yerr = 1 - (np.sum((ypred-y)**2)/np.sum((y - meany)**2))
    ax[a].plot(X, ypred, 'r-', alpha = 0.6, label = 'sig_func')
    #plt.yscale("log")
    ax[a].set_title("Variance: "+ str(round(v*100))+"%")
    ax[a].text(-1.2, 5.0, "L: "+str(round(popt_func[0],2)))
    ax[a].text(-1.2, 4.5, "b: "+str(round(popt_func[1],2)))
    ax[a].text(-1.2, 4.0, "R2: "+str(round(yerr,2)))
plt.savefig(os.path.join(project_dir, "logistic_func_predict_impact_decay_constant_compet_adapt_B2.png"), dpi = 300, bbox_inches = 'tight')
## the response of all communities seem to be similar
##(coefficients are similar with varying performance
## (R2))
#%%
X = B2_df['exp_x_var']
y = B2_df["dec_const_b2_norm"]
mean = np.mean(y)
sns.scatterplot(x="exp_x_var", y = "dec_const_b2_norm", alpha = 0.5, hue = "DOC_initial_int", style="Variance", data = B2_df)
plt.ylabel("Normalized decay constant")
plt.xlabel("Change in activity of microbial community given\nShannon diversity from 50% active microbes baseline")
popt_func, pcov_e = curve_fit(logistic_func, xdata=X, ydata=y, p0=[2, 1.5])
ypred = logistic_func(X, popt_func[0], popt_func[1])
yerr = 1 - (np.sum((ypred-y)**2)/np.sum((y - meany)**2))
plt.plot(X, ypred, 'r-', alpha = 0.6, label = 'sig_func')
plt.text(-1.2, 8.0, "L: "+str(round(popt_func[0],2)))
plt.text(-1.2, 6.5, "b: "+str(round(popt_func[1],2)))
plt.text(-1.2, 5.0, "R2: "+str(round(yerr,2)))
plt.legend(bbox_to_anchor = (1,1))
plt.savefig(os.path.join(project_dir, "logistic_func_predict_impact_decay_constant_compet_adapt_B2_condensed.png"), dpi = 300, bbox_inches = 'tight')
#%%


#%%
sns.scatterplot(x="vmax_ratio", y = "norm_dec_const", hue = "DOC_initial_int_x",style = "Variance", data = B2_paras_fd)
plt.xscale("log")
plt.yscale("log")
plt.legend(bbox_to_anchor=(1,1))
#%%
#%%
compl_paras = pd.merge(compl, para_data[cols_to_merge], on = ["Seed", "Variance", "biomass_species", "carbon_species", "Sim_series"])
compl_paras.drop_duplicates()
compl_paras_fd = pd.merge(compl_paras, all_data, on = ["Seed", "Variance", "biomass_species", "carbon_species", "Sim_series"])
compl_paras_fd.drop_duplicates()
compl_paras_fd["y_var"] = 1/compl_paras_fd['ratio_t_50'].to_numpy(dtype = float)
#%%
sns.scatterplot(x="vmax_ratio", y = "y_var", hue = "DOC_initial_int_x",style = "activity_x", data = compl_paras_fd)
plt.xscale("log")
plt.yscale("log")
plt.legend(bbox_to_anchor=(1,1))
#%%
actoff_fd_paras = pd.merge(actoff_fdiv_tim, para_data[cols_to_merge], on = ["Seed", "Variance", "biomass_species", "carbon_species", "Sim_series"])
actoff_fd_paras.drop_duplicates()
#%%
sns.scatterplot(x="vmax_ratio", y = "y_var", hue = "Variance",style = "activity_x", data = actoff_fd_paras)
plt.xscale("log")
plt.yscale("log")
plt.legend(bbox_to_anchor=(1,1))
#%%
sns.scatterplot(y="FD_initial_ratio", x = "vmax_ratio", hue = "activity_x", style = "Variance", data = actoff_fd_paras)
plt.yscale("log")
plt.xscale("log")
plt.xlabel("Normalized vmax")
plt.ylabel("Normalized functional diversity")
#%%
g = sns.FacetGrid(actoff_fd_paras, col="Variance", row = "activity_x")
g.map(sns.scatterplot,'FD_initial_ratio','vmax_ratio')
g.set(yscale="log")
g.set(xscale="log")
g.add_legend()
#%%
actoff_fd_paras["FD_cov"] = np.sqrt(actoff_fd_paras.FD_initial)/actoff_fd_paras.vmax_mean
#%%
sns.scatterplot(y="FD_initial", x = "vmax_median", hue = "activity_x", style = "Variance",data =actoff_fd_paras)
plt.yscale("log")
plt.xscale("log")
#%%
sns.scatterplot(y="FD_ratio", x = "FD_cov", hue = "DOC_initial_int", size = "activity_x", style = "Variance",data =actoff_fd_paras)
plt.xscale("log")
plt.legend(bbox_to_anchor=(1,1))
#%%
#%%
sns.scatterplot(y="Biomass_ratio", x = "FD_cov", hue = "DOC_initial_int", size = "activity_x", style = "Variance",data =actoff_fd_paras)
plt.xscale("log")
plt.legend(bbox_to_anchor=(1,1))
#%%
actoff_fd_paras["y_var"] = 1/actoff_fd_paras['ratio_t_50'].to_numpy(dtype = float)
sns.scatterplot(y="y_var", x = "FD_cov", hue = "DOC_initial_int", size = "activity_x", style = "Variance",data =actoff_fd_paras)
plt.yscale("log")
plt.xscale("log")
plt.legend(bbox_to_anchor=(1,1))
#%%
g = sns.FacetGrid(actoff_fd_paras, col="activity_x", row="Variance",hue = 'DOC_initial_int', palette = "Blues")
g.map(sns.scatterplot,'FD_cov','y_var')
g.set(yscale="log")
g.set(xscale="log")
g.add_legend()

#%%
faster = compl_paras_fd[compl_paras_fd.y_var>1]
#%%
low_variance_faster = faster[faster.Variance<0.5]
high_variance_faster = faster[faster.Variance>1.0]
mid_variance_faster = faster[(faster.Variance>=0.5) & (faster.Variance<=1.0)]

fig, axes = plt.subplots(3,1,figsize =(4,10))#sharex = True,
sns.scatterplot('vmax_ratio','y_var', data = low_variance_faster, size = "activity_x", alpha = 0.5, hue = "DOC_initial_int_x", ax = axes[0])
sns.scatterplot('vmax_ratio','y_var', data = mid_variance_faster, size = "activity_x", alpha = 0.5, hue = "DOC_initial_int_x", ax = axes[1])
sns.scatterplot('vmax_ratio','y_var', data = high_variance_faster, size = "activity_x", alpha = 0.5, hue = "DOC_initial_int_x", ax = axes[2])
axes[0].set_title ("Similar communities")
axes[1].set_title ("Transitioning communites")
axes[2].set_title ("Dissimilar communities")
#plt.ylim(bottom=1)
#plt.yscale("log")
for a in axes.flatten():
    a.legend(bbox_to_anchor = (1,1))

#%%
low_variance_faster = B2_paras_fd[B2_paras_fd.Variance<0.5]
high_variance_faster = B2_paras_fd[B2_paras_fd.Variance>1.0]
mid_variance_faster = B2_paras_fd[(B2_paras_fd.Variance>=0.5) & (B2_paras_fd.Variance<=1.)]

fig, axes = plt.subplots(3,1,figsize =(4,10),sharex = True)#,sharey=True)
sns.scatterplot('exp_x_var','dec_const_b2_norm', data = low_variance_faster, size = "FD_initial", alpha = 0.5, hue = "DOC_initial_x", ax = axes[0])
sns.scatterplot('exp_x_var','dec_const_b2_norm', data = mid_variance_faster, size = "FD_initial", alpha = 0.5, hue = "DOC_initial_x", ax = axes[1])
sns.scatterplot('exp_x_var','dec_const_b2_norm', data = high_variance_faster, size = "FD_initial", alpha = 0.5, hue = "DOC_initial_x", ax = axes[2])
axes[0].set_title ("Similar communities")
axes[1].set_title ("Transitioning communites")
axes[2].set_title ("Dissimilar communities")
#plt.ylim(bottom=1)
#plt.yscale("log")
for a in [0,1]:
    axes[a].get_legend().remove()
axes[2].legend(bbox_to_anchor = (1,1))

#%%
sns.scatterplot('exp_x_var','dec_const_b2_norm', data = B2_paras_fd, size = "FD_initial", alpha = 0.5, hue = "DOC_initial_x")
plt.legend(bbox_to_anchor=(1,1))
### OLD
#%%
gen_act['NOSC_ratio'] = gen_act.NOSC_maxbio/gen_act.NOSC_initial
sns.scatterplot(x = "NOSC_initial", y = "NOSC_maxbio", size = "S_initial_int", alpha = 0.5,hue = 'DOC_initial_int', style ="Variance", data = gen_act)
plt.xlabel("Initial NOSC")
plt.ylabel("NOSC at peak biomass")
plt.legend(bbox_to_anchor=(1, 1))
#%%
sns.scatterplot(x = "NOSC_initial", y = "NOSC_maxbio", size = "S_initial_int",alpha = 0.5, hue ="Variance", data = gen_act[gen_act.DOC_initial_int==15000])
plt.xlabel("Initial NOSC")
plt.ylabel("End state NOSC")
plt.legend(bbox_to_anchor=(1, 1))
#%%
#subset = all_data[all_data.DOC_initial==10000]
eqline = np.sort(all_data['FD_initial'])
plt.plot(eqline, eqline, 'r-')
all_data['hue_var'] = np.log10(all_data.DOC_initial/all_data.carbon_species)
sns.scatterplot(x = "FD_initial", y = "FD_maxbio", hue = 'DOC_initial', data = all_data)
#plt.xscale("log")
#plt.yscale("log")
#%%
subset = all_data[all_data.DOC_initial==10000]
eqline = np.sort(subset['FD_initial'])
plt.plot(eqline, eqline, 'r-')
sns.scatterplot(x = "FD_initial", y = "FD_final", hue = 'carbon_species', style ="Variance", data = subset)
plt.xscale("log")
plt.yscale("log")
#%%
eqline = np.sort(all_data["NOSC_initial"])
sns.scatterplot(x = "NOSC_initial", y = "NOSC_maxbio", hue = 'DOC_initial', style ="Variance", data = all_data)
plt.plot(eqline, eqline, 'r-')
#%%
all_data['hue_var'] = all_data.biomass_species*all_data.activity/100#*all_data.carbon_species))
sns.scatterplot(x = "NOSC_initial", y = "NOSC_maxbio", hue = 'hue_var', size = 'DOC_initial', data = all_data)
plt.plot(eqline, eqline, 'r-')
#%%
#sns.scatterplot(x = "Biomass_reldel_maxbio", y = "NOSC_reldel_maxbio", hue = 'DOC_initial', data = all_data)
sns.scatterplot(x = 'FD_reldel_maxbio', y = 'Biomass_reldel_maxbio', hue = 'S_initial', data = all_data)
#%%
all_data['x_var'] = all_data.S_initial*all_data.activity*all_data.DOC_initial
all_data['y_var'] = all_data.Biomass_reldel_maxbio
all_data['hue_var'] = 10**(all_data.FD_initial)
sns.scatterplot(x = 'x_var', y = 'y_var', hue = 'hue_var', size = 'hue_var', data = all_data)
plt.xscale("log")
#%%
subset = all_data[(all_data['DOC_initial']==2000)&(all_data['S_initial']<4)]
subset['hue_var'] = subset.activity#*subset.DOC_initial*subset.S_initial
subset['y_var'] = subset.Biomass_reldel_maxbio
subset['x_var'] = subset.FD_initial
sns.scatterplot(x = 'x_var', y = 'y_var', hue = 'hue_var', size = 'hue_var', data = subset)
plt.xscale("log")

#%%
all_data['hue_var'] = np.log10(all_data.DOC_initial/(all_data.carbon_species*all_data.biomass_species))
all_data["size_var"] = all_data.activity
sns.scatterplot(x = "FD_initial", y = "Biomass_reldel_maxbio", hue = "hue_var", size = "size_var",data = all_data)
plt.xscale("log")
#plt.yscale("log")
#%%

#%%
sns.scatterplot('carbon_species', 'NOSC_initial', hue = "Seed",data = all_data)
#%%
subset = all_data[all_data.carbon_species==12]
eqline = np.sort(subset['NOSC_initial'])
plt.plot(eqline, eqline, 'r-')
sns.scatterplot(x = "NOSC_initial", y = "NOSC_maxbio", hue = 'FD_maxbio', data = subset)

#%%



#%%
reg_line = mlines.Line2D([], [], linestyle = '-', color = "red", marker = None, label='Regression')

funcres = []
compl["X_regr"] = compl['S_initial']*compl['carbon_species']*compl["activity"]/100
compl = compl.sort_values(by = "X_regr")
X = compl.X_regr
y = 1/compl['t_50_days'].to_numpy(dtype = float)
size_var = 1/compl['t_50_days'].to_numpy(dtype = float)
fig, ax1 = plt.subplots(1,1)
plt.scatter(X,y, c = compl["DOC_initial"], cmap = 'YlGnBu')
plt.ylabel("Decay constant (1/day)")
plt.xlabel("Active S x #C pools")
for doci in init_doc_list:
    sub = compl[(compl['DOC_initial_int']==doci)]
    X_sub = sub.X_regr
    y_sub = 1/sub['t_50_days'].to_numpy(dtype = float)
    meany = np.mean(y_sub)
    size_var = 1/sub['t_50_days'].to_numpy(dtype = float)
    try:
        popt_func, pcov_e = curve_fit(powlaw, xdata=X_sub, ydata=y_sub, p0=[0.0001, 0.5])
        ypred = powlaw(X_sub, popt_func[0], popt_func[1])
        yerr = 1-np.sum((ypred-y_sub)**2)/np.sum((y_sub - meany)**2)
        plt.plot(X_sub, ypred, 'r-', alpha = 0.6)
    except:
        popt_func = [np.nan, np.nan]
        yerr = np.nan
    funcres.append([doci, popt_func[0], popt_func[1], yerr])
funcresdf = pd.DataFrame.from_records(funcres, columns = ["DOC_initial", "a", "b", "r2"])

sm = plt.cm.ScalarMappable(cmap='YlGnBu')
sm.set_array([])
cbar = plt.colorbar(sm, pad = 0.25, orientation = 'horizontal')
cbar.ax.get_xaxis().set_ticks([])
for j, lab in enumerate(init_doc_list):
    cbar.ax.text(j/4, -0.5, lab, ha='center', va='center')
cbar.ax.get_xaxis().labelpad = -27
cbar.ax.set_xlabel('Initial available carbon (uM)')#, rotation = 180)
legend1 = plt.legend(handles = [reg_line], bbox_to_anchor = (0.2,-0.2))
ax1.add_artist(legend1)
plt.savefig(os.path.join(figures_dir, "powlaw_decay_const_compet_adapt_condense.png"), dpi = 300, bbox_inches = 'tight')

fig, ax = plt.subplots(1,1)
g = sns.scatterplot(data = funcresdf, x = "DOC_initial", y = "a", size = "r2", hue = "b")
leg_handles, leg_labels = g.get_legend_handles_labels()
new_leg = leg_labels
for t in leg_labels:
    if leg_labels.index(t) in [1,2,3,4,5,7,8,9,10,11]:
        new_leg[leg_labels.index(t)] = str(np.round(np.float64(t),2))
    else:
        pass
new_leg[0] = "Exponent"
new_leg[6] = "Score"
plt.legend(leg_handles, new_leg, bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
g.set_ylabel("Scaling factor")
g.set_xlabel("Initial available carbon (uM)")
plt.savefig(os.path.join(figures_dir, "powlaw_metrics_decay_const_compet_adapt_condense.png"), dpi = 300, bbox_inches = 'tight')


#%%
fig, axes = plt.subplots(5,1, sharex=True, figsize = (8,20), sharey = True)
axe = axes.flatten()
for doc_i in init_doc_list:
    sub = generalist_act[generalist_act.DOC_initial_int==doc_i]
    ax1 = axe[list(init_doc_list).index(doc_i)]
    sns.boxplot(x = sub.carbon_species, y = sub.decay_const, hue = sub.S_initial_int, palette = "Greens", ax = ax1)
    ax1.set_ylabel("Decay constant")
    ax1.set_title(str(doc_i) + " uM")
    ax1.set_xlabel("")
    if list(init_doc_list).index(doc_i)<4:
        ax1.legend("")
ax1.legend(bbox_to_anchor=(0.75,-0.2), title = "Shannon diversity", ncol = 4)
plt.xlabel("Number of Carbon species")
plt.savefig(os.path.join(figures_dir, "decay_constant_compet_adapt_initial_DOC_#c_#S_box_full_active.png"), dpi = 300)
#%%
h = sns.boxplot(x = generalist_act.S_initial_int, y = generalist_act.decay_const, hue = generalist_act.DOC_initial_int)
plt.xlabel("Shannon diversity")
plt.ylabel("Decay constant")
#plt.savefig(os.path.join(figures_dir, "decay_constant_null_initial_S_box_full_active.png"), dpi = 300)

#%%
plt.figure(figsize = (6,2))
h = sns.violinplot(x = generalist_act.S_initial_int, y = generalist_act.decay_const)#, hue = generalist_act.DOC_initial_int)
plt.xlabel("Shannon diversity")
plt.ylabel("Decay constant")
#plt.savefig(os.path.join(figures_dir, "decay_constant_null_initial_S_violin_full_active.png"), dpi = 300)

#%%
#%%
func_div = pd.read_csv(os.path.join(results_dir, "competition_adaptation_parameters.csv"))

plt.figure()
sns.scatterplot(data = func_div, x = "vmax_mean", y= "k_mean", hue = "carbon_species")
#%%
for s in seed_list:
    for c in [3,6,12,18]:
        for b in [4,8,16,32]:
            sub_scb = func_div[(func_div["Seed"].astype(int)==s)&(func_div["carbon_species"].astype(int)==c)&(func_div["biomass_species"].astype(int)==b)]
            v_mean_base = sub_scb[sub_scb.Sim_series.isin(['b_4_a_', 'b_4_b_','b_4_c_','b_4_d_','b_4_e_'])]['vmax_mean'].median()
            func_div.loc[(func_div["Seed"].astype(int)==s)&(func_div["carbon_species"].astype(int)==c)&(func_div["biomass_species"].astype(int)==b), "vmax_mean_base"]=v_mean_base
#%%
cols_to_merge = func_div.columns.to_list()
B2_paras = pd.merge(B2_df, func_div[cols_to_merge], on = ["Seed", "biomass_species", "carbon_species", "Sim_series"])
B2_paras.drop_duplicates()
#%%
sns.jointplot(data = B2_paras, x = 'k_mean', y = 'vmax_mean', kind = 'kde', hue ='carbon_species')
#%%
B2_paras["x_var"] = B2_paras.DOC_initial*B2_paras.carbon_species*B2_paras.vmax_mean
B2_paras["decay_const"] = 1/B2_df.t_50_days
#%%
plt.figure()
sns.scatterplot(y = 'decay_const', x = "x_var", hue = "biomass_species", data = B2_paras)
plt.xscale("log")
plt.xlabel("Available carbon x carbon species number x average vmax")
plt.legend( bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
#%%
ylineqbase = np.sort(B2_paras.vmax_mean_base)
sns.scatterplot(y = "vmax_mean", x = "vmax_mean_base", hue = "activity", style = "biomass_species", data = B2_paras)
plt.plot(ylineqbase, ylineqbase, 'r-')
plt.ylabel("Vmax_mean")
plt.legend( bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.savefig(os.path.join(figures_dir, "vmax_base_scenario_comp_adapt.png"), dpi = 300, bbox_inches = 'tight')
#%%
B2_paras["reldelvmax"] = B2_paras.vmax_mean/B2_paras.vmax_mean_base
sub_off = B2_paras#[B2_paras['Sim_series']!='b_1_all_']
sub_off["ratio_xvar"] = (1 - sub_off.reldelvmax)
sub_off["yplot"]  = sub_off.dec_const_b2_norm#1/sub_off.ratio_t_50
sns.scatterplot(y = 'yplot', x = "ratio_xvar", hue = "DOC_initial_int", style = "S_initial_int",data = sub_off)
plt.xlabel("Reduction in decomposition potential (difference in vmean)\nnormalized by decomposition potential of base case(vmean)")
plt.ylabel("Decay constant normalized by\nthat of the base case")
plt.legend( bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.savefig(os.path.join(figures_dir, "decay_const_vmax_scenario_comp_adapt.png"), dpi = 300, bbox_inches = 'tight')

#%%
sns.scatterplot(y = 'decay_const', x = "x_var", hue = "DOC_initial_int", style = "biomass_species",data = B2_paras)
plt.xlabel("Available carbon x carbon species number x biomass weighted vmax")
plt.ylabel("Decay constant")
plt.legend( bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
#%%
sns.scatterplot(y = "ratio_xvar", x = "activity", hue = "DOC_initial_int", style = "biomass_species", data = sub_off)
plt.axhline(y=0, c = 'r', label = "baseline")
plt.ylabel("Reduction in decomposition potential (difference in vmean)\nnormalized by decomposition potential of base case(vmean)")
plt.legend( bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.savefig(os.path.join(figures_dir, "decay_const_activity_scenario_comp_adapt.png"), dpi = 300, bbox_inches = 'tight')