#%%
import os
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

import matplotlib.lines as mlines

## LOAD RESULTS

project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient","activity_loss_-02")
simulations_dir = os.path.join(project_dir, "simulations")
results_dir = os.path.join(project_dir, "results")
figures_dir = os.path.join(project_dir, "figures")

all_data = pd.read_pickle(os.path.join(results_dir, "competition_adaptation_carbon__loss_0.9_combined_dataset.pkl"))
print(all_data.columns)
all_data['DOC_initial_int'] = round(all_data.DOC_initial, -3)
all_data['S_initial_int'] = round(all_data.S_initial, 1)
all_data['ratio_t_50'] = all_data.T_50/all_data.T_50_B1
compl = all_data.dropna(subset = ['ratio_t_50'])
init_doc_list = np.sort(list(compl.DOC_initial_int.unique()))
activity_list = np.sort(list(compl.activity.unique()))
s_initial_list = np.sort(list(compl.S_initial_int.unique()))
non_100 = compl[compl.activity<100]

#%%
from mpl_toolkits import mplot3d
fig = plt.figure(figsize=(7,7))
ax = plt.axes(projection = '3d')
X = compl["DOC_initial"]
y = compl['carbon_species']*compl['S_initial']*compl["activity"]/100
z = 1/compl['t_50_days'].to_numpy(dtype = float)
size_var = 1/compl['t_50_days'].to_numpy(dtype = float)
ax.scatter(X,y, z, c = np.log10(size_var), cmap = 'YlGnBu')
ax.set_xlabel("Initial DOC (uM)")
ax.set_ylabel("Active B-C connection")
ax.set_zlabel("Decay constant")
plt.legend()

#%%
fig, axes = plt.subplots(1,5, figsize = (10,2), sharex = True, sharey = True)
for i in list(range(5)):
    sub = compl[compl.DOC_initial_int == init_doc_list[i]]
    ax = axes.flatten()[i]
    X = sub['S_initial']*sub['carbon_species']*sub["activity"]/100
    y = 1/sub['t_50_days'].to_numpy(dtype = float)
    size_var = 1/sub['t_50_days'].to_numpy(dtype = float)
    ax.scatter(X,y, c = np.log10(size_var), cmap = 'YlGnBu')
    ax.set_title(str(init_doc_list[i])+" (uM)")
    ax.set_xlabel("Active B-C connection \n(Sx#Cxactivity%)")
axes.flatten()[0].set_ylabel("Decay constant (1/day)")
plt.legend()
#plt.savefig(os.path.join(figures_dir, "decay_constant_predict_null_panel.png"), dpi = 300, bbox_inches = 'tight')

#%%
X = (compl['S_initial']*compl['carbon_species']*compl["activity"]/100)#compl["DOC_initial"]/
y = 1/compl['t_50_days'].to_numpy(dtype = float)
size_var = 1/compl['t_50_days'].to_numpy(dtype = float)
plt.scatter(X,y, c = compl["DOC_initial"], cmap = 'YlGnBu')
plt.ylabel("Decay constant (1/day)")
plt.xlabel("Active Sx #C pools")
#plt.xscale("log")
#plt.legend(title = "Available C")
sm = plt.cm.ScalarMappable(cmap="YlGnBu")
sm.set_array([])
cbar = plt.colorbar(sm)
cbar.ax.get_yaxis().set_ticks([])
for j, lab in enumerate(init_doc_list):
    cbar.ax.text(3.0, j/4.2, lab, ha='center', va='center')
    cbar.ax.get_yaxis().labelpad = -20
    cbar.ax.set_ylabel('Initial available carbon (uM) ', rotation=90)
#plt.savefig(os.path.join(figures_dir, "decay_constant_predict_null_condense.png"), dpi = 300)
#%%
fig, axes = plt.subplots(1,5, figsize = (10,2), sharex = True, sharey = True)
for i in list(range(5)):
    sub = non_100[non_100.DOC_initial_int == init_doc_list[i]]
    ax = axes.flatten()[i]
    X = sub['S_initial']*sub['carbon_species']*sub["activity"]/100
    y = 1/sub['ratio_t_50'].to_numpy(dtype = float)
    ax.scatter(X,y, c = np.log10(1/sub['t_50_days']), cmap = 'YlGnBu')
    ax.set_xlabel("Active B-C connection \n(Sx#Cxactivity%)")
    ax.set_title(str(init_doc_list[i])+" (uM)")
axes.flatten()[0].set_ylabel("Normalized decay \nconstant (1/day)")
#plt.savefig(os.path.join(figures_dir, "impact_decay_constant_comp_adapt_panel.png"), dpi = 300, bbox_inches = 'tight')
#%%
fig, axes = plt.subplots(1,5, figsize = (10,2), sharex = True, sharey = True)
for i in list(range(5)):
    sub = non_100[non_100.DOC_initial_int == init_doc_list[i]]
    ax = axes.flatten()[i]
    X = sub['S_initial']*sub['carbon_species']*sub["activity"]/100
    y = sub['ratio_t_50'].to_numpy(dtype = float)
    ax.scatter(X,y, c = np.log10(1/sub['t_50_days']), cmap = 'YlGnBu')
    ax.set_xlabel("Active B-C connection \n(Sx#Cxactivity%)")
    ax.set_title(str(init_doc_list[i])+" (uM)")
axes.flatten()[0].set_ylabel("Normalized reaction \ntime scale (day)")
plt.legend()
#plt.savefig(os.path.join(figures_dir, "impact_reaction_time_scale_null_panel.png"), dpi = 300, bbox_inches = 'tight')

#%%
X = (non_100['S_initial']*non_100['carbon_species']*non_100["activity"]/100)
y = 1/non_100['ratio_t_50'].to_numpy(dtype = float)
plt.scatter(X,y, c = non_100["DOC_initial"], cmap = 'YlGnBu')
plt.ylabel("Normalised decay constant (1/day)")
plt.xlabel("Active B-C connection")
sm = plt.cm.ScalarMappable(cmap="YlGnBu")
sm.set_array([])
cbar = plt.colorbar(sm)
cbar.ax.get_yaxis().set_ticks([])
for j, lab in enumerate(init_doc_list):
    cbar.ax.text(3.0, j/4.2, lab, ha='center', va='center')
    cbar.ax.get_yaxis().labelpad = -20
    cbar.ax.set_ylabel('Initial available carbon (uM) ', rotation=90)
#plt.savefig(os.path.join(figures_dir, "impact_decay_constant_null_condense.png"), dpi = 300)

#%%
X = (compl['S_initial']*compl['carbon_species']*compl["activity"]/100)
y = compl['ratio_t_50'].to_numpy(dtype = float)
plt.scatter(X,y, c = compl["DOC_initial"], cmap = 'YlGnBu')
plt.ylabel("Normalised reaction time scale (day)")
plt.xlabel("Active B-C connections")
sm = plt.cm.ScalarMappable(cmap="YlGnBu")
sm.set_array([])
cbar = plt.colorbar(sm)
cbar.ax.get_yaxis().set_ticks([])
for j, lab in enumerate(init_doc_list):
    cbar.ax.text(3.0, j/4.2, lab, ha='center', va='center')
    cbar.ax.get_yaxis().labelpad = -20
    cbar.ax.set_ylabel('Initial available carbon (uM) ', rotation=90)
#plt.savefig(os.path.join(figures_dir, "impact_reaction_time_scale_null_condense.png"), dpi = 300)

#%%
X_comb = (compl['carbon_species']*compl['biomass_species']*compl["activity"]/100)
y_dist = (1/compl['t_50_days']).to_numpy(dtype = float)
sns.displot(x = X_comb, y = y_dist, hue = compl["DOC_initial_int"],kind = "kde")
plt.xlabel("Active B-C connections")
plt.ylabel("Decay constant")
#plt.savefig(os.path.join(figures_dir, "decay_constant_null_kde.png"), dpi = 300)
#%%
X_comb = (non_100['carbon_species']*non_100['biomass_species']*non_100["activity"]/100)
y_dist = (1/non_100['ratio_t_50']).to_numpy(dtype = float)
sns.displot(x = X_comb, y = y_dist, hue = non_100["DOC_initial_int"],kind = "kde")
plt.xlabel("Active B-C connections")
plt.ylabel("Mormalised decay constant")
#plt.savefig(os.path.join(figures_dir, "impact_decay_constant_null_kde.png"), dpi = 300)
#%%
X_comb = (compl['carbon_species']*compl['biomass_species']*compl["activity"]/100)
y_dist = (1/compl['t_50_days']).to_numpy(dtype = float)
h = sns.jointplot(x=X_comb, y=y_dist, hue=compl["DOC_initial_int"])
h.set_axis_labels("Active B-C connections", "Decay constant")

#%%
X_comb = (non_100['carbon_species']*non_100['biomass_species']*non_100["activity"]/100)
y_dist = (1/non_100['ratio_t_50']).to_numpy(dtype = float)
h = sns.jointplot(x=X_comb, y=y_dist, hue=non_100["DOC_initial_int"])
h.set_axis_labels("Active B-C connections", "Normalized decay constant")

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

## PREDICT FUNCTION PARAMETERS FOR IMPACT ON DECAY CONSTANT ##
#%%
non_100["exp_x_var"] = non_100['S_initial']*non_100["activity"]/100
non_100 = non_100.sort_values(by=["exp_x_var"])
X = non_100['exp_x_var']
y = 1/non_100['ratio_t_50'].to_numpy(dtype = float)
meany = np.mean(y)
plt.scatter(X,y, marker = '.', label = "data")
plt.ylabel("Normalized decay constant")
plt.xlabel("Active Shannon diversity")
popt_func, pcov_e = curve_fit(powlaw, xdata=X, ydata=y, p0=[-0.8,-0.1])
ypred = powlaw(X, popt_func[0], popt_func[1])
yerr = 1 - (np.sum((ypred-y)**2)/np.sum((y - meany)**2))
plt.plot(X, ypred, 'r-', alpha = 0.6, label = 'exp_func')
plt.text(1.5, 0.3, "p1: "+str(round(popt_func[0],2)))
plt.text(1.5, 0.24, "p2: "+str(round(popt_func[1],2)))
#plt.text(1.5, 0.18, "p3: "+str(round(popt_func[2],2)))
plt.text(1.5, 0.12, "R2: "+str(round(yerr,2)))
#plt.savefig(os.path.join(figures_dir, "exponential_predict_impact_decay_constant_compet_adapt.png"), dpi = 300)

#%%
generalist_act = compl[compl.activity==100].reset_index()
generalist_act["decay_const"] = 1/generalist_act["t_50_days"]
print(generalist_act.t_50_days.min(), generalist_act.t_50_days.max())
print(generalist_act.S_initial.min(), generalist_act.S_initial.max())
#%%
compl["decay_const"] = 1/compl.t_50_days
h = sns.violinplot(x = compl.DOC_initial_int, y = compl.decay_const, hue = compl.carbon_species)
plt.xlabel("Initial available carbon (uM)")
plt.ylabel("Decay constant")
#plt.savefig(os.path.join(figures_dir, "decay_constant_comp_adapt_initial_DOC_violin.png"), dpi = 300)

#%%
h = sns.violinplot(x = generalist_act.DOC_initial_int, y = generalist_act.decay_const, hue = generalist_act.carbon_species)
plt.xlabel("Initial available carbon (uM)")
plt.ylabel("Decay constant")
#plt.savefig(os.path.join(figures_dir, "decay_constant_null_initial_DOC_violin_full_active.png"), dpi = 300)

#%%
plt.figure(figsize = (8,4))
h = sns.boxplot(x = generalist_act.DOC_initial_int, y = generalist_act.decay_const, hue = generalist_act.S_initial_int)
plt.xlabel("Initial available carbon (uM)")
plt.ylabel("Decay constant")

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
#plt.savefig(os.path.join(figures_dir, "decay_constant_compet_adapt_initial_DOC_#c_#S_box_full_active.png"), dpi = 300)
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
## ALTERNATIVE BASELINE?
#%%
seed_list = list(all_data.Seed.unique())
c_sp_list = list(all_data.carbon_species.unique())
bio_sp_list = list(all_data.biomass_species.unique())
row = []
for seed in seed_list:
    sub1 = all_data[all_data.Seed==seed]
    for doci in init_doc_list:
        sub2 = sub1[sub1.DOC_initial_int==doci]
        for c_sp in c_sp_list:
            sub3 = sub2[sub2.carbon_species==c_sp]
            for b_sp in bio_sp_list:
                sub4 = sub3[sub3.biomass_species==b_sp]
                #sub = all_data[(all_data['Seed']==seed)&]
                sub50 = sub4[sub4.activity==50]
                median_t_50_B2 = np.median(sub50["T_50"])
                row.append([seed, doci, c_sp, b_sp, median_t_50_B2])
T_50_B2 = pd.DataFrame.from_records(row, columns = ["Seed", "DOC_initial_int", "carbon_species", "biomass_species", "T_50_B2"])

#%%
B2_df = pd.merge(all_data, T_50_B2, on = ["Seed", "DOC_initial_int", "carbon_species", "biomass_species"]).reset_index()
B2_df = B2_df.drop(['index'], axis = 1)
B2_df["ratio_t_50_b2"] = B2_df["T_50"]/B2_df["T_50_B2"]
B2_df["dec_const_b2_norm"] = 1/B2_df.ratio_t_50_b2

#%%
B2_df["x_b2"] = B2_df.S_initial*(B2_df.activity/100 - 50/100)
B2_df["x_b2_v2"] = B2_df.S_initial
B2_df["x_b2_v3"] = (B2_df.activity/100 - 50/100)
plt.scatter(data = B2_df, x="x_b2", y="dec_const_b2_norm")
plt.ylabel("Normalized decay constant")
plt.xlabel("Change in active Shannon diversity")
#plt.savefig(os.path.join(figures_dir, "impact_decay_constant_compet_adapt_S_active_B2.png"), dpi = 300, bbox_inches='tight')

#%%
B2_df = B2_df.dropna()
### PREDICT IMPACT OF CHANGING DIVERSITY on DECAY CONSTANT B2
#%%

B2_df["exp_x_var"] = B2_df.S_initial*(B2_df.activity/100 - 50/100)#/B2_df.DOC_initial#
B2_df = B2_df.sort_values(by=["exp_x_var"])
X = B2_df['exp_x_var']
y = B2_df["dec_const_b2_norm"]
meany = np.mean(y)
plt.scatter(X,y, marker = '.', label = "data")
plt.ylabel("Normalized decay constant")
plt.xlabel("Change in activity of microbial community given\nShannon diversity from 50% active microbes baseline")
#try:
popt_func, pcov_e = curve_fit(glf, xdata=X, ydata=y, p0=[0.1, 2, 1, 1, 0.1, 0.5])
ypred = glf(X, popt_func[0], popt_func[1], popt_func[2], popt_func[3], popt_func[4], popt_func[5])
yerr = 1 - (np.sum((ypred-y)**2)/np.sum((y - meany)**2))
plt.plot(X, ypred, 'r-', alpha = 0.6, label = 'sig_func')
plt.text(-1.2, 4.6, "a: "+str(round(popt_func[0],2)))
plt.text(-1.2, 4.3, "k: "+str(round(popt_func[1],2)))
plt.text(-1.2, 4.0, "c: "+str(round(popt_func[2],2)))
plt.text(-1.2, 3.7, "q: "+str(round(popt_func[3],2)))
plt.text(-1.2, 3.4, "b: "+str(round(popt_func[4],2)))
plt.text(-1.2, 3.1, "v: "+str(round(popt_func[5],2)))
plt.text(-1.2, 2.8, "R2: "+str(round(yerr,2)))
#plt.savefig(os.path.join(figures_dir, "glf_predict_impact_decay_constant_compet_adapt_B2.png"), dpi = 300, bbox_inches = 'tight')
#%%
def logistic_func(x, L, b):
    return L/(1+np.exp(-b*(x)))

plt.scatter(X,y, marker = '.', label = "data", alpha = 0.1)
plt.ylabel("Normalized decay constant")
plt.xlabel("Change in activity of microbial community given\nShannon diversity from 50% active microbes baseline")
popt_func, pcov_e = curve_fit(logistic_func, xdata=X, ydata=y, p0=[2, 0.5])
ypred = logistic_func(X, popt_func[0], popt_func[1])
yerr = 1 - (np.sum((ypred-y)**2)/np.sum((y - meany)**2))
plt.plot(X, ypred, 'r-', alpha = 0.6, label = 'sig_func')
plt.text(-1.2, 4.6, "L: "+str(round(popt_func[0],2)))
plt.text(-1.2, 4.3, "b: "+str(round(popt_func[1],2)))
plt.text(-1.2, 3.7, "R2: "+str(round(yerr,2)))
#plt.savefig(os.path.join(figures_dir, "logistic_func_predict_impact_decay_constant_compet_adapt_B2.png"), dpi = 300, bbox_inches = 'tight')

#%%
base_sims = B2_df[B2_df['exp_x_var'] == 0.0]
print(base_sims.Sim_series.unique())
#%%
sns.boxplot(y = 'dec_const_b2_norm', x = "S_initial_int", data = base_sims)
#%%
wide_var = B2_df[B2_df['S_initial_int']==1.4]
sns.boxplot(y = 'dec_const_b2_norm', x = "Seed", data = base_sims)
plt.xticks(rotation=90)
#%%
seed_var = B2_df[B2_df['Seed'].isin([610229235, 643338060, 714504443, 983307757])]
sns.boxplot(y = 'dec_const_b2_norm', x = "Sim_series", data = base_sims)
#%%
import h5py
def gather_paras(sim_data):
    dom_n, bio_n = sim_data['species_number']
    init_bio = sim_data['initial_conditions']['biomass']
    init_bio_sum = np.sum(init_bio)
    paras = sim_data['parameters']
    dying = np.tile(init_bio*np.asarray(paras['mortality_const'])/init_bio_sum, (dom_n, 1)).flatten().reshape(-1,1)
    exo_enz = np.tile(init_bio*np.asarray(paras['exo_enzyme_rate'])/init_bio_sum, (dom_n, 1)).flatten().reshape(-1,1)
    v_params = (init_bio*np.asarray((paras['max_rate']))/init_bio_sum).flatten().reshape(-1,1)*np.asarray(paras['carbon_uptake']).flatten().reshape(-1,1)
    k_params = v_params*np.asarray(paras['half_saturation']).flatten().reshape(-1,1)
    ox = np.asarray(list([x]*bio_n for x in paras['oxidation_state'])).flatten().reshape(-1,1)
    paras_arr = np.append(np.append(np.append(dying, exo_enz, axis=1),np.append(v_params, k_params, axis = 1), axis = 1), ox, axis = 1)
    return paras_arr
#%%
seed_sim_list = seed_list#[610229235, 643338060, 714504443, 983307757]
cn_list = [3,6,12,18]
bio_n_series = [4]
init_dom_list = [1000,2000,5000,10000,15000]
filestring = "competition_adaptation_carbon_"

all_para_arr = np.zeros((0,9))
for c_n in cn_list:
    row = []
    c_b_row = []
    results_filename = os.path.join(results_dir, filestring + str(c_n))

    for seed_sim in seed_sim_list:
        # Load all datasets and save their Shannon and diversity indices in a dataframe
        seed_all = 'seed_'+str(seed_sim)
        
        details_subfolder = filestring + str(c_n) + '_'+str(seed_sim) + '_ip_0'
        simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
        hrw = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'r+')
        
        for b_n in bio_n_series:
            for t_dom_initial in init_dom_list:
                c_b = "bio_n_"+ str(b_n)
                dom_init = "dom_initial_" + str(t_dom_initial)
                doc_input = (t_dom_initial)
                seed_arr = np.zeros((c_n*b_n,1))+seed_sim
                carbon_arr = np.zeros((c_n*b_n,1))+c_n
                seed_carbon = np.append(seed_arr, carbon_arr, axis = 1)
                for base, act in zip(["b_4"], [0.5]):
                    act_arr = np.zeros((c_n*b_n,1))+act
                    seed_sim_arr = np.append(seed_carbon, act_arr, axis = 1)
                    if base == "b_1":
                        base_case = "b_1_all_"
                        sim_data = hrw[base_case][c_b][dom_init][seed_all]
                        paras_sim = np.append(seed_sim_arr, gather_paras(sim_data), axis = 1)
                        all_para_arr = np.append(all_para_arr, paras_sim, axis = 0)
                    else:
                        for label,label_id in zip(["a", "b", "c","d","e"], [0,1,2,3,4]):
                            sim = base + "_" + label + "_"
                            level_arr = np.zeros((c_n*b_n,1))+label_id
                            sim_data = hrw[sim][c_b][dom_init][seed_all]
                            seed_sim_level = np.append(seed_sim_arr, level_arr, axis = 1)
                            paras_sim = np.append(seed_sim_level, gather_paras(sim_data), axis = 1)
                            all_para_arr = np.append(all_para_arr, paras_sim, axis = 0)
        hrw.close()
#%%
para_df = pd.DataFrame(all_para_arr, columns = ['Seed', 'carbon_species','Activity', 'level_id', 'm', 'exo_enz','vmax', 'k', 'oxidation_state'])
para_df.loc[para_df["level_id"] == 0.0, "Sim_series"] = 'b_4_a_'
para_df.loc[para_df["level_id"] == 1.0, "Sim_series"] = 'b_4_b_'
para_df.loc[para_df["level_id"] == 2.0, "Sim_series"] = 'b_4_c_'
para_df.loc[para_df["level_id"] == 3.0, "Sim_series"] = 'b_4_d_'
para_df.loc[para_df["level_id"] == 4.0, "Sim_series"] = 'b_4_e_'
#%%
func_div = para_df.groupby(['Seed', 'carbon_species', 'Sim_series', 'level_id'])['vmax', 'k', 'm', 'exo_enz'].var().reset_index()
func_div["Seed"] = func_div["Seed"].astype(int)
func_div["carbon_species"] = func_div["carbon_species"].astype(int)
#%%
seed_paras = pd.merge(B2_df, func_div[["Seed", "carbon_species", "Sim_series", "level_id","m","exo_enz", "vmax", "k"]], on = ["Seed", "carbon_species", "Sim_series"])
seed_paras.drop_duplicates()
#%%
seed_paras['num_cum'] = seed_paras.exo_enz*seed_paras.vmax
sns.boxplot(y = 'num_cum', x = "Seed", data = seed_paras)
#%%
sns.jointplot(data = seed_paras, x = 'k', y = 'num_cum', kind = 'kde', hue ='carbon_species')
#%%
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
seed_paras["log_doc_initial_int"] = np.log(seed_paras["DOC_initial"])
seed_paras["vmax100000"] = seed_paras["vmax"]*100000
X = seed_paras[['vmax100000', 'm', 'exo_enz']]#,'log_doc_initial_int']]
X = StandardScaler().fit_transform(X)
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(X)
#%%
principalDf = pd.DataFrame(data = principalComponents
             , columns = ['principal component 1', 'principal component 2'])
pca_df = pd.concat([seed_paras, principalDf], axis = 1)
#%%
sns.scatterplot(data = pca_df.sort_values(by=["ratio_t_50"]), x = 'principal component 1', y = 'principal component 2',
hue = pca_df['dec_const_b2_norm']*10)
#%%
from mpl_toolkits.axes_grid1 import make_axes_locatable
components = pca.components_.T
colmax = np.abs(components).max()
fig,ax = plt.subplots(1,1)
im = ax.imshow(components, cmap="RdBu_r", vmax=colmax, vmin=-colmax)
ax.set_yticks(np.arange(0,components.shape[0],1).astype(int))
ax.set_yticklabels(['vmax', 'm', 'exo_enz'])#,'DOC_initial_int']) #'DOC_initial_int',
ax.set_xticklabels(['', "Comp. 1", "Comp. 2"])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
#%%

#%%
sns.scatterplot(y = 'dec_const_b2_norm', x = "carbon_species", hue = "Sim_series", data = seed_paras)
#%%
seed_paras["v_k"] = seed_paras.vmax/seed_paras.k
g = sns.FacetGrid(seed_paras, col="Sim_series",  row="DOC_initial_int", hue = "carbon_species", palette = "tab10")
g.map(sns.scatterplot, "num_cum",  "dec_const_b2_norm")
g.add_legend()
#%%
vkl2 = seed_paras[seed_paras['v_k']<2e-6]
g = sns.FacetGrid(vkl2, col="Sim_series",  row="DOC_initial_int", hue = "carbon_species", palette = "tab10")
g.map(sns.scatterplot, "v_k",  "T_50")
g.add_legend()
#%%
sns.jointplot(para_df["m"], para_df["vmax"])
#%%
sns.jointplot(para_df["m"], para_df["vmax"], kind = 'kde', hue = para_df['Activity'])