#%%
import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

import matplotlib.lines as mlines

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient", "gen_spec_lognorm_1_5x")
#project_dir = os.path.join("C:/", "Users", "swkh9804", "Documents", "Project_data", "HoliSoils", "transient", "gen_spec_lognorm_1_5x")
simulations_dir = os.path.join(project_dir, "simulations")
results_dir = os.path.join(project_dir, "results")
figures_dir = os.path.join(project_dir, "figures")

all_data = pd.read_pickle(os.path.join(results_dir, "competition_adaptation_carbon__loss_0.9_cue_combined_dataset.pkl"))
print(all_data.columns)
all_data['DOC_initial_int'] = round(all_data.DOC_initial, -3)
all_data['S_initial_int'] = round(all_data.S_initial, 1)
all_data['NOSC_reldel_maxbio'] = (all_data.NOSC_maxbio/all_data.NOSC_initial - 1) * 100
all_data['FD_reldel_maxbio'] = (all_data.FD_maxbio/all_data.FD_initial - 1) * 100
all_data['Biomass_reldel_maxbio'] = (all_data.Biomass_maxbio/all_data.Biomass_initial - 1) * 100
init_doc_list = np.sort(list(all_data.DOC_initial_int.unique()))
activity_list = np.sort(list(all_data.activity.unique()))
s_initial_list = np.sort(list(all_data.S_initial_int.unique()))

#%%
sns.jointplot(y = 'CUE_maxbio', x = 'NOSC_maxbio', hue = 'DOC_initial', data = all_data)
#%%
sns.scatterplot(x = "FD_initial", y = "FD_maxbio", hue = 'DOC_initial', data = all_data)
#%%
subset = all_data[all_data.DOC_initial==10000]
eqline = np.sort(subset['FD_initial'])
plt.plot(eqline, eqline, 'r-')
sns.scatterplot(x = "FD_initial", y = "FD_maxbio", hue = 'carbon_species', data = subset)
plt.xscale("log")
plt.yscale("log")
#%%
subset = all_data[all_data.DOC_initial==10000]
eqline = np.sort(subset['FD_initial'])
plt.plot(eqline, eqline, 'r-')
sns.scatterplot(x = "FD_initial", y = "FD_final", hue = 'carbon_species', data = subset)
plt.xscale("log")
plt.yscale("log")
#%%
eqline = np.sort(all_data["NOSC_initial"])
sns.scatterplot(x = "NOSC_initial", y = "NOSC_maxbio", hue = 'DOC_initial', data = all_data)
plt.plot(eqline, eqline, 'r-')
#%%
#sns.scatterplot(x = "Biomass_reldel_maxbio", y = "NOSC_reldel_maxbio", hue = 'DOC_initial', data = all_data)
sns.scatterplot(x = 'FD_reldel_maxbio', y = 'Biomass_reldel_maxbio', hue = 'S_initial', data = all_data)
#%%
all_data['x_var'] = all_data.FD_initial#/all_data.carbon_species
all_data['y_var'] = all_data.Biomass_reldel_maxbio#/all_data.DOC_initial
sns.scatterplot(x = 'x_var', y = 'y_var', hue = 'carbon_species', data = all_data)
plt.yscale("log")
plt.xscale("log")
#plt.ylim((-500,500))
#%%
sns.scatterplot('carbon_species', 'NOSC_initial', hue = 'Seed', data = all_data)
#%%
subset = all_data[all_data.carbon_species==12]
eqline = np.sort(subset['NOSC_initial'])
plt.plot(eqline, eqline, 'r-')
sns.scatterplot(x = "NOSC_initial", y = "NOSC_maxbio", hue = 'FD_maxbio', data = subset)

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
#plt.savefig(os.path.join(figures_dir, "decay_constant_predict_comp_adapt_panel.png"), dpi = 300, bbox_inches = 'tight')

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
plt.savefig(os.path.join(figures_dir, "decay_constant_predict_comp_adapt_condense.png"), dpi = 300)
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
#plt.savefig(os.path.join(figures_dir, "impact_reaction_time_scale_comp_ada_panel.png"), dpi = 300, bbox_inches = 'tight')

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
plt.savefig(os.path.join(figures_dir, "impact_decay_constant_comp_ada_condense.png"), dpi = 300)

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
plt.savefig(os.path.join(figures_dir, "impact_reaction_time_scale_comp_ada_condense.png"), dpi = 300)

#%%
X_comb = (compl['carbon_species']*compl['biomass_species']*compl["activity"]/100)
y_dist = (1/compl['t_50_days']).to_numpy(dtype = float)
sns.displot(x = X_comb, y = y_dist, hue = compl["DOC_initial_int"],kind = "kde")
plt.xlabel("Active B-C connections")
plt.ylabel("Decay constant")
plt.savefig(os.path.join(figures_dir, "decay_constant_comp_adapt_kde.png"), dpi = 300)
#%%
X_comb = (non_100['carbon_species']*non_100['biomass_species']*non_100["activity"]/100)
y_dist = (1/non_100['ratio_t_50']).to_numpy(dtype = float)
sns.displot(x = X_comb, y = y_dist, hue = non_100["DOC_initial_int"],kind = "kde")
plt.xlabel("Active B-C connections")
plt.ylabel("Mormalised decay constant")
plt.savefig(os.path.join(figures_dir, "impact_decay_constant_comp_adapt_kde.png"), dpi = 300)
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
plt.savefig(os.path.join(figures_dir, "exponential_predict_impact_decay_constant_compet_adapt.png"), dpi = 300)

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
plt.savefig(os.path.join(figures_dir, "decay_constant_comp_adapt_initial_DOC_violin.png"), dpi = 300)

#%%
h = sns.violinplot(x = generalist_act.DOC_initial_int, y = generalist_act.decay_const, hue = generalist_act.carbon_species)
plt.xlabel("Initial available carbon (uM)")
plt.ylabel("Decay constant")
plt.savefig(os.path.join(figures_dir, "decay_constant_comp_adapt_initial_DOC_violin_full_active.png"), dpi = 300)

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
B2_df = B2_df.dropna()
#%%
B2_df["x_b2"] = B2_df.S_initial*(B2_df.activity/100 - 50/100)
B2_df["x_b2_v2"] = B2_df.S_initial
B2_df["x_b2_v3"] = (B2_df.activity/100 - 50/100)
plt.scatter(data = B2_df, x="x_b2", y="dec_const_b2_norm")
plt.ylabel("Normalized decay constant")
plt.xlabel("Change in active Shannon diversity")
#plt.savefig(os.path.join(figures_dir, "impact_decay_constant_compet_adapt_S_active_B2.png"), dpi = 300, bbox_inches='tight')

### PREDICT IMPACT OF CHANGING DIVERSITY on DECAY CONSTANT B2
#%%
B2_df["exp_x_var"] = B2_df.S_initial*(B2_df.activity/100 - 50/100)#/B2_df.DOC_initial#
B2_df = B2_df.sort_values(by=["exp_x_var"])
X = B2_df['exp_x_var']
y = B2_df["dec_const_b2_norm"]
meany = np.mean(y)
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
plt.savefig(os.path.join(figures_dir, "logistic_func_predict_impact_decay_constant_compet_adapt_B2.png"), dpi = 300, bbox_inches = 'tight')

#%%
func_div = pd.read_csv(os.path.join(results_dir, "competition_adaptation_parameters.csv"))
fig, axes = plt.subplots(3,3, sharex = 'col', figsize = (12,6))
sns.scatterplot(data = func_div, x = 'vmax_mean', y = 'vmax_median', hue = 'biomass_species', ax = axes[0,0])
axes[0,0].set_title("vmax")
sns.scatterplot(data = func_div, x = 'vmax_mean', y = 'vmax_std', hue = 'biomass_species', ax = axes[1,0])
sns.scatterplot(data = func_div, x = 'vmax_mean', y = 'vmax_skew', hue = 'biomass_species', ax = axes[2,0])
sns.scatterplot(data = func_div, x = 'k_mean', y = 'k_median', hue = 'biomass_species', ax = axes[0,1])
axes[0,1].set_title("Half-sat-const")
sns.scatterplot(data = func_div, x = 'k_mean', y = 'k_std', hue = 'biomass_species', ax = axes[1,1])
sns.scatterplot(data = func_div, x = 'k_mean', y = 'k_skew', hue = 'biomass_species', ax = axes[2,1])
sns.scatterplot(data = func_div, x = 'm_mean', y = 'm_median', hue = 'biomass_species', ax = axes[0,2])
axes[0,2].set_title("Mortality")
sns.scatterplot(data = func_div, x = 'm_mean', y = 'm_std', hue = 'biomass_species', ax = axes[1,2])
sns.scatterplot(data = func_div, x = 'm_mean', y = 'm_skew', hue = 'biomass_species', ax = axes[2,2])
for ax in axes.flatten()[:-1]:
    ax.legend([],[], frameon=False)
for ax in axes.flatten()[:]:
    ax.set_ylabel('')
    ax.set_xlabel('')
for ax in axes[2,:]:
    ax.set_xlabel("Mean")
axes[0,0].set_ylabel("Median")
axes[1,0].set_ylabel("Sdev")
axes[2,0].set_ylabel("Skew")
plt.savefig(os.path.join(figures_dir, "paras_mean_distribution.png"), dpi = 300, bbox_inches = 'tight')

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