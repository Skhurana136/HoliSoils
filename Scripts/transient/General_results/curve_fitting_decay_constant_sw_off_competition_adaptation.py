#%%
import os
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit
from scipy.stats import gamma, beta
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
def linres(x, a, c):
    return a * x + c

def expres(x, a, b, c):
    return a * np.exp(b * x) + c

def sig(x, k, c):
    return c + 1/(1+((x)**k))**(1/k)

def glf (x, a, k, c, q, b, v):
    return a + (k-a)/(c + q* np.exp(-b*x))**(1/v)

def sqf (x, a, c):
    return a * np.sqrt(x) + c

def sig_exp(x,k,c):
    return c + 1/(1+np.exp(-k*x))

#%%
funcres = []
fig, axes = plt.subplots(1,5, figsize = (10,2), sharex = True, sharey = True)
ax = axes.flatten()
idx = 0
for doci in init_doc_list:
    sub = compl[(compl['DOC_initial_int']==doci)].sort_values(by=['S_initial'])
    X = sub['S_initial']*sub['carbon_species']*sub["activity"]/100#.to_numpy(dtype = float)
    y = 1/ sub['t_50_days'].to_numpy(dtype = float)
    meany = np.mean(y)
    b = sub["Biomass_initial"].to_numpy(dtype = float)[0]
    axes[idx].scatter(X,y, marker = '.', label = "data")
    axes[idx].set_title(str(doci) + " uM")
    try:
        popt_func, pcov_e = curve_fit(linres, xdata=X, ydata=y, p0=[0.0001,0.0001])
        ypred = linres(X, popt_func[0], popt_func[1])
        yerr = 1 - (np.sum((ypred-y)**2)/np.sum((y - meany)**2))
        axes[idx].plot(X, ypred, 'r-', alpha = 0.6, label = 'lin_func')
    except:
        popt_func = [np.nan, np.nan]
        yerr = np.nan
    funcres.append([doci, b, popt_func[0], popt_func[1],yerr])
    idx=idx+1
ax[0].set_ylabel("Decay constant\n(1/day)")
for a in ax[:]:
    a.set_xlabel("Active B-C connections")
funcresdf = pd.DataFrame.from_records(funcres, columns = ["DOC_initial", "Biomass_initial", "a", "c", "r2"])
plt.show()

plt.figure()
g = sns.scatterplot(data = funcresdf, x = "DOC_initial", y = "a", size = "r2", hue = "c")
g.set_ylabel("Slope")
g.set_xlabel("Initial available carbon (uM)")
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
#plt.savefig(os.path.join(figures_dir, "decay_const_reg_paras_null.png"), dpi = 300, bbox_inches = 'tight')

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
    X_sub =sub.X_regr
    y_sub = 1/sub['t_50_days'].to_numpy(dtype = float)
    meany = np.mean(y_sub)
    size_var = 1/sub['t_50_days'].to_numpy(dtype = float)
    try:
        popt_func, pcov_e = curve_fit(expres, xdata=X_sub, ydata=y_sub, p0=[-0.01,-1, -0.0001])
        X_sub_sort = X_sub.sort_values()
        ypred = expres(X_sub_sort, popt_func[0], popt_func[1], popt_func[2])
        yerr = 1-np.sum((ypred-y_sub)**2)/np.sum((y_sub - meany)**2)
        plt.plot(X_sub_sort, ypred, 'r-', alpha = 0.6)
    except:
        popt_func = [np.nan, np.nan, np.nan]
        yerr = np.nan
    funcres.append([doci,b, popt_func[0], popt_func[1],popt_func[2], yerr])
funcresdf = pd.DataFrame.from_records(funcres, columns = ["DOC_initial", "Biomass_initial", "a", "b", "c", "r2"])
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
#plt.savefig(os.path.join(figures_dir, "decay_const_null_condense.png"), dpi = 300, bbox_inches = 'tight')

fig, ax = plt.subplots(1,1)
g = sns.scatterplot(data = funcresdf, x = "DOC_initial", y = "a", size = "b", hue = "c")
leg_handles, leg_labels = g.get_legend_handles_labels()
new_leg = leg_labels
for t in leg_labels:
    if leg_labels.index(t) in [1,2,3,4,5,7,8,9,10,11]:
        new_leg[leg_labels.index(t)] = str(np.round(np.float64(t),3))
    else:
        pass
plt.legend(leg_handles, new_leg, bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.text(1300,-0.002, "r2: "+str(round(funcresdf[funcresdf.DOC_initial==1000]['r2'].values[0],2)))
plt.text(2100,-0.003, "r2: "+str(round(funcresdf[funcresdf.DOC_initial==2000]['r2'].values[0],2)))
plt.text(5100,-0.007, "r2: "+str(round(funcresdf[funcresdf.DOC_initial==5000]['r2'].values[0],2)))
plt.text(10100,-0.011, "r2: "+str(round(funcresdf[funcresdf.DOC_initial==10000]['r2'].values[0],2)))
plt.text(13000,-0.013, "r2: "+str(round(funcresdf[funcresdf.DOC_initial==15000]['r2'].values[0],2)))
g.set_ylabel("Slope")
g.set_xlabel("Initial available carbon (uM)")
#plt.savefig(os.path.join(figures_dir, "decay_const_exp_reg_paras_null_condense.png"), dpi = 300, bbox_inches = 'tight')
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
        popt_func, pcov_e = curve_fit(sqf, xdata=X_sub, ydata=y_sub, p0=[0.0001, 1])
        ypred = sqf(X_sub, popt_func[0], popt_func[1])
        yerr = 1-np.sum((ypred-y_sub)**2)/np.sum((y_sub - meany)**2)
        plt.plot(X_sub, ypred, 'r-', alpha = 0.6)
    except:
        popt_func = [np.nan, np.nan]
        yerr = np.nan
    funcres.append([doci,b, popt_func[0], popt_func[1], yerr])
funcresdf = pd.DataFrame.from_records(funcres, columns = ["DOC_initial", "Biomass_initial",  "a", "c", "r2"])
#plt.xticks(ticks = [0.5,1,5,10,50], labels = [0.5,1,5,10,50])
#plt.yticks(ticks = [0.0002, 0.0005,0.001,0.002,0.004], labels = [0.0002, 0.0005,0.001,0.002,0.004])
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
#plt.savefig(os.path.join(figures_dir, "decay_const_null_condense.png"), dpi = 300, bbox_inches = 'tight')

fig, ax = plt.subplots(1,1)
g = sns.scatterplot(data = funcresdf, x = "DOC_initial", y = "a", size = "r2", hue = "c")
leg_handles, leg_labels = g.get_legend_handles_labels()
new_leg = leg_labels
for t in leg_labels:
    if leg_labels.index(t) in [1,2,3,4,5,7,8,9,10,11]:
        new_leg[leg_labels.index(t)] = str(np.round(np.float64(t),3))
    else:
        pass
plt.legend(leg_handles, new_leg, bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
g.set_ylabel("Slope")
g.set_xlabel("Initial available carbon (uM)")

## PREDICT FUNCTION PARAMETERS FOR IMPACT ON DECAY CONSTANT ##
#%%
X = non_100['S_initial']*non_100["activity"]/100#.to_numpy(dtype = float)
y = 1/ non_100['ratio_t_50'].to_numpy(dtype = float)
meany = np.mean(y)
b = non_100["Biomass_initial"].to_numpy(dtype = float)[0]
plt.scatter(X,y, marker = '.', label = "data")
plt.ylabel("Normalized decay constant")
plt.xlabel("Active Shannon diversity")
try:
    popt_func, pcov_e = curve_fit(linres, xdata=X, ydata=y, p0=[0.01,0.1])
    ypred = linres(X, popt_func[0], popt_func[1])
    yerr = 1 - (np.sum((ypred-y)**2)/np.sum((y - meany)**2))
    plt.plot(X, ypred, 'r-', alpha = 0.6, label = 'lin_func')
except:
    popt_func = [np.nan, np.nan]
    yerr = np.nan
plt.text(1.5, 0.3, "slope: "+str(round(popt_func[0],2)))
plt.text(1.5, 0.2, "intercept: "+str(round(popt_func[1],2)))
plt.text(1.5, 0.1, "R2: "+str(round(yerr,2)))
#plt.savefig(os.path.join(figures_dir, "linear_predict_impact_decay_constant_null.png"), dpi = 300)
#%%
non_100["exp_x_var"] = non_100['S_initial']*non_100["activity"]/100#non_100['carbon_species']*non_100["activity"]/100
non_100 = non_100.sort_values(by=["exp_x_var"])
X = non_100['exp_x_var']
y = 1/non_100['ratio_t_50'].to_numpy(dtype = float)
meany = np.mean(y)
#size_var = 1/non_100['ratio-t_50'].to_numpy(dtype = float)
plt.scatter(X,y, marker = '.', label = "data")
plt.ylabel("Normalized decay constant")
plt.xlabel("Active Shannon diversity")
#try:
popt_func, pcov_e = curve_fit(expres, xdata=X, ydata=y, p0=[-0.8,-0.1, 0.9])
ypred = expres(X, popt_func[0], popt_func[1], popt_func[2])
yerr = 1 - (np.sum((ypred-y)**2)/np.sum((y - meany)**2))
plt.plot(X, ypred, 'r-', alpha = 0.6, label = 'exp_func')
plt.text(1.5, 0.3, "p1: "+str(round(popt_func[0],2)))
plt.text(1.5, 0.24, "p2: "+str(round(popt_func[1],2)))
plt.text(1.5, 0.18, "p3: "+str(round(popt_func[2],2)))
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
#size_var = 1/non_100['ratio-t_50'].to_numpy(dtype = float)
plt.scatter(X,y, marker = '.', label = "data")
plt.ylabel("Normalized decay constant")
plt.xlabel("Chaange in active Shannon diversity from\n50% active microbes baseline")
#try:
popt_func, pcov_e = curve_fit(sig_exp, xdata=X, ydata=y, p0=[8,0.5])
ypred = sig_exp(X, popt_func[0], popt_func[1])
yerr = 1 - (np.sum((ypred-y)**2)/np.sum((y - meany)**2))
plt.plot(X, ypred, 'r-', alpha = 0.6, label = 'sig_func')
plt.text(-1.2, 2.5, "p1: "+str(round(popt_func[0],2)))
plt.text(-1.2, 2.3, "p2: "+str(round(popt_func[1],2)))
plt.text(-1.2, 2.1, "R2: "+str(round(yerr,2)))
#plt.savefig(os.path.join(figures_dir, "sigmoid_predict_impact_decay_constant_compet_adapt_B2.png"), dpi = 300, bbox_inches = 'tight')