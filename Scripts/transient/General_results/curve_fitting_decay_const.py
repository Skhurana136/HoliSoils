#%%
import os
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
simulations_dir = os.path.join(project_dir, "simulations")
results_dir = os.path.join(project_dir, "results")
figures_dir = os.path.join(project_dir, "figures")

filename = os.path.join(results_dir, "null_cont_combined_dataset.pkl")
diversity_data = pd.read_pickle(filename)
diversity_data = diversity_data.drop_duplicates()
diversity_data['DOC_initial_int'] = round(diversity_data.DOC_initial, -3)
diversity_data['decay_const_base'] = 1/diversity_data.T_50_B1
diversity_data['decay_const'] = 1/diversity_data.T_50
diversity_data['decay_ratio'] = diversity_data.decay_const/diversity_data.decay_const_base
diversity_data['S_initial_int'] = round(diversity_data.S_initial, 1)
compl = diversity_data.dropna(subset = ['decay_ratio'])
init_doc_list = list(compl.DOC_initial_int.unique())
activity_list = list(compl.activity.unique())
s_initial_list = list(compl.S_initial_int.unique())

#%%
def abxc(x, a, b, c):
    return a * np.exp(-b * x) + c

activity_list_plot = [10.0, 25.0, 50.0, 75.0]
funcres = []
fig, ax = plt.subplots(5,4, figsize = (10,10), sharex = True, sharey = "row")
axes = ax.flatten()
idx = 0
for doci in init_doc_list:
    for act in activity_list_plot:
        sub = compl[(compl['DOC_initial_int']==doci) & (compl['activity']==act)].sort_values(by=['S_initial'])
        X = sub['S_initial'].to_numpy(dtype = float)
        y = sub['decay_ratio'].to_numpy(dtype = float)
        meany = np.mean(y)
        b = sub["Biomass_initial"].to_numpy(dtype = float)[0]
        axes[idx].scatter(X,y, marker = '.', label = "data")
        try:
            popt_func, pcov_e = curve_fit(abxc, xdata=X, ydata=y, p0=[1,0.5,0.2])
            ypred = abxc(X, popt_func[0], popt_func[1], popt_func[2])
            yerr = 1 - (np.sum((ypred-y)**2)/np.sum((y - meany)**2))
            axes[idx].plot(X, ypred, 'r-', alpha = 0.6, label = 'exp_func')
        except:
            popt_func = [np.nan, np.nan, np.nan]
            yerr = np.nan
        funcres.append([doci, act, b, popt_func[0], popt_func[1], popt_func[2], yerr])
        idx=idx+1

for a,act in zip(ax[0,:], activity_list_plot):
    a.set_title("Activity:\n"+ str(int(act))+"%")
for a,doci in zip (ax[:,0], init_doc_list):
    a.set_ylabel("Initial DOC:\n"+str(int(doci))+ " (uM)")
for a in ax[-1,:]:
    a.set_xlabel("Initial\nShannon diversity")

plt.savefig(os.path.join(figures_dir, "decay_curv_fit.pdf"), dpi = 300)
plt.savefig(os.path.join(figures_dir, "decay_curv_fit.png"), dpi = 300)

funcresdf = pd.DataFrame.from_records(funcres, columns = ["DOC_initial", "activity", "Biomass_initial", "a", "b", "c", "r2"])

#%%
resdf=funcresdf[funcresdf.activity<100]
resdf["y_variable"] = resdf["DOC_initial"]
resdf["x_variable"] = resdf["Biomass_initial"]*resdf["activity"]/100
resdf["loga"] = np.log10(abs(resdf.a))
resdf.loc[resdf.a<0, "mark_sign"] = 0
resdf.loc[resdf.a>=0, "mark_sign"] = 1
pos_a = resdf#[resdf.mark_sign>0]

fig, axes = plt.subplots(2,2, figsize = (6,6), sharex = True, sharey = True)
sns.scatterplot("x_variable", "y_variable", alpha = 0.6, size = "loga", ax = axes[0,0], data = pos_a, style = "mark_sign")
#sns.scatterplot("x_variable", "y_variable", alpha = 0.6, size = "a", ax = axes[0,0], data = pos_a)
axes[0,0].set_ylabel ("Initial DOC")
axes[0,0].legend(title = "a")
sns.scatterplot("x_variable", "y_variable", alpha = 0.6, size= "b",ax = axes[0,1], data = pos_a)
axes[0,1].legend(title = "b")
sns.scatterplot("x_variable", "y_variable", alpha = 0.6, size= "c",ax = axes[1,0], data = pos_a)
axes[1,0].legend(title = "c")
axes[1,0].set_ylabel ("Initial DOC")
axes[1,0].set_xlabel ("Active biomass")
sns.scatterplot("x_variable", "y_variable", alpha = 0.6, size= "r2",ax = axes[1,1], data = pos_a)
axes[1,1].legend(title = "R2")
axes[1,1].set_xlabel ("Active biomass")
fig.tight_layout()
plt.savefig(os.path.join(figures_dir, "null_exp_func_decay.png"), dpi = 300)

#%%
def abxc(x, a, b, c):
    return a * np.exp(-b * x) + c

activity_list_plot = [10.0, 25.0, 50.0, 75.0]
funcres = []
fig, ax = plt.subplots(5,4, figsize = (10,10), sharex = True, sharey = "row")
axes = ax.flatten()
idx = 0
for doci in init_doc_list:
    for act in activity_list_plot:
        for c_n in [3,6,12,18]:
            sub = compl[(compl['DOC_initial_int']==doci) & (compl['activity']==act) & (compl['carbon_species']==c_n)].sort_values(by=['S_initial'])
            X = sub['S_initial'].to_numpy(dtype = float)
            y = sub['decay_ratio'].to_numpy(dtype = float)
            meany = np.mean(y)
            if y.size>0:
                b = sub["Biomass_initial"].to_numpy(dtype = float)[0]
                axes[idx].scatter(X,y, marker = '.', label = c_n)
                try:
                    popt_func, pcov_e = curve_fit(abxc, xdata=X, ydata=y, p0=[2,1,10])
                    ypred = abxc(X, popt_func[0], popt_func[1], popt_func[2])
                    yerr = 1 - (np.sum((ypred-y)**2)/np.sum((y - meany)**2))
                    axes[idx].plot(X, ypred, alpha = 0.6, label = c_n)
                except:
                    popt_func = [np.nan, np.nan, np.nan]
                    yerr = np.nan
                funcres.append([doci, act, b, c_n, popt_func[0], popt_func[1], popt_func[2], yerr])
        idx=idx+1

for a,act in zip(ax[0,:], activity_list_plot):
    a.set_title("Activity:\n"+ str(int(act))+"%")
for a,doci in zip (ax[:,0], init_doc_list):
    a.set_ylabel("Initial DOC:\n"+str(int(doci))+ " (uM)")
for a in ax[-1,:]:
    a.set_xlabel("Initial\nShannon diversity")

plt.savefig(os.path.join(figures_dir, "curv_fit_c_n_null_decay.pdf"), dpi = 300)
plt.savefig(os.path.join(figures_dir, "curv_fit_c_n_null_decay.png"), dpi = 300)

funcresdf = pd.DataFrame.from_records(funcres, columns = ["DOC_initial", "activity", "Biomass_initial", "carbon_species", "a", "b", "c", "r2"])

#%%
resdf=funcresdf[funcresdf.activity<100].sort_values(by = ["carbon_species"], ascending=False)
resdf["y_variable"] = resdf["DOC_initial"]
resdf["x_variable"] = resdf["Biomass_initial"]*resdf["activity"]/100
resdf["loga"] = np.log10(abs(resdf.a))
resdf.loc[resdf.a<0, "mark_sign"] = 0
resdf.loc[resdf.a>=0, "mark_sign"] = 1
pos_a = resdf#[resdf.mark_sign>0]

fig, axes = plt.subplots(2,2, figsize = (8,8), sharex = True, sharey = True)
sns.scatterplot(x = "x_variable", y = "y_variable", alpha = 0.7, size = "loga", hue = "carbon_species", ax = axes[0,0], data = pos_a, style = "mark_sign")
#sns.scatterplot(x = "x_variable", y = "y_variable", alpha = 0.6, size = "a", hue = "carbon_species", ax = axes[0,0], data = pos_a)
axes[0,0].set_ylabel ("Initial DOC")
axes[0,0].legend(title = "a")
sns.scatterplot(x = "x_variable", y = "y_variable",  alpha = 0.7, size= "b",hue = "carbon_species", ax = axes[0,1], data = pos_a)
axes[0,1].legend(title = "b")
sns.scatterplot(x = "x_variable", y = "y_variable", alpha = 0.7, size= "c",hue = "carbon_species", ax = axes[1,0], data = pos_a)
axes[1,0].legend(title = "c")
axes[1,0].set_ylabel ("Initial DOC")
axes[1,0].set_xlabel ("Active biomass")
sns.scatterplot(x = "x_variable", y = "y_variable",  alpha = 0.7, size= "r2",hue = "carbon_species", ax = axes[1,1], data = pos_a)
axes[1,1].legend(title = "R2")
axes[1,1].set_xlabel ("Active biomass")
fig.tight_layout()
plt.savefig(os.path.join(figures_dir, "null_exp_func_decay_ratio_c_n.png"), dpi = 300)

#%%
def linres(x, a, c):
    return a * x + c

activity_list_plot = [10.0, 25.0, 50.0, 75.0]
funcres = []
fig, ax = plt.subplots(5,4, figsize = (10,10), sharex = True, sharey = "row")
axes = ax.flatten()
idx = 0
for doci in init_doc_list:
    for act in activity_list_plot:
        sub = compl[(compl['DOC_initial_int']==doci) & (compl['activity']==act)].sort_values(by=['S_initial'])
        X = sub['S_initial'].to_numpy(dtype = float)
        y = sub['decay_ratio'].to_numpy(dtype = float)
        meany = np.mean(y)
        b = sub["Biomass_initial"].to_numpy(dtype = float)[0]
        axes[idx].scatter(X,y, marker = '.', label = "data")
        if y.size>0:
            try:
                popt_func, pcov_e = curve_fit(linres, xdata=X, ydata=y, p0=[-2,10])
                ypred = linres(X, popt_func[0], popt_func[1])
                yerr = 1 - (np.sum((ypred-y)**2)/np.sum((y - meany)**2))
                axes[idx].plot(X, ypred, 'r-', alpha = 0.6, label = 'exp_func')
            except:
                popt_func = [np.nan, np.nan]
                yerr = np.nan
        else:
            popt_func = [np.nan, np.nan]
            yerr = np.nan
        funcres.append([doci, act, b, popt_func[0], popt_func[1],yerr])
        idx=idx+1

for a,act in zip(ax[0,:], activity_list_plot):
    a.set_title("Activity:\n"+ str(int(act))+"%")
for a,doci in zip (ax[:,0], init_doc_list):
    a.set_ylabel("Initial DOC:\n"+str(int(doci))+ " (uM)")
for a in ax[-1,:]:
    a.set_xlabel("Initial\nShannon diversity")

plt.savefig(os.path.join(figures_dir, "linres_fit_null_decay.pdf"), dpi = 300)
plt.savefig(os.path.join(figures_dir, "linres_fit_null_decay.png"), dpi = 300)

funcresdf = pd.DataFrame.from_records(funcres, columns = ["DOC_initial", "activity", "Biomass_initial", "a", "c", "r2"])

#%%
resdf=funcresdf[funcresdf.activity<100]
resdf["y_variable"] = resdf["DOC_initial"]
resdf["x_variable"] = resdf["Biomass_initial"]*resdf["activity"]/100
resdf["loga"] = np.log10(abs(resdf.a))
resdf.loc[resdf.a>0, "mark_sign"] = 1
resdf.loc[resdf.a<=0, "mark_sign"] = 0
pos_a = resdf#[resdf.mark_sign>0]

fig, axes = plt.subplots(3,1, figsize = (5,10), sharex = True, sharey = True)
sns.scatterplot("x_variable", "y_variable", alpha = 0.6, size = abs(resdf["a"]), ax = axes[0], data = pos_a, style = "mark_sign")
#sns.scatterplot("x_variable", "y_variable", alpha = 0.6, size = "a", ax = axes[0], data = pos_a)
axes[0].set_ylabel ("Initial DOC")
axes[0].legend(title = "a")
sns.scatterplot("x_variable", "y_variable", alpha = 0.6, size= "c",ax = axes[1], data = pos_a)
axes[1].legend(title = "c")
axes[1].set_ylabel ("Initial DOC")
sns.scatterplot("x_variable", "y_variable", alpha = 0.6, size= "r2",ax = axes[2], data = pos_a)
axes[2].legend(title = "R2")
axes[2].set_ylabel ("Initial DOC")
axes[2].set_xlabel ("Active biomass")
fig.tight_layout()
plt.savefig(os.path.join(figures_dir, "lin_func_decay_ratio_null.png"), dpi = 300)

#%%
resdf=funcresdf[funcresdf.activity<100]
resdf["y_variable"] = resdf["DOC_initial"]
resdf["x_variable"] = resdf["Biomass_initial"]*resdf["activity"]/100
resdf["loga"] = np.log10(abs(resdf.a))
resdf.loc[resdf.a>0, "mark_sign"] = 1
resdf.loc[resdf.a<=0, "mark_sign"] = 0
markstyle = ['o', '+']
lablist = ["negative", "positive"]

fig,ax = plt.subplots()

for m in [0.0, 1.0]:
    subset = resdf[resdf.mark_sign == m]
    xpl = subset["x_variable"]
    ypl = subset["y_variable"]
    spl = abs(subset["a"])
    ax.scatter(xpl, ypl, alpha = 0.8, s = 2000*spl, c = "steelblue", marker = markstyle[int(m)], label = lablist[int(m)])

plt.ylabel ("Initial C [N $L^{-3}$]")
plt.xlabel ("Active biomass [N $L^{-3}$]")
legend1 = ax.legend(title = "Slope sign", loc = "lower right")
ax.add_artist(legend1)

sizescatter = ax.scatter(xpl, ypl, s = 2000*spl, alpha = 0.01)

handles, labels = sizescatter.legend_elements(prop="sizes", alpha = 0.5)
labels2 = [40/2000, 60/2000, 80/2000, 100/2000, 120/2000]
legend2 = ax.legend(handles[::2], labels2, loc = "upper right", title = "Slope magnitude")
plt.savefig(os.path.join(figures_dir, "null_lin_func_decay_ratio_slope.png"), dpi = 300)
#%%
plt.figure()
sizescatter = plt.scatter(resdf["x_variable"], resdf["y_variable"], alpha = 0.8, s= 200*resdf["c"])
handles, labels = sizescatter.legend_elements(prop="sizes", c = "steelblue", alpha=0.8)
labels2 = [20/200, 40/200, 60/200, 80/200, 100/200, 120/200, 140/200]
plt.legend(handles[::2], labels2, loc="lower right", title="Intercept")
plt.ylabel ("Initial C [N $L^{-3}$]")
plt.xlabel ("Active biomass [N $L^{-3}$]")
plt.savefig(os.path.join(figures_dir, "null_lin_func_decay_ratio_intercept.png"), dpi = 300)
#%%
plt.figure()
sizescatter = plt.scatter(resdf["x_variable"], resdf["y_variable"], alpha = 0.8, s= 200*resdf["r2"])
handles, labels = sizescatter.legend_elements(prop="sizes", c = "steelblue", alpha=0.8)
labels2 = [20/200, 40/200, 60/200, 80/200, 100/200, 120/200, 140/200]
plt.legend(handles[::2], labels2, loc="lower right", title="R2")
plt.ylabel ("Initial C [N $L^{-3}$]")
plt.xlabel ("Active biomass [N $L^{-3}$]")
plt.savefig(os.path.join(figures_dir, "null_lin_func_decay_ratio_r2.png"), dpi = 300)