#%%
import os
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit
from scipy.stats import gamma, beta
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

## LOAD RESULTS

project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
simulations_dir = os.path.join(project_dir, "simulations")
results_dir = os.path.join(project_dir, "results")
figures_dir = os.path.join(project_dir, "figures")

null_prefixes = ["expo_", "gamma_", "1c_", ""]
adapt_prefixes = ["expo_", "gamma_", "1c_", "gen_"]

null_files = []
adapt_files = []
for p,a in zip(null_prefixes, adapt_prefixes):
    filename = os.path.join(results_dir, p + "null_combined_dataset.pkl")
    n_data = pd.read_pickle(filename)
    null_files.append(n_data)
    filename = os.path.join(results_dir, a + "adaptation_combined_dataset.pkl")
    a_data = pd.read_pickle(filename)
    adapt_files.append(a_data)

null_data = pd.concat(null_files)
adapt_data = pd.concat(adapt_files)

null_data = null_data.drop_duplicates()
null_data['DOC_initial_int'] = round(null_data.DOC_initial, -3)
null_data['S_initial_int'] = round(null_data.S_initial, 1)
null_data['ratio_t_50'] = null_data.T_50/null_data.T_50_B1

adapt_data = adapt_data.drop_duplicates()
adapt_data['DOC_initial_int'] = round(adapt_data.DOC_initial, -3)
adapt_data['S_initial_int'] = round(adapt_data.S_initial, 1)
adapt_data['ratio_t_50'] = adapt_data.T_50/adapt_data.T_50_B1

a_n_data = [null_data, adapt_data]
all_data = pd.concat(a_n_data)

#%%
compl = adapt_data.dropna(subset = ['ratio_t_50'])
init_doc_list = list(compl.DOC_initial_int.unique())
activity_list = list(compl.activity.unique())
s_initial_list = list(compl.S_initial_int.unique())

#%%
gammares = []
#for doci in init_doc_list:
for act in activity_list:
    sub = compl[(compl['activity']==act)].sort_values(by=['S_initial']) #(compl['DOC_initial_int']==doci) & 
    #for c_n in [3,6,12,18]:
    #    X = sub[sub.carbon_species == c_n]['S_initial'].to_numpy(dtype = float)
    #    y = sub[sub.carbon_species == c_n]['ratio_t_50'].to_numpy(dtype = float)
    X = sub['S_initial'].to_numpy(dtype = float)
    y = sub['ratio_t_50'].to_numpy(dtype = float)
    b = sub["Biomass_initial"].to_numpy(dtype = float)[0]

    if y.size>0:
        popt_gamma, pcov_e = curve_fit(gamma.pdf, xdata=X, ydata=y, p0=[0.05,1.0])
    else:
        popt_gamma = [np.nan, np.nan, np.nan]
    gammares.append([act, b, popt_gamma[0], popt_gamma[1]])

gammaresdf = pd.DataFrame.from_records(gammares, columns = ["activity", "Biomass_initial", "a", "b"])#"DOC_initial", 

# Check how the plot looks like against the data points for specific cases
resdf = gammaresdf
#for doci in init_doc_list:
for act in activity_list[1:]:
    #for c_n in [3,6,12,18]:
    sub = compl[(compl['activity']==act)].sort_values(by=["S_initial"]) #'DOC_initial_int']==doci) &#& (diversity_data['carbon_species'] == c_n)
    X = sub['S_initial'].to_numpy(dtype = float)
    y = sub['ratio_t_50'].to_numpy(dtype = float)
    b = sub["Biomass_initial"].to_numpy(dtype = float)[0]
    ressub = resdf[(resdf['activity']==act)].reset_index() #(resdf['DOC_initial']==doci) & 
    ypred = gamma.pdf(X, ressub.a[0], loc = ressub.b[0])#, scale  = 1)
    plt.figure()
    plt.title(str(act))# +  " " + str(c_n)) str(doci) + " " + 
    plt.plot(X, ypred, 'r-', alpha = 0.6, label = 'gamma')
    plt.scatter(X,y, marker = '.', label = "data")
    plt.legend()

#%%
x = np.linspace(1.4,4,999)
gamma0505 = gamma.pdf(x, 0.05, 0.05)
gamma00505 = gamma.pdf(x, 0.005, 0.05)
plt.plot(x, gamma0505, label = "0505")
plt.plot(x, gamma00505, label = "00505")
plt.legend()

#%%
gamma0505 = gamma.pdf(x, 0.05, 0.05, scale = 100)
gamma00505 = gamma.pdf(x, 0.005, 0.05, scale = 100)
plt.plot(x, gamma0505, label = "0505")
plt.plot(x, gamma00505, label = "00505")
plt.legend()

#%%
gamma0505 = gamma.pdf(x, 0.05, 2, scale = 10)
gamma00505 = gamma.pdf(x, 0.005, 2, scale = 10)
plt.plot(x, gamma0505, label = "0505")
plt.plot(x, gamma00505, label = "00505")
plt.legend()

#%%
gamma0505 = gamma.pdf(x, 0.05, 1.4, scale = 10)
gamma00505 = gamma.pdf(x, 0.005, 1.4, scale = 10)
plt.plot(x, gamma0505, label = "0505")
plt.plot(x, gamma00505, label = "00505")
plt.legend()

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
        y = sub['ratio_t_50'].to_numpy(dtype = float)
        meany = np.mean(y)
        b = sub["Biomass_initial"].to_numpy(dtype = float)[0]
        axes[idx].scatter(X,y, marker = '.', label = "data")
        try:
            popt_func, pcov_e = curve_fit(abxc, xdata=X, ydata=y, p0=[2,1,10])
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

plt.savefig(os.path.join(figures_dir, "curv_fit_adapt.pdf"), dpi = 300)
plt.savefig(os.path.join(figures_dir, "curv_fit_adapt.png"), dpi = 300)

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
plt.savefig(os.path.join(figures_dir, "adapt_exp_func_char_rxn_tim_ratio.png"), dpi = 300)

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
            y = sub['ratio_t_50'].to_numpy(dtype = float)
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

plt.savefig(os.path.join(figures_dir, "curv_fit_c_n_adapt.pdf"), dpi = 300)
plt.savefig(os.path.join(figures_dir, "curv_fit_c_n_adapt.png"), dpi = 300)

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
plt.savefig(os.path.join(figures_dir, "adapt_exp_func_char_rxn_tim_ratio_c_n.png"), dpi = 300)

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
        y = sub['ratio_t_50'].to_numpy(dtype = float)
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

plt.savefig(os.path.join(figures_dir, "linres_fit_adapt.pdf"), dpi = 300)
plt.savefig(os.path.join(figures_dir, "linres_fit_adapt.png"), dpi = 300)

funcresdf = pd.DataFrame.from_records(funcres, columns = ["DOC_initial", "activity", "Biomass_initial", "a", "c", "r2"])

#%%
resdf=funcresdf[funcresdf.activity<100]
resdf["y_variable"] = resdf["DOC_initial"]
resdf["x_variable"] = resdf["Biomass_initial"]*resdf["activity"]/100
resdf["loga"] = np.log10(abs(resdf.a))
resdf.loc[resdf.a<0, "mark_sign"] = 0
resdf.loc[resdf.a>=0, "mark_sign"] = 1
pos_a = resdf#[resdf.mark_sign>0]

fig, axes = plt.subplots(3,1, figsize = (5,10), sharex = True, sharey = True)
sns.scatterplot("x_variable", "y_variable", alpha = 0.6, size = abs(resdf["a"]), ax = axes[0], data = pos_a, style = "mark_sign")
#sns.scatterplot("x_variable", "y_variable", alpha = 0.6, size = "a", ax = axes[0,0], data = pos_a)
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
plt.savefig(os.path.join(figures_dir, "lin_func_char_rxn_tim_ratio_adapt.png"), dpi = 300)
