#%%
import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

import matplotlib.lines as mlines

from scipy.optimize import curve_fit

## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient")
all_data = pd.read_csv(os.path.join(project_dir,"simulation_results_diversity_parameters.csv"))
print(all_data.shape)
print(all_data.dtypes)
#%%
compl = all_data.dropna(subset = ['ratio_t_50'])
act_full = all_data[all_data.Sim_series=="b_1_all_"]
act_off = all_data[all_data.Sim_series!="b_1_all_"]
extr_full = act_full[(act_full['DOC_initial_int']==2000.)|(act_full['DOC_initial_int']==10000.)]
extr_off = act_off[(act_off['DOC_initial_int']==2000.)|(act_off['DOC_initial_int']==10000.)]
compl = compl.sort_values(by = 'active_H_c_connections')
compl_act_off = compl[compl.activity<100]
#%%
plt.figure(figsize=(8,4))
h = sns.scatterplot(x = "Variance", y = "FD_initial", size = "carbon_species", hue = "biomass_species", data = act_full)
h.set_xlabel("Log normal variance\n(x mean)", fontsize = 16)
h.set_ylabel("Initial functional\ndiversity: Variance", fontsize = 16)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.yscale("log")
plt.legend(bbox_to_anchor=(1,1.1), fontsize = 14)
plt.tight_layout()
plt.savefig(os.path.join(project_dir, "imposed_func_div.png"), dpi = 300)
#%%
## OR
sns.boxplot('biomass_species', 'FD_initial', hue = 'Variance',data = act_full)
plt.ylabel("Initial functional diversity: Variance")
plt.yscale("log")

#%%
fig, axes = plt.subplots(3,1, figsize = (8,10), sharex = True)
sns.scatterplot(x = "FD_initial", y = "FD_ratio", hue = 'DOC_initial_int', style = "Variance", data = extr_full, ax = axes[0])
axes[0].set_ylabel("Functional diversity:\nPeak/Initial", fontsize = 16)
sns.scatterplot(x = "FD_initial", y = "Biomass_ratio", hue = 'DOC_initial_int', style = "Variance", data = extr_full, ax = axes[1])
axes[1].set_ylabel("Gain in biomass:\nPeak/Initial", fontsize = 16)
sns.scatterplot(x = "FD_initial", y = "Decay_constant", hue = 'DOC_initial_int', style = "Variance", data = extr_full, ax = axes[2])
axes[2].set_xscale("log")
axes[2].set_xlabel("Initial functional diversity: Variance", fontsize = 16)
axes[2].set_ylabel("Decay constant (1/day)", fontsize = 16)
for a in axes.flatten():
    a.tick_params(axis='both', which='major', labelsize=14)
    a.legend().set_visible(False)
axes[0].legend(bbox_to_anchor=(1,1), fontsize = 14)
plt.tight_layout()
plt.savefig(os.path.join(project_dir, "ecosystem_func_div_extr.png"), dpi = 300)
#%%
fig, axes = plt.subplots(3,1, figsize = (8,14), sharex = True)
sns.scatterplot(x = "FD_initial", y = "FD_ratio", hue = 'DOC_initial_int', style = "Variance", data = act_full, ax = axes[0])
axes[0].set_ylabel("Functional diversity:\nPeak/Initial", fontsize = 16)
sns.scatterplot(x = "FD_initial", y = "Biomass_ratio", hue = 'DOC_initial_int', style = "Variance", data = act_full, ax = axes[1])
axes[1].set_ylabel("Gain in biomass:\nPeak/Initial", fontsize = 16)
sns.scatterplot(x = "FD_initial", y = "Decay_constant", hue = 'DOC_initial_int', style = "Variance", data = act_full, ax = axes[2])
axes[2].set_xscale("log")
axes[2].set_xlabel("Initial functional diversity: Variance", fontsize = 16)
axes[2].set_ylabel("Decay constant (1/day)", fontsize = 16)
for a in axes.flatten():
    a.tick_params(axis='both', which='major', labelsize=14)
    a.legend().set_visible(False)
axes[2].legend(bbox_to_anchor=(0.75,-0.2), fontsize = 14, ncol = 2)
plt.tight_layout()
plt.savefig(os.path.join(project_dir, "ecosystem_func_div_all.png"), dpi = 300)
#%%
sns.scatterplot(x = "S_initial", y = "FD_initial", hue = 'carbon_species', style = "Variance", data = act_full)
plt.yscale("log")
plt.xlabel("Initial Shannon diversity index")
plt.ylabel("Initial functional diversity: Vatiance")
plt.legend(bbox_to_anchor=(1, 1))
#%%
sns.scatterplot(x = "Biomass_ratio", y = "FD_ratio", hue = 'DOC_initial_int', style = "Variance", size = "carbon_species",data = extr_full)
plt.axhline(y=1, c = 'red', label = "Same community")
plt.axvline(x=1, c = 'black', label = 'Dying community')
plt.yscale("log")
plt.xlabel("Ratio: Biomass")
plt.ylabel("Ratio: Functional diversity")
plt.legend(bbox_to_anchor=(1, 1))
#%%
### LOSS IN ACTIVITY
sns.scatterplot(x = "FD_initial", y = "Biomass_ratio", hue = 'DOC_initial', style = "Variance", data = extr_off)
plt.xscale("log")
plt.xlabel("Initial functional diversity: Variance")
plt.ylabel("Gain in biomass: Peak/Initial")
plt.legend(bbox_to_anchor=(1, 1))
#%%
g = sns.FacetGrid(extr_off, col="DOC_initial_int", row="activity",hue = 'carbon_species', palette = "Blues")
g.map(sns.scatterplot,'FD_initial','Biomass_ratio')
g.set(xscale="log")
g.set(xlabel="Initial functional diversity\n(Variance)")
g.set(ylabel="Biomass ratio: Peak/initial")
g.add_legend()
#%%
sns.scatterplot(x = "Biomass_ratio", y = "FD_ratio", hue = 'activity', size = "Variance", style = "DOC_initial_int",data = extr_off)
plt.axhline(y=1, c = 'red', label = "Same community")
plt.axvline(x=1, c = 'black', label = 'Dying community')
plt.xlabel("Ratio: Biomass")
plt.ylabel("Ratio: Functional diversity (Variance)")
plt.legend(bbox_to_anchor=(1, 1))
#%%
g = sns.FacetGrid(act_off, row="Variance", col="activity",hue = 'DOC_initial_int', palette = "Blues")
g.map(sns.scatterplot,'Biomass_ratio','FD_ratio')
g.set(xlabel="Biomass ratio: Peak/initial")
g.set(ylabel="Functional diversity\nratio: Peak/initial")
g.add_legend()
## Community threshold for gain in biomass and functional diversity in times of disturbance
#%%
extr_off['log_fd_initial'] = np.log(extr_off.FD_initial)
sns.scatterplot(x = "Biomass_ratio", y = "FD_ratio", hue = 'log_fd_initial', size = "DOC_initial", alpha = 0.5,data = extr_off)
plt.axhline(y=1, c = 'red', label = "Same community")
plt.axvline(x=1, c = 'black', label = 'Dying community')
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Ratio: Biomass")
plt.ylabel("Ratio: Functional diversity (Variance)")
plt.legend(bbox_to_anchor=(1, 1))
#%%
## Predicting decay constant
g = sns.FacetGrid(all_data, row="Variance", col="activity",hue = 'DOC_initial_int', palette = "magma")
g.map(sns.scatterplot,'active_H_c_connections','Decay_constant')
#%%
## Decay constant as a function of initial diversity?
## Fully active communities first
#%%
sns.scatterplot('FD_initial', 'Decay_constant', hue = 'DOC_initial_int', style = 'Variance', data = extr_full)
plt.xscale("log")
plt.ylabel("Decay constant (/day)")
plt.xlabel("Initial functional diversity: Variance")
plt.legend(bbox_to_anchor=(1,1))
#%%
sns.scatterplot('FD_cov', 'Decay_constant', hue = 'DOC_initial_int', style = 'Variance', data = extr_full)
plt.xscale("log")
plt.ylabel("Decay constant (/day)")
plt.xlabel("Initial functional diversity: CoV")
plt.legend(bbox_to_anchor=(1,1))
#%%
#%%
sns.scatterplot('vmax_mean', 'Decay_constant', hue = 'DOC_initial_int', style = 'Variance', data = act_full)
plt.xscale("log")
plt.ylabel("Decay constant (/day)")
plt.xlabel("Mean vmax")
plt.legend(bbox_to_anchor=(1,1))
#%%
sns.scatterplot('DOC_initial', 'Decay_constant', hue = 'vmax_mean', style = 'Variance', data = act_full)
#plt.xscale("log")
plt.ylabel("Decay constant (/day)")
plt.xlabel("Initial carbon available")
plt.legend(bbox_to_anchor=(1,1))
#%%
act_full["xocxcarbon"] = act_full.FD_cov*act_full.DOC_initial
sns.scatterplot('xocxcarbon', 'Decay_constant', hue = 'DOC_initial_int', style = 'Variance', data = act_full)
plt.xscale("log")
plt.ylabel("Decay constant (/day)")
plt.xlabel("Coeff of variation in functional capacity of community x carbon available")
plt.legend(bbox_to_anchor=(1,1))
#%%
sns.scatterplot('FD_initial', 'Decay_constant', hue = 'DOC_initial_int', style = 'Variance', data = extr_off)
plt.xscale("log")
plt.ylabel("Decay constant (/day)")
plt.xlabel("Initial functional diversity: Variance")
plt.legend(bbox_to_anchor=(1,1))
#%%
#%%
sns.scatterplot('FD_cov', 'Decay_constant', hue = 'DOC_initial_int', style = 'Variance', data = extr_off)
plt.xscale("log")
plt.ylabel("Decay constant (/day)")
plt.xlabel("Initial functional diversity: CoV")
plt.legend(bbox_to_anchor=(1,1))
#%%
sns.scatterplot('FD_initial','Norm_b1_decay_const', data = extr_off, style = "Variance", alpha = 0.5, hue = "DOC_initial_int")
plt.xscale("log")
plt.yscale("log")
plt.ylabel("Normalized decay constant (-)")
plt.xlabel("Initial functional diversity:Variance")
plt.legend(bbox_to_anchor = (1,1))
#%%
sns.scatterplot('FD_cov','Norm_b1_decay_const', data = extr_off, style = "Variance", alpha = 0.5, hue = "DOC_initial_int")
plt.xscale("log")
plt.yscale("log")
plt.ylabel("Normalized decay constant (-)")
plt.xlabel("Initial functional diversity: CoV")
plt.legend(bbox_to_anchor = (1,1))
#%%
g = sns.FacetGrid(extr_off, row="Variance", col="activity",hue = 'DOC_initial_int', palette = "magma")
g.map(sns.scatterplot,'FD_cov','Norm_b1_decay_const')
g.set(xscale="log")
g.add_legend()
#%%
sns.scatterplot('FD_initial_ratio','Norm_b1_decay_const', data = extr_off, style = "Variance", alpha = 0.5, hue = "DOC_initial_int")
plt.axhline(y=1, c = 'red', label = 'No change in\ndecay constant')
plt.axvline(x=1, c = 'black', label = 'No change in\nfunctional diversity')
plt.xscale("log")
plt.ylabel("Normalized decay constant (-)")
plt.xlabel("Normalized functional diversity (Variance)")
plt.legend(bbox_to_anchor = (1,1))
#%%
g = sns.FacetGrid(extr_off, col="Variance",col_wrap=2,hue = 'DOC_initial_int', palette = "Blues")
g.map(sns.scatterplot,'FD_initial_ratio','Norm_b1_decay_const')
#g.set(yscale="log")
g.set(xscale="log")
g.add_legend()
#%%
g = sns.FacetGrid(extr_off, col="Variance",col_wrap=2,hue = 'DOC_initial_int', palette = "Blues")
g.map(sns.scatterplot,'FD_initial_ratio','Norm_b1_decay_const')
g.set(xscale="log")
g.add_legend()
#%%
sns.scatterplot('active_H_c_connections','Norm_b1_decay_const', data = compl, style = "Variance", alpha = 0.5, hue = "DOC_initial_int")
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

def logistic_func(x, L, b):
    return L/(1+np.exp(-b*(x)))
#%%
sns.scatterplot('active_H_c_connections','Norm_b1_decay_const', data = compl_act_off, style = "Variance", alpha = 0.5, hue = "DOC_initial_int")
#plt.ylim(bottom=1)
plt.yscale("log")
plt.legend(bbox_to_anchor = (1,1))
for doc in [1000.,2000.,5000.,10000.,15000.]:
    X = compl_act_off[compl_act_off.DOC_initial_int==doc]["active_H_c_connections"]
    y = compl_act_off[compl_act_off.DOC_initial_int==doc]['Norm_b1_decay_const']
    meany = np.mean(y)
    popt_func, pcov_e = curve_fit(powlaw, xdata=X, ydata=y, p0=[0.4,0.7])
    ypred = powlaw(X, popt_func[0], popt_func[1])
    yerr = 1 - (np.sum((ypred-y)**2)/np.sum((y - meany)**2))
    plt.plot(X, ypred, 'r-', alpha = 0.6, label = doc)
    print(yerr)
plt.legend(bbox_to_anchor=(1,1))

#%%
## PREDICT FUNCTION PARAMETERS FOR IMPACT ON DECAY CONSTANT ##
#%%
X = compl_act_off['active_H']
y =compl_act_off['Norm_b1_decay_const'].to_numpy(dtype = float)
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
#plt.savefig(os.path.join(project_dir, "exponential_predict_impact_decay_constant_compet_adapt.png"), dpi = 300)
#%%
plt.scatter(data = all_data, x="H_diff_from_50_act", y="Normalized_decay_constant_b2",marker = '.', alpha = 0.1)
plt.yscale("log")
plt.ylabel("Normalized decay constant")
plt.xlabel("Change in active Shannon diversity")
#plt.savefig(os.path.join(figures_dir, "impact_decay_constant_compet_adapt_S_active_B2.png"), dpi = 300, bbox_inches='tight')

### PREDICT IMPACT OF CHANGING DIVERSITY on DECAY CONSTANT B2
#%%
all_data = all_data.sort_values(by=["H_diff_from_50_act"])
compl = all_data.dropna(subset = ['Normalized_decay_constant_b2'])
#%%
fig, axes = plt.subplots(3,2, figsize=(8,10), sharex = True)
ax = axes.flatten()
for v, a in zip([0.01, 0.1,0.5,1.,1.5],[0,1,2,3,4]):
    X = compl[compl.Variance==v]['H_diff_from_50_act']
    y = compl[compl.Variance==v]["Normalized_decay_constant_b2"]
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
#plt.savefig(os.path.join(project_dir, "logistic_func_predict_impact_decay_constant_compet_adapt_B2.png"), dpi = 300, bbox_inches = 'tight')
## the response of all communities seem to be similar
##(coefficients are similar with varying performance
## (R2))
#%%
X = compl['H_diff_from_50_act']
y = compl["Normalized_decay_constant_b2"]
mean = np.mean(y)
sns.scatterplot(x="H_diff_from_50_act", y = "Normalized_decay_constant_b2", alpha = 0.5, hue = "DOC_initial_int", style="Variance", data = compl)
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
sns.scatterplot(x="vmax_ratio", y = "Norm_b1_decay_const", hue = "activity", size = "S_initial_int", style = "Variance", data = extr_off)
plt.xscale("log")
plt.yscale("log")
plt.legend(bbox_to_anchor=(1,1))
#%%
sns.scatterplot(x="vmax_ratio", y = "Norm_b1_decay_const", hue = "DOC_initial_int",style = "activity", data = extr_off)
plt.xscale("log")
plt.yscale("log")
plt.legend(bbox_to_anchor=(1,1))
#%%
sns.scatterplot(x="vmax_ratio", y = "Norm_b1_decay_const", hue = "Variance",style = "activity", data = extr_off)
plt.xscale("log")
plt.yscale("log")
plt.legend(bbox_to_anchor=(1,1))
#%%
sns.scatterplot(y="FD_initial_ratio", x = "vmax_ratio", hue = "activity", style = "Variance", data = act_off)
plt.yscale("log")
plt.xscale("log")
plt.xlabel("Normalized vmax")
plt.ylabel("Normalized functional diversity")
#%%
g = sns.FacetGrid(act_off, col="Variance", row = "activity")
g.map(sns.scatterplot,'FD_initial_ratio','vmax_ratio')
g.set(yscale="log")
g.set(xscale="log")
g.add_legend()
#%%
sns.scatterplot(y="FD_initial", x = "vmax_median", hue = "activity", style = "Variance",data =act_off)
plt.yscale("log")
plt.xscale("log")
#%%
sns.scatterplot(y="FD_ratio", x = "FD_cov", hue = "DOC_initial_int", size = "activity", style = "Variance",data =act_off)
plt.xscale("log")
plt.legend(bbox_to_anchor=(1,1))
#%%
#%%
sns.scatterplot(y="Biomass_ratio", x = "FD_cov", hue = "DOC_initial_int", size = "activity_x", style = "Variance",data =act_off)
plt.xscale("log")
plt.legend(bbox_to_anchor=(1,1))
#%%
sns.scatterplot(y="Norm_b1_decay_const", x = "FD_cov", hue = "DOC_initial_int", size = "activity", style = "Variance",data =act_off)
plt.yscale("log")
plt.xscale("log")
plt.legend(bbox_to_anchor=(1,1))
#%%
g = sns.FacetGrid(act_off, col="activity", row="Variance",hue = 'DOC_initial_int', palette = "Blues")
g.map(sns.scatterplot,'FD_cov','Norm_b1_decay_const')
g.set(yscale="log")
g.set(xscale="log")
g.add_legend()

#%%
faster = all_data[all_data.Normalized_decay_constant>1]
low_variance_faster = faster[faster.Variance<0.5]
high_variance_faster = faster[faster.Variance>1.0]
mid_variance_faster = faster[(faster.Variance>=0.5) & (faster.Variance<=1.0)]
#%%
fig, axes = plt.subplots(3,1,figsize =(4,10))#sharex = True,
sns.scatterplot('vmax_ratio','Normalized_decay_constant', data = low_variance_faster, size = "activity_x", alpha = 0.5, hue = "DOC_initial_int_x", ax = axes[0])
sns.scatterplot('vmax_ratio','Normalized_decay_constant', data = mid_variance_faster, size = "activity_x", alpha = 0.5, hue = "DOC_initial_int_x", ax = axes[1])
sns.scatterplot('vmax_ratio','Normalized_decay_constant', data = high_variance_faster, size = "activity_x", alpha = 0.5, hue = "DOC_initial_int_x", ax = axes[2])
axes[0].set_title ("Similar communities")
axes[1].set_title ("Transitioning communites")
axes[2].set_title ("Dissimilar communities")
#plt.ylim(bottom=1)
#plt.yscale("log")
for a in axes.flatten():
    a.legend(bbox_to_anchor = (1,1))

#%%
faster_b2 = all_data[all_data.Normalized_decay_constant_b2>1]
low_variance_faster = faster[faster.Variance<0.5]
high_variance_faster = faster_b2[faster_b2.Variance>1.0]
mid_variance_faster = faster_b2[(faster_b2.Variance>=0.5) & (faster_b2.Variance<=1.)]

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
act_full['NOSC_ratio'] = act_full.NOSC_maxbio/act_full.NOSC_initial
sns.scatterplot(x = "NOSC_initial", y = "NOSC_maxbio", size = "S_initial_int", alpha = 0.5,hue = 'DOC_initial_int', style ="Variance", data = act_full)
plt.xlabel("Initial NOSC")
plt.ylabel("NOSC at peak biomass")
plt.legend(bbox_to_anchor=(1, 1))
#%%
sns.scatterplot(x = "NOSC_initial", y = "NOSC_maxbio", size = "S_initial_int",alpha = 0.5, hue ="Variance", data = act_full[act_full.DOC_initial_int==15000])
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