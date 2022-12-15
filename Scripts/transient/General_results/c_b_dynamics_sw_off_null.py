#%%
import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

## LOAD RESULTS

project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient","gen_spec_skew")
simulations_dir = os.path.join(project_dir, "simulations")
results_dir = os.path.join(project_dir, "results")
figures_dir = os.path.join(project_dir, "figures")

all_data = pd.read_pickle(os.path.join(results_dir, "competition_adaptation_carbon__loss_0.9_c_b_combined_dataset.pkl"))
print(all_data.columns)
all_data['DOC_initial_int'] = round(all_data.DOC_initial, -3)
compl = all_data
all_act = all_data[all_data.activity==100]
sw_off = all_data[all_data.activity<100]
init_doc_list = list(compl.DOC_initial_int.unique())
activity_list = list(compl.activity.unique())
all_data = all_data.drop_duplicates()
all_data["mean_os_ratio"] = all_data.mean_os_end/all_data.mean_os_initial
base = all_data[all_data.activity==100]
#%%
sns.violinplot(data = all_data, x = "biomass_species", y = "DOC_initial_max_os", hue = "carbon_species", palette = "Greens")
#%%
sns.displot(data = all_data, x = "DOC_initial_max_os", kind = "kde", hue = "carbon_species")#, palette = "Greens")
#%%
sns.displot(data = all_data, x = "mean_os_initial", kind = "kde", hue = "carbon_species")#, palette = "Greens")

#%%
sns.violinplot(data = all_data, x = "biomass_species", y = "DOC_end_max_os", hue = "carbon_species", palette = "Greens")
#%%
sns.displot(data = all_data, x = "DOC_end_max_os_pc", kind = "kde", hue = "carbon_species")#, palette = "Greens")
#%%
#%%

#all_data_n["os_ratio_ratio"] = all_data_n.mean_os_ratio/all_data_n.mean_os_ratio_base
#%%
activity_list = [100.0, 75.0, 50.0, 25.0, 10.0]#list(all_data.activity.unique())
all_data["hue_var"] = 1 - all_data.S_in                          itial#*all_data_n.carbon_species#/all_data.vmax_sum_base##all_data.vmax_mean*(all_data.DOC_initial_int/(all_data.carbon_species*all_data.biomass_species))
fig, axes = plt.subplots(5,5, figsize = (12,12), sharex = True, sharey = True)
for i in list(range(5)):
    sub1 = all_data[all_data.DOC_initial_int == init_doc_list[i]]
    for j in list(range(5)):
        sub2 = sub1[sub1.activity == activity_list[j]]
        sub2 = sub2.sort_values(by = ["vmax_ratio"])#, ascending=False)
        #print(i,j)
        #print(sub2.hue_var.describe())
        ax = axes[i,j]#.flatten()[i]
        X = sub2.mean_os_initial#mean_os_ratio_base#mean_os_initial
        y = sub2.mean_os_end
        yline = np.sort(X)
        #if j==0:
        #    ax.scatter(X,y, c = 'g', alpha = 0.6)
        #else:
        ax.scatter(X,y, c = sub2.hue_var, cmap = 'YlGnBu', alpha = 0.6)
        ax.plot(yline, yline, c = 'red')
        ax.set_xticks([-0.4, -0.2, 0, 0.2, 0.4])
        ax.set_yticks([-0.4, -0.2, 0, 0.2, 0.4])
        if i==0:
            ax.set_title(str(activity_list[j])+"%")
        if j ==0:
            ax.set_ylabel(str(init_doc_list[i]))
        if i==4:
            ax.set_xlabel("Weighted Initial OS")
plt.legend()
#%%
#all_data["hue_var"] =all_data.activity/100*all_data.biomass_species*all_data.activity/100
all_data["hue_var"] = np.log(all_data.vmax_ratio*all_data.DOC_initial*all_data.S_initial)#*all_data.vmax_ratio)
all_data["size_var"] = 5#all_data.DOC_initial_int
all_data["y_var"] = all_data.mean_os_end
all_data["x_var"] = all_data.mean_os_initial
yline = np.sort(all_data["x_var"])
f,a = plt.subplots(1,1)
sns.scatterplot(data = all_data, x = "x_var", y = "y_var", ax = a, hue = "hue_var")#, size = "size_var")
a.spines['left'].set_position(('data', 0))
a.spines['bottom'].set_position(('data', 0))
a.spines['top'].set_visible(False)
a.spines['right'].set_visible( False)
a.plot(yline, yline, c = "red")
plt.ylabel("")
plt.xlabel("")
plt.legend(bbox_to_anchor = (0.9,-0.2), ncol = 2)

#%%
base["mean_os_ratio_base"] = base.mean_os_ratio
all_data_n = pd.merge(all_data, base[["Seed", "carbon_species", "biomass_species" ,"DOC_initial_int", "mean_os_ratio_base"]], on = ["Seed", "carbon_species", "biomass_species" ,"DOC_initial_int"])
#%%
high_DOC = all_data[all_data.DOC_initial<9999]
high_DOC["hue_var"] = np.log(1000*high_DOC.vmax_mean)#*high_DOC.S_initial#*all_data.vmax_ratio)
high_DOC["size_var"] = 5#all_data.DOC_initial_int
high_DOC["y_var"] = high_DOC.mean_os_end
high_DOC["x_var"] = high_DOC.mean_os_initial
yline = np.sort(high_DOC["x_var"])
f,a = plt.subplots(1,1)
sns.scatterplot(data = high_DOC, x = "x_var", y = "y_var", ax = a, hue = "hue_var")#, size = "size_var")
a.spines['left'].set_position(('data', 0))
a.spines['bottom'].set_position(('data', 0))
a.spines['top'].set_visible(False)
a.spines['right'].set_visible( False)
a.plot(yline, yline, c = "red")
plt.ylabel("")
plt.xlabel("")
plt.legend(bbox_to_anchor = (0.9,-0.2), ncol = 2)
#%%
oxic_sys = all_data[all_data.mean_os_initial>0]
oxic_sys["hue_var"] = np.log(oxic_sys.DOC_initial_int/(oxic_sys.carbon_species*oxic_sys.biomass_species))#*oxic_sys.activity/100))
oxic_sys["size_var"] = oxic_sys.activity
oxic_sys["y_var"] = oxic_sys.mean_os_end
oxic_sys["x_var"] = oxic_sys.mean_os_initial
yline = np.sort(oxic_sys["x_var"])
f,a = plt.subplots(1,1)
sns.scatterplot(data = oxic_sys, x = "x_var", y = "y_var", ax = a, hue = "hue_var", size = "size_var")
a.spines['left'].set_position(('data', 0))
a.spines['bottom'].set_position(('data', 0))
a.spines['top'].set_visible(False)
a.spines['right'].set_visible( False)
a.plot(yline, yline, c = "red")
plt.ylabel("")
plt.xlabel("")
plt.legend(bbox_to_anchor = (0.9,-0.2), ncol = 2)
#%%
anoxic_sys = all_data[all_data.mean_os_initial<=0]
anoxic_sys["hue_var"] = np.log(anoxic_sys.DOC_initial_int/(anoxic_sys.carbon_species*anoxic_sys.biomass_species))#*oxic_sys.activity/100))
anoxic_sys["size_var"] = anoxic_sys.activity
anoxic_sys["y_var"] = anoxic_sys.mean_os_end
anoxic_sys["x_var"] = anoxic_sys.mean_os_initial
yline = np.sort(anoxic_sys["x_var"])
f,a = plt.subplots(1,1)
sns.scatterplot(data = anoxic_sys, x = "x_var", y = "y_var", ax = a, hue = "hue_var", size = "size_var")
a.spines['left'].set_position(('data', 0))
a.spines['bottom'].set_position(('data', 0))
a.spines['top'].set_visible(False)
a.spines['right'].set_visible( False)
a.plot(yline, yline, c = "red")
plt.ylabel("")
plt.xlabel("")
plt.legend(bbox_to_anchor = (0.9,-0.2), ncol = 2)
#%%

#%%
high_DOC = all_data[all_data.DOC_initial>5000]
high_DOC["hue_var"] = np.log(high_DOC.DOC_initial_int/(high_DOC.carbon_species*high_DOC.biomass_species))#*oxic_sys.activity/100))
high_DOC["size_var"] = high_DOC.activity
high_DOC["y_var"] = high_DOC.mean_os_end
high_DOC["x_var"] = high_DOC.mean_os_initial
yline = np.sort(high_DOC["x_var"])
f,a = plt.subplots(1,1)
sns.scatterplot(data = high_DOC, x = "x_var", y = "y_var", ax = a, hue = "hue_var", size = "size_var")
#a.spines['left'].set_position(('data', 0))
#a.spines['bottom'].set_position(('data', 0))
#a.spines['top'].set_visible(False)
#a.spines['right'].set_visible( False)
a.plot(yline, yline, c = "red")
plt.ylabel("")
plt.xlabel("")
plt.legend(bbox_to_anchor = (0.9,-0.2), ncol = 2)
#%%
#%%
high_DOC = all_data[all_data.activity<50]
high_DOC["hue_var"] = np.log(high_DOC.DOC_initial_int/(high_DOC.carbon_species*high_DOC.biomass_species))#*oxic_sys.activity/100))
high_DOC["size_var"] = high_DOC.activity
high_DOC["y_var"] = high_DOC.mean_os_end
high_DOC["x_var"] = high_DOC.mean_os_initial
yline = np.sort(high_DOC["x_var"])
f,a = plt.subplots(1,1)
sns.scatterplot(data = high_DOC, x = "x_var", y = "y_var", ax = a, hue = "hue_var", size = "size_var")
#a.spines['left'].set_position(('data', 0))
#a.spines['bottom'].set_position(('data', 0))
#a.spines['top'].set_visible(False)
#a.spines['right'].set_visible( False)
a.plot(yline, yline, c = "red")
plt.ylabel("")
plt.xlabel("")
plt.legend(bbox_to_anchor = (0.9,-0.2), ncol = 2)
#%%
activity_list = [100.0, 75.0, 50.0, 25.0, 10.0]#list(all_data.activity.unique())
all_data["hue_var"] = all_data.biomass_species*all_data.carbon_species#/all_data.vmax_sum_base##all_data.vmax_mean*(all_data.DOC_initial_int/(all_data.carbon_species*all_data.biomass_species))
fig, axes = plt.subplots(5,5, figsize = (12,12), sharex = True, sharey = True)
for i in list(range(5)):
    sub1 = all_data[all_data.DOC_initial_int == init_doc_list[i]]
    for j in list(range(5)):
        sub2 = sub1[sub1.activity == activity_list[j]]
        sub2 = sub2.sort_values(by = ["biomass_species", "carbon_species"])#, ascending=False)
        ax = axes[i,j]#.flatten()[i]
        X = sub2.mean_os_initial
        y = sub2.mean_os_end
        yline = np.sort(X)
        ax.scatter(X,y, c = sub2.hue_var, cmap = 'YlGnBu', alpha = 0.6)
        ax.plot(yline, yline, c = 'red')
        ax.set_xticks([-0.4, -0.2, 0, 0.2, 0.4])
        ax.set_yticks([-0.4, -0.2, 0, 0.2, 0.4])
        if i==0:
            ax.set_title(str(activity_list[j])+"%")
        if j ==0:
            ax.set_ylabel(str(init_doc_list[i]))
        if i==4:
            ax.set_xlabel("Weighted Initial OS")

plt.legend()
#%%
## Tracing C/B ratios
print(base.C_B_median.describe())