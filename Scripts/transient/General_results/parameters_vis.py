#%%
## Import libraries
import os
import pandas as pd
import seaborn as sns


## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient", "gen_spec_skew")
results_dir = os.path.join(project_dir, "results")
filestring = "competition_adaptation"#sys.argv[1]
df = pd.read_csv(os.path.join(results_dir, filestring+"_parameters.csv"))
cases = list(df.Sim_series.unique())
for c, a in zip (cases, [100, 10, 10, 10, 10, 10, 25, 25, 25, 25, 25, 50, 50, 50, 50, 50, 75, 75, 75, 75, 75]):        
    df.loc[df['Sim_series']==c, 'activity'] = a
base = df[df.activity==100]
non_100 = df[df.activity<100]
#%%
sns.jointplot(data = base, y = "vmax_mean", x = "k_mean")#, hue = "Oxidation_state", color = "YlGnBu")
#%%
sns.jointplot(data = non_100, y = "vmax_mean", x = "k_mean")#, hue = "Oxidation_state", color = "YlGnBu")
#%%
sns.jointplot(data = non_100[non_100['activity']==10], y = "vmax_mean", x = "k_mean")#, hue = "Oxidation_state")
#%%
sns.jointplot(data = non_100[non_100['activity']==25], y = "vmax_mean", x = "k_mean")#, hue = "Oxidation_state")
#%%
sns.jointplot(data = non_100[non_100['activity']==50], y = "vmax_mean", x = "k_mean")#, hue = "Oxidation_state")
#%%
sns.jointplot(data = non_100[non_100['activity']==75], y = "vmax_mean", x = "k_mean")#, hue = "Oxidation_state")

