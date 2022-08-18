#%%
## Import libraries
import os
import pandas as pd
import seaborn as sns


## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient", "activity_loss_-02")
results_dir = os.path.join(project_dir, "results")
filestring = "competition_adaptation"#sys.argv[1]
df = pd.read_csv(os.path.join(results_dir, filestring+"_parameters.csv"))
base = df[df.Activity==1]
non_100 = df[df.Activity<1]
#%%
sns.scatterplot(data = base, y = "vmax", x = "Ks", hue = "Oxidation_state", color = "YlGnBu")
#%%
sns.scatterplot(data = non_100, y = "vmax", x = "Ks", hue = "Oxidation_state", color = "YlGnBu")
#%%
sns.scatterplot(data = non_100[non_100['Activity']==0.1], y = "vmax", x = "Ks", hue = "Oxidation_state")
#%%
sns.scatterplot(data = non_100[non_100['Activity']==0.25], y = "vmax", x = "Ks", hue = "Oxidation_state")
#%%
sns.scatterplot(data = non_100[non_100['Activity']==0.5], y = "vmax", x = "Ks", hue = "Oxidation_state")
#%%
sns.scatterplot(data = non_100[non_100['Activity']==0.75], y = "vmax", x = "Ks", hue = "Oxidation_state")
#%%
sns.jointplot(data = non_100[non_100['Activity']==0.75], y = "vmax", x = "Ks", hue = "Oxidation_state")
#%%
sns.displot(data = non_100[non_100['Activity']==0.25], y = "vmax", x = "Oxidation_state", kind = "kde")
#%%
