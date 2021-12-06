#%%
## Import libraries
import os
import pandas as pd
import seaborn as sns

from DS.solvers.diff_eqn_system import diversity_carbon

## LOAD RESULTS
project_dir = "C:/Users/swami/Documents/Projects/HoliSoils/data"
results_dir = os.path.join(project_dir, "results")
figures_dir = os.path.join(project_dir, "figures")

# Shannon diversity vs carbon stock
data = pd.read_pickle(os.path.join(results_dir, "diversity_data.pkl"))

#%%
sns.scatterplot(x="Shannon", y = "DOC", data = data)