#%%
import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns
import scipy

import statsmodels.api as sm
import statsmodels.formula.api as smf

## LOAD RESULTS
project_dir = "C:/Users/swkh9804/Documents/Scripts/HoliSoils/Data"#os.path.join("D:/", "Projects", "HoliSoils","data","transient")
all_data = pd.read_csv(os.path.join(project_dir,"simulation_results_temporal_initial_conditions_decay_const_with_paras.csv"))
print(all_data.shape)
print(all_data.dtypes)
#%%
all_data['C_avail_x'] = all_data['DOC_initial_int']/1000.
all_data['C_avail'] = all_data['DOC_initial_int']/all_data['k_mean']
k10 = all_data[all_data['%C']==100.]
act_full = k10[k10.Sim_series=="b_1_all_"]
doc_full = act_full[act_full.C_pool=='DOC']
toc_full = act_full[act_full.C_pool=='TOC']
rec_full = act_full[act_full.C_pool=='reduced_C']
oxc_full = act_full[act_full.C_pool=='oxidized_C']
#no_na = all_data.dropna()
#%%
para_features = ['Variance','carbon_species', 'biomass_species', 'activity', 'active_H_c_connections', 'FD_initial', 'DOC_initial_int','NOSC_initial','vmax_mean', 'm_mean', 'k_mean', 'exo_enz_mean']
all_features = ['Variance','carbon_species', 'biomass_species', 'activity', 'active_H_c_connections', 'FD_initial', 'DOC_initial_int','NOSC_initial','vmax_mean', 'm_mean', 'k_mean', 'exo_enz_mean','decay_const', 'FD_ratio', 'Biomass_ratio']
para_corr = doc_full[para_features].corr()
para_mask = np.triu(np.ones_like(para_corr, dtype=bool))
sns.heatmap(para_corr, vmin=-1, vmax = 1, center = 0, mask = para_mask, cmap = 'vlag')
#%%
doc_corr = doc_full[all_features].corr()
doc_mask = np.triu(np.ones_like(doc_corr, dtype=bool))
toc_corr = toc_full[all_features].corr()
toc_mask = np.triu(np.ones_like(toc_corr, dtype=bool))
rec_corr = rec_full[all_features].corr()
rec_mask = np.triu(np.ones_like(rec_corr, dtype=bool))
oxc_corr = oxc_full[all_features].corr()
oxc_mask = np.triu(np.ones_like(oxc_corr, dtype=bool))
fig, axes = plt.subplots(nrows=2, ncols=2, sharex = True, sharey = True, figsize = (10,8))
sns.heatmap(toc_corr, vmin = -1, vmax = 1, center = 0, mask = toc_mask, cmap='vlag', ax = axes[0][0])
sns.heatmap(doc_corr, vmin = -1, vmax = 1, center = 0, mask = doc_mask, cmap='vlag', ax = axes[0][1])
sns.heatmap(rec_corr, vmin = -1, vmax = 1, center = 0, mask = oxc_mask, cmap='vlag', ax = axes[1][0])
sns.heatmap(oxc_corr, vmin = -1, vmax = 1, center = 0, mask = rec_mask, cmap='vlag', ax = axes[1][1])
axes[0][0].set_title("TOC")
axes[0][1].set_title("DOC")
axes[1][0].set_title("Oxidized_C")
axes[1][1].set_title("Reduced_C")
#%%
features_to_test = ['Variance','carbon_species', 'biomass_species', 'activity', 'active_H_c_connections', 'FD_initial', 'DOC_initial_int']
feature_data = doc_full[features_to_test]
print(feature_data.shape)
#%%
corr = feature_data.corr()
mask = np.triu(np.ones_like(corr, dtype=bool))
sns.heatmap(corr, vmin = -1, vmax = 1, center = 0, mask = mask, cmap='vlag')
#Limited linear relationships
#%%
def run_anova_test_biomass_species(df, cpool_test,col_to_test):
    df_sub = df[df.C_pool==cpool_test]
    bio_1_4 = df_sub[df_sub.biomass_species==4][col_to_test].values
    bio_2_1 = df_sub[df_sub.biomass_species==8][col_to_test].values
    bio_2_5 = df_sub[df_sub.biomass_species==12][col_to_test].values
    bio_2_8 = df_sub[df_sub.biomass_species==16][col_to_test].values
    bio_3_0 = df_sub[df_sub.biomass_species==20][col_to_test].values
    bio_3_2 = df_sub[df_sub.biomass_species==24][col_to_test].values
    bio_3_3 = df_sub[df_sub.biomass_species==28][col_to_test].values
    bio_3_5 = df_sub[df_sub.biomass_species==32][col_to_test].values
    #Test for assumptions
    #Test for normality
    for dataset in [bio_1_4, bio_2_1, bio_2_5,bio_2_8,bio_3_0, bio_3_2, bio_3_3, bio_3_5]:
        print(scipy.stats.shapiro(dataset))
    #Test for variance
    print(scipy.stats.bartlett(bio_1_4, bio_2_1, bio_2_5,bio_2_8,bio_3_0, bio_3_2, bio_3_3, bio_3_5))
    #Test for variation
    kruskal = scipy.stats.kruskal(bio_1_4, bio_2_1, bio_2_5,bio_2_8,bio_3_0, bio_3_2, bio_3_3, bio_3_5, nan_policy = 'omit')
    print(kruskal)
    

#%%
#Run for decay constant of DOC governed by biomass species number
run_anova_test_biomass_species(act_full, 'DOC', 'decay_const')
#ShapiroResult(statistic=0.6732597947120667, pvalue=2.2119217400972837e-38)
#ShapiroResult(statistic=0.7106255292892456, pvalue=1.0307665414478878e-36)
#ShapiroResult(statistic=0.7325906753540039, pvalue=1.1936939937663602e-35)
#ShapiroResult(statistic=0.7100738286972046, pvalue=9.711703266950476e-37)
#ShapiroResult(statistic=0.7641540765762329, pvalue=5.440868356632077e-34)
#ShapiroResult(statistic=0.7748679518699646, pvalue=2.183646427887468e-33)
#ShapiroResult(statistic=0.7468857765197754, pvalue=6.420260288406586e-35)
#ShapiroResult(statistic=0.761171817779541, pvalue=3.7289562146774917e-34)
#The datasets are not normal because p-value is less than 0.01.
#BartlettResult(statistic=59.11793812774482, pvalue=2.2638821640309805e-10)
#p-value is less than 0.01, so we can reject the null hypothesis and say that not all the groups have the same variance.
#Assumptions to run factorial ANOVA are not met
#KruskalResult(statistic=10.79211166272051, pvalue=0.14794757951324053)
#pvalue is greater than 0.05, so the decay constant is not different for each microbial group
#%%
#Run for decay constant of TOC governed by biomass species number
run_anova_test_biomass_species(act_full, 'TOC', 'decay_const')
#ShapiroResult(statistic=0.669887900352478, pvalue=1.592018828695432e-38)
#ShapiroResult(statistic=0.7107539176940918, pvalue=1.0451753873021139e-36)
#ShapiroResult(statistic=0.7302412986755371, pvalue=9.116510421068601e-36)
#ShapiroResult(statistic=0.702407956123352, pvalue=4.2845716207858225e-37)
#ShapiroResult(statistic=0.758575439453125, pvalue=2.6919029559584143e-34)
#ShapiroResult(statistic=0.7716890573501587, pvalue=1.4380325652127989e-33)
#ShapiroResult(statistic=0.7423242330551147, pvalue=3.7228706781273973e-35)
#ShapiroResult(statistic=0.7569100856781006, pvalue=2.1873516228154584e-34)
#The datasets are not normal because p-value is less than 0.01.
#BartlettResult(statistic=28.315072313585084, pvalue=0.00019277915340305105)
#p-value is less than 0.01, so we can reject the null hypothesis and say that not all the groups have the same variance.
#KruskalResult(statistic=4.468785148671871, pvalue=0.7244727371662647)
#pvalue is greater than 0.05, so the decay constant is not different for each microbial group
#%%
#Run for decay constant of reduced C governed by biomass species number
run_anova_test_biomass_species(act_full, 'reduced_C', 'decay_const')
#ShapiroResult(statistic=nan, pvalue=1.0)
#ShapiroResult(statistic=0.7228778600692749, pvalue=3.9642528181446143e-36)
#ShapiroResult(statistic=0.7255913019180298, pvalue=5.376806870467295e-36)
#ShapiroResult(statistic=0.6947033405303955, pvalue=1.9146853443629944e-37)
#ShapiroResult(statistic=0.7681847810745239, pvalue=9.122155749907581e-34)
#ShapiroResult(statistic=0.768735945224762, pvalue=9.79569013569975e-34)
#ShapiroResult(statistic=0.7536840438842773, pvalue=1.4679280727424723e-34)
#ShapiroResult(statistic=0.7702281475067139, pvalue=1.188688907568283e-33)
#The datasets are not normal because p-value is less than 0.01.
#BartlettResult(statistic=nan, pvalue=nan)
#No interpretation
#KruskalResult(statistic=11.255232048047572, pvalue=0.12786690900690803)
#pvalue is greater than 0.05, so the decay constant is not different for each microbial group
#%%
#Run for decay constant of oxidized C governed by biomass species number
run_anova_test_biomass_species(act_full, 'oxidized_C', 'decay_const')
#ShapiroResult(statistic=nan, pvalue=1.0)
#ShapiroResult(statistic=nan, pvalue=1.0)
#ShapiroResult(statistic=nan, pvalue=1.0)
#ShapiroResult(statistic=nan, pvalue=1.0)
#ShapiroResult(statistic=nan, pvalue=1.0)
#ShapiroResult(statistic=nan, pvalue=1.0)
#ShapiroResult(statistic=nan, pvalue=1.0)
#ShapiroResult(statistic=nan, pvalue=1.0)
#BartlettResult(statistic=nan, pvalue=nan)
#KruskalResult(statistic=17.003638890830814, pvalue=0.017372733361234155)
#%%
#Run for initial FD of governed by biomass species number
run_anova_test_biomass_species(act_full, 'DOC', 'FD_initial')
#ShapiroResult(statistic=0.2230616807937622, pvalue=0.0)
#ShapiroResult(statistic=0.36607033014297485, pvalue=0.0)
#ShapiroResult(statistic=0.3079666495323181, pvalue=0.0)
#ShapiroResult(statistic=0.27479273080825806, pvalue=0.0)
#ShapiroResult(statistic=0.49384385347366333, pvalue=9.80908925027372e-45)
#ShapiroResult(statistic=0.31869786977767944, pvalue=0.0)
#ShapiroResult(statistic=0.3972741365432739, pvalue=0.0)
#ShapiroResult(statistic=0.372001051902771, pvalue=0.0)
#The datasets are not normal because p-value is less than 0.01.
#BartlettResult(statistic=698.6882923161544, pvalue=1.322612533302575e-146)
#p-value is less than 0.01, so we can reject the null hypothesis and say that not all the groups have the same variance.
#Assumptions to run factorial ANOVA are not met
#KruskalResult(statistic=30.5478994717306, pvalue=7.531621981869381e-05)
#pvalue is less than 0.05, so the initial functional diversity is different for some of the microbial groups
#statistically different distributions
model = smf.ols('FD_initial ~ C(biomass_species)', data=doc_full).fit()
sm.stats.anova_lm(model, typ=2)
#%%
#Run for gain in biomass of governed by biomass species number
run_anova_test_biomass_species(act_full, 'DOC', 'Biomass_ratio')
#ShapiroResult(statistic=0.9840590953826904, pvalue=2.4595943060035097e-08)
#ShapiroResult(statistic=0.9428347945213318, pvalue=3.880958553350782e-18)
#ShapiroResult(statistic=0.9753326177597046, pvalue=3.272739371174005e-11)
#ShapiroResult(statistic=0.9739757776260376, pvalue=1.3437987701758747e-11)
#ShapiroResult(statistic=0.9514113068580627, pvalue=1.1796871378382976e-16)
#ShapiroResult(statistic=0.9558745622634888, pvalue=8.283763054371955e-16)
#ShapiroResult(statistic=0.9788958430290222, pvalue=3.9720876587878706e-10)
#ShapiroResult(statistic=0.9722157716751099, pvalue=4.4238085417092066e-12)
#The datasets are not normal because p-value is less than 0.01.
#BartlettResult(statistic=4710.571972951083, pvalue=0.0)
#p-value is less than 0.01, so we can reject the null hypothesis and say that not all the groups have the same variance.
#Assumptions to run factorial ANOVA are not met
#KruskalResult(statistic=277.48493185411644, pvalue=3.860994978463975e-56)
#pvalue is less than 0.05, so the gain in biomass is different for some of the microbial groups
#statistically different distributions
#%%
def run_anova_test_variance(df, cpool_test,col_to_test):
    df_sub = df[df.C_pool==cpool_test]
    var_0_01 = df_sub[df_sub.Variance==0.01][col_to_test].values
    var_0_1 = df_sub[df_sub.Variance==0.1][col_to_test].values
    var_0_5 = df_sub[df_sub.Variance==0.5][col_to_test].values
    var_1_0 = df_sub[df_sub.Variance==1.][col_to_test].values
    var_1_5 = df_sub[df_sub.Variance==1.5][col_to_test].values
    #Test for assumptions
    #Test for normality
    for dataset in [var_0_01, var_0_1, var_0_5, var_1_0, var_1_5]:
        print(scipy.stats.shapiro(dataset))
    #Test for variance
    print(scipy.stats.bartlett(var_0_01, var_0_1, var_0_5, var_1_0, var_1_5))
    #Test for variation
    kruskal = scipy.stats.kruskal(var_0_01, var_0_1, var_0_5, var_1_0, var_1_5, nan_policy = 'omit')
    print(kruskal)

#%%
#Run for community evolution of governed by community type
run_anova_test_variance(act_full, 'DOC', 'FD_ratio')
#ShapiroResult(statistic=0.7545769810676575, pvalue=7.47312471024425e-42)
#ShapiroResult(statistic=0.989211916923523, pvalue=7.892042575008418e-09)
#ShapiroResult(statistic=0.936612606048584, pvalue=2.7043895213576823e-24)
#ShapiroResult(statistic=0.8900383710861206, pvalue=7.743664486183106e-31)
#ShapiroResult(statistic=0.904102623462677, pvalue=3.904050784708043e-29)
#The datasets are not normal because p-value is less than 0.01.
#BartlettResult(statistic=18387.35096372093, pvalue=0.0)
#p-value is less than 0.01, so we can reject the null hypothesis and say that not all the groups have the same variance.
#ANOVA conditions are not met
#KruskalResult(statistic=5061.740927581108, pvalue=0.0)
#pvalue is less than 0.05, so the community evolution is different for some of the community types
#statistically different distributions
#%%
#Run for gain in biomass of governed by community type
run_anova_test_variance(act_full, 'DOC', 'Biomass_ratio')
#ShapiroResult(statistic=0.809991717338562, pvalue=3.762289073888357e-38)
#ShapiroResult(statistic=0.8097764253616333, pvalue=3.6262801659797616e-38)
#ShapiroResult(statistic=0.8347630500793457, pvalue=3.317383672222509e-36)
#ShapiroResult(statistic=0.8736051321029663, pvalue=1.2701367805532485e-32)
#ShapiroResult(statistic=0.8618928790092468, pvalue=8.721082587536786e-34)
#The datasets are not normal because p-value is less than 0.01.
#BartlettResult(statistic=21.48082286153952, pvalue=0.00025420371832127696)
#p-value is less than 0.01, so we can reject the null hypothesis and say that not all the groups have the same variance.
#ANOVA conditions are not met
#KruskalResult(statistic=1010.9629129772984, pvalue=1.502298701881689e-217)
#pvalue is less than 0.05, so the gain in biomass is different for some of the community types
#statistically different distributions
#%%
#Run for decay constant of carbon pool governed by community type
run_anova_test_variance(act_full, 'TOC', 'decay_const')
#ShapiroResult(statistic=0.8764711022377014, pvalue=2.520263414646037e-32)
#ShapiroResult(statistic=0.8781024813652039, pvalue=3.7435953896958005e-32)
#ShapiroResult(statistic=0.8878821134567261, pvalue=4.39772796204774e-31)
#ShapiroResult(statistic=0.8766567707061768, pvalue=2.63577454252394e-32)
#ShapiroResult(statistic=0.8700128197669983, pvalue=5.4737856882039554e-33)
#The datasets are not normal because p-value is less than 0.01.
#BartlettResult(statistic=4306.440607601319, pvalue=0.0)
#p-value is less than 0.01, so we can reject the null hypothesis and say that not all the groups have the same variance.
#ANOVA conditions are not met
#KruskalResult(statistic=1569.2956960092843, pvalue=0.0)
#pvalue is less than 0.05, so the gain in biomass is different for some of the community types
#statistically different distributions
#%%
#Run for decay constant of carbon pool governed by community type
run_anova_test_variance(act_full, 'DOC', 'decay_const')
#ShapiroResult(statistic=0.8738792538642883, pvalue=1.3554610665652412e-32)
#ShapiroResult(statistic=0.8757646083831787, pvalue=2.126098999917009e-32)
#ShapiroResult(statistic=0.8854638934135437, pvalue=2.355179809202997e-31)
#ShapiroResult(statistic=0.8733806610107422, pvalue=1.20440328339206e-32)
#ShapiroResult(statistic=0.8696289658546448, pvalue=5.008382495458457e-33)
#The datasets are not normal because p-value is less than 0.01.
#BartlettResult(statistic=4150.333195869615, pvalue=0.0)
#p-value is less than 0.01, so we can reject the null hypothesis and say that not all the groups have the same variance.
#ANOVA conditions are not met
#KruskalResult(statistic=1528.3276620689335, pvalue=0.0)
#pvalue is less than 0.05, so the gain in biomass is different for some of the community types
#statistically different distributions
#%%
#Run for decay constant of carbon pool governed by community type
run_anova_test_variance(act_full, 'reduced_C', 'decay_const')
#ShapiroResult(statistic=0.873421311378479, pvalue=1.216059264969623e-32)
#ShapiroResult(statistic=0.8751832842826843, pvalue=1.8495072362836388e-32)
#ShapiroResult(statistic=nan, pvalue=1.0)
#ShapiroResult(statistic=0.8724944591522217, pvalue=9.771352627191929e-33)
#ShapiroResult(statistic=0.8666689395904541, pvalue=2.5417051295548062e-33)
#Most of the datasets are not normal because p-value is less than 0.01.
#BartlettResult(statistic=nan, pvalue=nan)
#NAN
#ANOVA conditions are not met
#KruskalResult(statistic=1465.4582056701286, pvalue=0.0)
#pvalue is less than 0.05, so the gain in biomass is different for some of the community types
#statistically different distributions
#%%
#Run for decay constant of carbon pool governed by community type
run_anova_test_variance(act_full, 'oxidized_C', 'decay_const')
#ShapiroResult(statistic=nan, pvalue=1.0)
#ShapiroResult(statistic=nan, pvalue=1.0)
#ShapiroResult(statistic=nan, pvalue=1.0)
#ShapiroResult(statistic=nan, pvalue=1.0)
#ShapiroResult(statistic=nan, pvalue=1.0)
#BartlettResult(statistic=nan, pvalue=nan)
#NAN
#ANOVA conditions are not met
#KruskalResult(statistic=1305.939141962005, pvalue=1.7158369398394523e-281)
#pvalue is less than 0.05, so the gain in biomass is different for some of the community types
#statistically different distributions

#%%
def run_anova_test_c_avail(df, cpool_test,col_to_test):
    df_sub = df[df.C_pool==cpool_test]
    doc_1k = df_sub[df_sub.DOC_initial_int==1000][col_to_test].values
    doc_2k = df_sub[df_sub.DOC_initial_int==2000][col_to_test].values
    doc_5k = df_sub[df_sub.DOC_initial_int==5000][col_to_test].values
    doc_10k = df_sub[df_sub.DOC_initial_int==10000][col_to_test].values
    doc_15k = df_sub[df_sub.DOC_initial_int==15000][col_to_test].values
    #Test for assumptions
    #Test for normality
    for dataset in [doc_1k, doc_2k,doc_5k,doc_10k, doc_15k]:
        print(scipy.stats.shapiro(dataset))
    #Test for variance
    print(scipy.stats.bartlett(doc_1k, doc_2k,doc_5k,doc_10k, doc_15k))
    #Test for variation
    kruskal = scipy.stats.kruskal(doc_1k, doc_2k,doc_5k,doc_10k, doc_15k, nan_policy = 'omit')
    print(kruskal)

#%%
#Run for community evolution of governed by C availability
run_anova_test_c_avail(act_full, 'DOC', 'FD_ratio')
#ShapiroResult(statistic=0.7728613615036011, pvalue=1.028679189676205e-40)
#ShapiroResult(statistic=0.7754626870155334, pvalue=1.5146354841194083e-40)
#ShapiroResult(statistic=0.7792455554008484, pvalue=2.6761297422443194e-40)
#ShapiroResult(statistic=0.778817892074585, pvalue=2.5083803030799955e-40)
#ShapiroResult(statistic=0.7783592939376831, pvalue=2.34036461720745e-40)
#The datasets are not normal because p-value is less than 0.01.
#BartlettResult(statistic=129.20727713486076, pvalue=5.7533949334551644e-27)
#p-value is less than 0.01, so we can reject the null hypothesis and say that not all the groups have the same variance.
#ANOVA conditions are not met
#KruskalResult(statistic=28.11851543220837, pvalue=1.1801725794842952e-05)
#pvalue is less than 0.05, so the community evolution is different for some of the community types
#statistically different distributions
#%%
#Run for gain in biomass of governed by C availability
run_anova_test_variance(act_full, 'DOC', 'Biomass_ratio')
#ShapiroResult(statistic=0.809991717338562, pvalue=3.762289073888357e-38)
#ShapiroResult(statistic=0.8097764253616333, pvalue=3.6262801659797616e-38)
#ShapiroResult(statistic=0.8347630500793457, pvalue=3.317383672222509e-36)
#ShapiroResult(statistic=0.8736051321029663, pvalue=1.2701367805532485e-32)
#ShapiroResult(statistic=0.8618928790092468, pvalue=8.721082587536786e-34)
#The datasets are not normal because p-value is less than 0.01.
#BartlettResult(statistic=21.48082286153952, pvalue=0.00025420371832127696)
#p-value is less than 0.01, so we can reject the null hypothesis and say that not all the groups have the same variance.
#ANOVA conditions are not met
#KruskalResult(statistic=1010.9629129772984, pvalue=1.502298701881689e-217)
#pvalue is less than 0.05, so the gain in biomass is different for some of the community types
#statistically different distributions
#%%
#Run for decay constant of carbon pool governed by C availability
run_anova_test_variance(act_full, 'TOC', 'decay_const')
#ShapiroResult(statistic=0.8764711022377014, pvalue=2.520263414646037e-32)
#ShapiroResult(statistic=0.8781024813652039, pvalue=3.7435953896958005e-32)
#ShapiroResult(statistic=0.8878821134567261, pvalue=4.39772796204774e-31)
#ShapiroResult(statistic=0.8766567707061768, pvalue=2.63577454252394e-32)
#ShapiroResult(statistic=0.8700128197669983, pvalue=5.4737856882039554e-33)
#The datasets are not normal because p-value is less than 0.01.
#BartlettResult(statistic=4306.440607601319, pvalue=0.0)
#p-value is less than 0.01, so we can reject the null hypothesis and say that not all the groups have the same variance.
#ANOVA conditions are not met
#KruskalResult(statistic=1569.2956960092843, pvalue=0.0)
#pvalue is less than 0.05, so the gain in biomass is different for some of the community types
#statistically different distributions
#%%
#Run for decay constant of carbon pool governed by C availability
run_anova_test_variance(act_full, 'DOC', 'decay_const')
#ShapiroResult(statistic=0.8738792538642883, pvalue=1.3554610665652412e-32)
#ShapiroResult(statistic=0.8757646083831787, pvalue=2.126098999917009e-32)
#ShapiroResult(statistic=0.8854638934135437, pvalue=2.355179809202997e-31)
#ShapiroResult(statistic=0.8733806610107422, pvalue=1.20440328339206e-32)
#ShapiroResult(statistic=0.8696289658546448, pvalue=5.008382495458457e-33)
#The datasets are not normal because p-value is less than 0.01.
#BartlettResult(statistic=4150.333195869615, pvalue=0.0)
#p-value is less than 0.01, so we can reject the null hypothesis and say that not all the groups have the same variance.
#ANOVA conditions are not met
#KruskalResult(statistic=1528.3276620689335, pvalue=0.0)=
#pvalue is less than 0.05, so the gain in biomass is different for some of the community types
#statistically different distributions
#%%
#Run for decay constant of carbon pool governed by C availability
run_anova_test_variance(act_full, 'reduced_C', 'decay_const')
#ShapiroResult(statistic=0.873421311378479, pvalue=1.216059264969623e-32)
#ShapiroResult(statistic=0.8751832842826843, pvalue=1.8495072362836388e-32)
#ShapiroResult(statistic=nan, pvalue=1.0)
#ShapiroResult(statistic=0.8724944591522217, pvalue=9.771352627191929e-33)
#ShapiroResult(statistic=0.8666689395904541, pvalue=2.5417051295548062e-33)
#Most of the datasets are not normal because p-value is less than 0.01.
#BartlettResult(statistic=nan, pvalue=nan)
#NAN
#ANOVA conditions are not met
#KruskalResult(statistic=1465.4582056701286, pvalue=0.0)
#pvalue is less than 0.05, so the gain in biomass is different for some of the community types
#statistically different distributions
#%%
#Run for decay constant of carbon pool governed by community type
run_anova_test_variance(act_full, 'oxidized_C', 'decay_const')
#ShapiroResult(statistic=nan, pvalue=1.0)
#ShapiroResult(statistic=nan, pvalue=1.0)
#ShapiroResult(statistic=nan, pvalue=1.0)
#ShapiroResult(statistic=nan, pvalue=1.0)
#ShapiroResult(statistic=nan, pvalue=1.0)
#BartlettResult(statistic=nan, pvalue=nan)
#NAN
#ANOVA conditions are not met
#KruskalResult(statistic=1305.939141962005, pvalue=1.7158369398394523e-281)
#pvalue is less than 0.05, so the gain in biomass is different for some of the community types
#statistically different distributions

#%%
subset = act_full[act_full.biomass_species <12]
g = sns.displot(data= toc_full, x = "decay_const", hue = "DOC_initial_int", row = "carbon_species", col = "biomass_species", kind = 'kde')
g.set(xscale="log")
#%%
subset = toc_full[toc_full.biomass_species <8]
g = sns.jointplot(data= subset, x = "FD_initial", y = "decay_const", hue = "DOC_initial_int")
g.set(xscale="log")
#%%
plt.figure()
sm.qqplot(np.log10(toc_full['FD_initial']), line = 's')
plt.title("Initial functional diversity")
plt.figure()
sm.qqplot(np.log10(toc_full['decay_const']), line = 's')
plt.title("Decay constant")
#%%
#Log transformation of the dataframe
no_inf = act_full.drop(act_full[act_full.decay_const==np.inf].index)
lognona = no_inf.apply(lambda x: np.log10(x) if np.issubdtype(x.dtype, np.number) else x)
#%%
### Tests for normality ###
## Initial functional diversity
print("Shapiro test: ", scipy.stats.shapiro(lognona['FD_initial']))
print("KS test: ",scipy.stats.kstest(lognona['FD_initial'], scipy.stats.norm.cdf))
#%%
# Decay constant
print("Shapiro test: ", scipy.stats.shapiro(lognona['decay_const']))
print("KS test: ",scipy.stats.kstest(lognona['decay_const'], scipy.stats.norm.cdf))
### ORDINARY LINEAR REGRESSION ###
#%%
#Run prediction for decay constant of different pools of carbon:
def dec_cont_f_fd_doc (df, cpool_test):
    df_sub  = df[df.C_pool==cpool_test]
    model = smf.ols('decay_const ~ FD_initial + DOC_initial_int', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results

def dec_cont_f_cbvdoc (df, cpool_test):
    df_sub  = df[df.C_pool==cpool_test]
    model = smf.ols('decay_const ~ carbon_species + biomass_species + Variance + DOC_initial_int', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results

def dec_cont_f_fcbvdoc (df, cpool_test):
    df_sub  = df[df.C_pool==cpool_test]
    model = smf.ols('decay_const ~ FD_initial + carbon_species + biomass_species + Variance + DOC_initial_int', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results

def dec_stop_f_fd_doc (df, cpool_test):
    df_sub  = df[df.C_pool==cpool_test]
    model = smf.ols('decay_stop ~ FD_initial + DOC_initial_int', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results

def dec_stop_f_fcbvdoc (df, cpool_test):
    df_sub  = df[df.C_pool==cpool_test]
    model = smf.ols('decay_stop ~ FD_initial + carbon_species + biomass_species + Variance + DOC_initial_int', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results


def dec_stop_f_cbvdoc (df, cpool_test):
    df_sub  = df[df.C_pool==cpool_test]
    model = smf.ols('decay_stop ~ carbon_species + biomass_species + Variance + DOC_initial_int', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results

def dec_cont_f_cbvosdoc (df, cpool_test):
    df_sub  = df[df.C_pool==cpool_test]
    model = smf.ols('decay_const ~ carbon_species + biomass_species + Variance + NOSC_initial + DOC_initial_int', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results

def dec_stop_f_cbvosdoc (df, cpool_test):
    df_sub  = df[df.C_pool==cpool_test]
    model = smf.ols('decay_stop ~ carbon_species + biomass_species + Variance + NOSC_initial + DOC_initial_int', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results

#%%
dec_stop_toc_fd_doc_m,dec_stop_toc_fd_doc_r = dec_stop_f_fd_doc(lognona, 'TOC')
dec_stop_doc_fd_doc_m,dec_stop_doc_fd_doc_r = dec_stop_f_fd_doc(lognona, 'DOC')
dec_stop_rec_fd_doc_m,dec_stop_rec_fd_doc_r = dec_stop_f_fd_doc(lognona, 'reduced_C')
dec_stop_oxc_fd_doc_m,dec_stop_oxc_fd_doc_r = dec_stop_f_fd_doc(lognona, 'oxidized_C')
dec_stop_toc_fcbvdoc_m,dec_stop_toc_fcbvdoc_r = dec_stop_f_fcbvdoc(lognona, 'TOC')
dec_stop_doc_fcbvdoc_m,dec_stop_doc_fcbvdoc_r = dec_stop_f_fcbvdoc(lognona, 'DOC')
dec_stop_rec_fcbvdoc_m,dec_stop_rec_fcbvdoc_r = dec_stop_f_fcbvdoc(lognona, 'reduced_C')
dec_stop_oxc_fcbvdoc_m,dec_stop_oxc_fcbvdoc_r = dec_stop_f_fcbvdoc(lognona, 'oxidized_C')
dec_stop_toc_cbvdoc_m,dec_stop_toc_cbvdoc_r = dec_stop_f_cbvdoc(lognona, 'TOC')
dec_stop_doc_cbvdoc_m,dec_stop_doc_cbvdoc_r = dec_stop_f_cbvdoc(lognona, 'DOC')
dec_stop_rec_cbvdoc_m,dec_stop_rec_cbvdoc_r = dec_stop_f_cbvdoc(lognona, 'reduced_C')
dec_stop_oxc_cbvdoc_m,dec_stop_oxc_cbvdoc_r = dec_stop_f_cbvdoc(lognona, 'oxidized_C')
dec_toc_cbvdoc_m,dec_toc_cbvdoc_r = dec_cont_f_cbvdoc(lognona, 'TOC')
dec_doc_cbvdoc_m,dec_doc_cbvdoc_r = dec_cont_f_cbvdoc(lognona, 'DOC')
dec_rec_cbvdoc_m,dec_rec_cbvdoc_r = dec_cont_f_cbvdoc(lognona, 'reduced_C')
dec_oxc_cbvdoc_m,dec_oxc_cbvdoc_r = dec_cont_f_cbvdoc(lognona, 'oxidized_C')
dec_toc_cbvosdoc_m,dec_toc_cbvosdoc_r = dec_cont_f_cbvosdoc(lognona, 'TOC')
dec_doc_cbvosdoc_m,dec_doc_cbvosdoc_r = dec_cont_f_cbvosdoc(lognona, 'DOC')
dec_rec_cbvosdoc_m,dec_rec_cbvosdoc_r = dec_cont_f_cbvosdoc(lognona, 'reduced_C')
dec_oxc_cbvosdoc_m,dec_oxc_cbvosdoc_r = dec_cont_f_cbvosdoc(lognona, 'oxidized_C')
dec_stop_toc_cbvosdoc_m,dec_stop_toc_cbvosdoc_r = dec_stop_f_cbvosdoc(lognona, 'TOC')
dec_stop_doc_cbvosdoc_m,dec_stop_doc_cbvosdoc_r = dec_stop_f_cbvosdoc(lognona, 'DOC')
dec_stop_rec_cbvosdoc_m,dec_stop_rec_cbvosdoc_r = dec_stop_f_cbvosdoc(lognona, 'reduced_C')
dec_stop_oxc_cbvosdoc_m,dec_stop_oxc_cbvosdoc_r = dec_stop_f_cbvosdoc(lognona, 'oxidized_C')

#%%
#Run prediction for decay constant TOC
dec_toc_fd_doc_m,dec_toc_fd_doc_r = dec_cont_f_fd_doc(lognona, 'TOC')
#                            OLS Regression Results                            
#==============================================================================
#Dep. Variable:            decay_const   R-squared:                       0.751
#Model:                            OLS   Adj. R-squared:                  0.751
#Method:                 Least Squares   F-statistic:                 1.086e+04
#Date:                Thu, 19 Jan 2023   Prob (F-statistic):               0.00
#Time:                        17:09:37   Log-Likelihood:                 1877.4
#No. Observations:                7200   AIC:                            -3749.
#Df Residuals:                    7197   BIC:                            -3728.
#Df Model:                           2                                         
#Covariance Type:            nonrobust                                         
#===================================================================================
#                      coef    std err          t      P>|t|      [0.025      0.975]
#-----------------------------------------------------------------------------------
#Intercept          -4.4912      0.020   -228.546      0.000      -4.530      -4.453
#FD_initial          0.1038      0.001     88.792      0.000       0.102       0.106
#DOC_initial_int     0.5946      0.005    117.611      0.000       0.585       0.604
#==============================================================================
#Omnibus:                      625.707   Durbin-Watson:                   0.383
#Prob(Omnibus):                  0.000   Jarque-Bera (JB):             1041.059
#Skew:                          -0.639   Prob(JB):                    8.65e-227
#Kurtosis:                       4.355   Cond. No.                         64.2
#==============================================================================
#
#Notes:
#[1] Standard Errors assume that the covariance matrix of the errors is correctly specified.
#%%
dec_doc_fd_doc_m,dec_doc_fd_doc_r = dec_cont_f_fd_doc(lognona, 'DOC')
#                            OLS Regression Results                            
#==============================================================================
#Dep. Variable:            decay_const   R-squared:                       0.755
#Model:                            OLS   Adj. R-squared:                  0.755
#Method:                 Least Squares   F-statistic:                 1.107e+04
#Date:                Thu, 19 Jan 2023   Prob (F-statistic):               0.00
#Time:                        18:24:26   Log-Likelihood:                 1695.2
#No. Observations:                7200   AIC:                            -3384.
#Df Residuals:                    7197   BIC:                            -3364.
#Df Model:                           2                                         
#Covariance Type:            nonrobust                                         
#===================================================================================
#                      coef    std err          t      P>|t|      [0.025      0.975]
#-----------------------------------------------------------------------------------
#Intercept          -4.4839      0.020   -222.477      0.000      -4.523      -4.444
#FD_initial          0.1060      0.001     88.433      0.000       0.104       0.108
#DOC_initial_int     0.6206      0.005    119.693      0.000       0.610       0.631
#==============================================================================
#Omnibus:                      757.979   Durbin-Watson:                   0.375
#Prob(Omnibus):                  0.000   Jarque-Bera (JB):             1412.798
#Skew:                          -0.704   Prob(JB):                    1.64e-307
#Kurtosis:                       4.651   Cond. No.                         64.2
#==============================================================================
#
#Notes:
#[1] Standard Errors assume that the covariance matrix of the errors is correctly specified.
#%%
dec_rec_fd_doc_m, dec_rec_fd_doc_r = dec_cont_f_fd_doc(lognona, 'reduced_C')
#                            OLS Regression Results                            
#==============================================================================
#Dep. Variable:            decay_const   R-squared:                       0.694
#Model:                            OLS   Adj. R-squared:                  0.694
#Method:                 Least Squares   F-statistic:                     8154.
#Date:                Thu, 19 Jan 2023   Prob (F-statistic):               0.00
#Time:                        17:41:19   Log-Likelihood:                 494.76
#No. Observations:                7198   AIC:                            -983.5
#Df Residuals:                    7195   BIC:                            -962.9
#Df Model:                           2                                         
#Covariance Type:            nonrobust                                         
#===================================================================================
#                      coef    std err          t      P>|t|      [0.025      0.975]
#-----------------------------------------------------------------------------------
#Intercept          -4.4989      0.024   -188.914      0.000      -4.546      -4.452
#FD_initial          0.1089      0.001     76.876      0.000       0.106       0.112
#DOC_initial_int     0.6248      0.006    101.985      0.000       0.613       0.637
#==============================================================================
#Omnibus:                      997.585   Durbin-Watson:                   0.391
#Prob(Omnibus):                  0.000   Jarque-Bera (JB):             1965.408
#Skew:                          -0.865   Prob(JB):                         0.00
#Kurtosis:                       4.887   Cond. No.                         64.2
#==============================================================================
#
#Notes:
#[1] Standard Errors assume that the covariance matrix of the errors is correctly specified.
#%%
dec_oxc_fd_doc_m,dec_oxc_fd_doc_r = dec_cont_f_fd_doc(lognona, 'oxidized_C')
#                            OLS Regression Results                            
#==============================================================================
#Dep. Variable:            decay_const   R-squared:                       0.737
#Model:                            OLS   Adj. R-squared:                  0.737
#Method:                 Least Squares   F-statistic:                     9819.
#Date:                Thu, 19 Jan 2023   Prob (F-statistic):               0.00
#Time:                        17:41:28   Log-Likelihood:                 1426.1
#No. Observations:                6998   AIC:                            -2846.
#Df Residuals:                    6995   BIC:                            -2826.
#Df Model:                           2                                         
#Covariance Type:            nonrobust                                         
#===================================================================================
#                      coef    std err          t      P>|t|      [0.025      0.975]
#-----------------------------------------------------------------------------------
#Intercept          -4.4994      0.021   -213.265      0.000      -4.541      -4.458
#FD_initial          0.1020      0.001     81.163      0.000       0.100       0.104
#DOC_initial_int     0.6202      0.005    114.230      0.000       0.610       0.631
#==============================================================================
#Omnibus:                     1142.125   Durbin-Watson:                   0.412
#Prob(Omnibus):                  0.000   Jarque-Bera (JB):             3081.968
#Skew:                          -0.885   Prob(JB):                         0.00
#Kurtosis:                       5.726   Cond. No.                         64.0
#==============================================================================
#
#Notes:
#[1] Standard Errors assume that the covariance matrix of the errors is correctly specified.
#%%
dec_toc_fcbvdoc_m, dec_toc_fcbvdoc_r = dec_cont_f_fcbvdoc(lognona, 'TOC')
#                            OLS Regression Results                            
#==============================================================================
#Dep. Variable:            decay_const   R-squared:                       0.867
#Model:                            OLS   Adj. R-squared:                  0.867
#Method:                 Least Squares   F-statistic:                     9383.
#Date:                Thu, 19 Jan 2023   Prob (F-statistic):               0.00
#Time:                        17:43:56   Log-Likelihood:                 4135.1
#No. Observations:                7200   AIC:                            -8258.
#Df Residuals:                    7194   BIC:                            -8217.
#Df Model:                           5                                         
#Covariance Type:            nonrobust                                         
#===================================================================================
#                      coef    std err          t      P>|t|      [0.025      0.975]
#-----------------------------------------------------------------------------------
#Intercept          -2.8614      0.036    -79.004      0.000      -2.932      -2.790
#FD_initial          0.3933      0.005     78.283      0.000       0.383       0.403
#carbon_species     -0.3218      0.011    -28.202      0.000      -0.344      -0.299
#biomass_species    -0.1165      0.006    -19.775      0.000      -0.128      -0.105
#Variance           -0.6905      0.011    -61.427      0.000      -0.712      -0.668
#DOC_initial_int     0.5946      0.004    160.895      0.000       0.587       0.602
#==============================================================================
#Omnibus:                      434.992   Durbin-Watson:                   0.648
#Prob(Omnibus):                  0.000   Jarque-Bera (JB):             1526.677
#Skew:                          -0.226   Prob(JB):                         0.00
#Kurtosis:                       5.210   Cond. No.                         175.
#...
#==============================================================================
#
#Notes:
#[1] Standard Errors assume that the covariance matrix of the errors is correctly specified.
#%%
dec_doc_fcbvdoc_m,dec_doc_fcbvdoc_r = dec_cont_f_fcbvdoc(lognona, 'DOC')
#                            OLS Regression Results                            
#==============================================================================
#Dep. Variable:            decay_const   R-squared:                       0.872
#Model:                            OLS   Adj. R-squared:                  0.872
#Method:                 Least Squares   F-statistic:                     9783.
#Date:                Thu, 19 Jan 2023   Prob (F-statistic):               0.00
#Time:                        17:44:13   Log-Likelihood:                 4030.3
#No. Observations:                7200   AIC:                            -8049.
#Df Residuals:                    7194   BIC:                            -8007.
#Df Model:                           5                                         
#Covariance Type:            nonrobust                                         
#===================================================================================
#                      coef    std err          t      P>|t|      [0.025      0.975]
#-----------------------------------------------------------------------------------
#Intercept          -2.8696      0.037    -78.086      0.000      -2.942      -2.798
#FD_initial          0.4010      0.005     78.653      0.000       0.391       0.411
#carbon_species     -0.3170      0.012    -27.383      0.000      -0.340      -0.294
#biomass_species    -0.0888      0.006    -14.865      0.000      -0.101      -0.077
#Variance           -0.7054      0.011    -61.846      0.000      -0.728      -0.683
#DOC_initial_int     0.6206      0.004    165.511      0.000       0.613       0.628
#==============================================================================
#Omnibus:                      509.259   Durbin-Watson:                   0.642
#Prob(Omnibus):                  0.000   Jarque-Bera (JB):             1920.533
#Skew:                          -0.275   Prob(JB):                         0.00
#Kurtosis:                       5.470   Cond. No.                         175.
#...
#==============================================================================
#
#Notes:
#[1] Standard Errors assume that the covariance matrix of the errors is correctly specified.
#%%
dec_rec_fcbvdoc_m,dec_rec_fcbvdoc_r = dec_cont_f_fcbvdoc(lognona, 'reduced_C')
#                            OLS Regression Results                            
#==============================================================================
#Dep. Variable:            decay_const   R-squared:                       0.804
#Model:                            OLS   Adj. R-squared:                  0.804
#Method:                 Least Squares   F-statistic:                     5914.
#Date:                Thu, 19 Jan 2023   Prob (F-statistic):               0.00
#Time:                        17:44:28   Log-Likelihood:                 2106.2
#No. Observations:                7198   AIC:                            -4200.
#Df Residuals:                    7192   BIC:                            -4159.
#Df Model:                           5                                         
#Covariance Type:            nonrobust                                         
#===================================================================================
#                      coef    std err          t      P>|t|      [0.025      0.975]
#-----------------------------------------------------------------------------------
#Intercept          -2.8766      0.048    -59.878      0.000      -2.971      -2.782
#FD_initial          0.4048      0.007     60.681      0.000       0.392       0.418
#carbon_species     -0.3021      0.015    -19.950      0.000      -0.332      -0.272
#biomass_species    -0.1054      0.008    -13.501      0.000      -0.121      -0.090
#Variance           -0.7093      0.015    -47.525      0.000      -0.739      -0.680
#DOC_initial_int     0.6251      0.005    127.606      0.000       0.615       0.635
#==============================================================================
#Omnibus:                      736.074   Durbin-Watson:                   0.578
#Prob(Omnibus):                  0.000   Jarque-Bera (JB):             2549.678
#Skew:                          -0.498   Prob(JB):                         0.00
#Kurtosis:                       5.740   Cond. No.                         175.
#...
#==============================================================================
#
#Notes:
#[1] Standard Errors assume that the covariance matrix of the errors is correctly specified.
#%%
dec_oxc_fcbvdoc_m,dec_oxc_fcbvdoc_r = dec_cont_f_fcbvdoc(lognona, 'oxidized_C')
#                            OLS Regression Results                            
#==============================================================================
#Dep. Variable:            decay_const   R-squared:                       0.844
#Model:                            OLS   Adj. R-squared:                  0.844
#Method:                 Least Squares   F-statistic:                     7575.
#Date:                Thu, 19 Jan 2023   Prob (F-statistic):               0.00
#Time:                        17:44:40   Log-Likelihood:                 3252.4
#No. Observations:                6998   AIC:                            -6493.
#Df Residuals:                    6992   BIC:                            -6452.
#Df Model:                           5                                         
#Covariance Type:            nonrobust                                         
#===================================================================================
#                      coef    std err          t      P>|t|      [0.025      0.975]
#-----------------------------------------------------------------------------------
#Intercept          -2.9422      0.041    -71.738      0.000      -3.023      -2.862
#FD_initial          0.3900      0.006     68.572      0.000       0.379       0.401
#carbon_species     -0.3175      0.013    -24.638      0.000      -0.343      -0.292
#biomass_species    -0.0676      0.007    -10.093      0.000      -0.081      -0.054
#Variance           -0.6869      0.013    -53.947      0.000      -0.712      -0.662
#DOC_initial_int     0.6204      0.004    148.309      0.000       0.612       0.629
#==============================================================================
#Omnibus:                     1139.692   Durbin-Watson:                   0.629
#Prob(Omnibus):                  0.000   Jarque-Bera (JB):             6095.714
#Skew:                          -0.677   Prob(JB):                         0.00
#Kurtosis:                       7.367   Cond. No.                         175.
#...
#==============================================================================
#
#Notes:
#[1] Standard Errors assume that the covariance matrix of the errors is correctly specified.
#%%
#Run prediction for community evolution:
def fd_ratio_f_fd_doc (df):
    df_sub  = df[df.C_pool=='DOC']
    model = smf.ols('FD_ratio ~ FD_initial + DOC_initial_int', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results

def fd_ratio_f_cbvdoc (df):
    df_sub  = df[df.C_pool=='DOC']
    model = smf.ols('FD_ratio ~ carbon_species + biomass_species + Variance + DOC_initial_int', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results

def fd_ratio_f_cbvosdoc (df):
    df_sub  = df[df.C_pool=='DOC']
    model = smf.ols('FD_ratio ~ carbon_species + biomass_species + Variance + NOSC_initial + DOC_initial_int', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results

def fd_ratio_f_fcbvdoc (df):
    df_sub  = df[df.C_pool=='DOC']
    model = smf.ols('FD_ratio ~ FD_initial + carbon_species + biomass_species + Variance + DOC_initial_int', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results

def bio_ratio_f_fd_doc (df):
    df_sub  = df[df.C_pool=='DOC']
    model = smf.ols('Biomass_ratio ~ FD_initial + DOC_initial_int', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results

def bio_ratio_f_cbvdoc (df):
    df_sub  = df[df.C_pool=='DOC']
    model = smf.ols('Biomass_ratio ~ carbon_species + biomass_species + Variance + DOC_initial_int', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results

def bio_ratio_f_cbvosdoc (df):
    df_sub  = df[df.C_pool=='DOC']
    model = smf.ols('Biomass_ratio ~ carbon_species + biomass_species + Variance + NOSC_initial + DOC_initial_int', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results

def bio_ratio_f_fcbvdoc (df):
    df_sub  = df[df.C_pool=='DOC']
    model = smf.ols('Biomass_ratio ~ FD_initial + carbon_species + biomass_species + Variance + DOC_initial_int', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results
#%%
fd_ratio_cbvosdoc_m, fd_ratio_cbvosdoc_r = fd_ratio_f_cbvosdoc(lognona)
br_cbvosdoc_m, br_cbvosdoc_r = bio_ratio_f_cbvosdoc(lognona)
#%%
#%%
fd_ratio_fd_doc_m, fd_ratio_fd_doc_r = fd_ratio_f_fd_doc (lognona)
#                            OLS Regression Results                            
#==============================================================================
#Dep. Variable:               FD_ratio   R-squared:                       0.465
#Model:                            OLS   Adj. R-squared:                  0.465
#Method:                 Least Squares   F-statistic:                     3130.
#Date:                Thu, 19 Jan 2023   Prob (F-statistic):               0.00
#Time:                        17:49:51   Log-Likelihood:                 8674.6
#No. Observations:                7200   AIC:                        -1.734e+04
#Df Residuals:                    7197   BIC:                        -1.732e+04
#Df Model:                           2                                         
#Covariance Type:            nonrobust                                         
#===================================================================================
#                      coef    std err          t      P>|t|      [0.025      0.975]
#-----------------------------------------------------------------------------------
#Intercept           0.2201      0.008     28.787      0.000       0.205       0.235
#FD_initial          0.0357      0.000     78.585      0.000       0.035       0.037
#DOC_initial_int     0.0181      0.002      9.190      0.000       0.014       0.022
#==============================================================================
#Omnibus:                     1728.877   Durbin-Watson:                   0.229
#Prob(Omnibus):                  0.000   Jarque-Bera (JB):             4829.437
#Skew:                           1.266   Prob(JB):                         0.00
#Kurtosis:                       6.113   Cond. No.                         64.2
#==============================================================================

#Notes:
#[1] Standard Errors assume that the covariance matrix of the errors is correctly specified.
#%%
fd_ratio_cbvdoc_m, fd_ratio_cbvdoc_r = fd_ratio_f_cbvdoc (lognona)
#                            OLS Regression Results                            
#==============================================================================
#Dep. Variable:               FD_ratio   R-squared:                       0.494
#Model:                            OLS   Adj. R-squared:                  0.494
#Method:                 Least Squares   F-statistic:                     1755.
#Date:                Thu, 19 Jan 2023   Prob (F-statistic):               0.00
#Time:                        17:50:04   Log-Likelihood:                 8872.6
#No. Observations:                7200   AIC:                        -1.774e+04
#Df Residuals:                    7195   BIC:                        -1.770e+04
#Df Model:                           4                                         
#Covariance Type:            nonrobust                                         
#===================================================================================
#                      coef    std err          t      P>|t|      [0.025      0.975]
#-----------------------------------------------------------------------------------
#Intercept           0.0248      0.008      3.014      0.003       0.009       0.041
#carbon_species     -0.0261      0.003     -9.324      0.000      -0.032      -0.021
#biomass_species     0.0588      0.003     20.195      0.000       0.053       0.064
#Variance            0.0838      0.001     80.218      0.000       0.082       0.086
#DOC_initial_int     0.0181      0.002      9.445      0.000       0.014       0.022
#==============================================================================
#Omnibus:                     1510.336   Durbin-Watson:                   0.272
#Prob(Omnibus):                  0.000   Jarque-Bera (JB):             3704.691
#Skew:                           1.158   Prob(JB):                         0.00
#Kurtosis:                       5.643   Cond. No.                         42.2
#==============================================================================

#Notes:
#[1] Standard Errors assume that the covariance matrix of the errors is correctly specified.
#%%
fd_ratio_fcbvdoc_m,fd_ratio_fcbvdoc_r = fd_ratio_f_fcbvdoc (lognona)
#                            OLS Regression Results                            
#==============================================================================
#Dep. Variable:               FD_ratio   R-squared:                       0.644
#Model:                            OLS   Adj. R-squared:                  0.644
#Method:                 Least Squares   F-statistic:                     2604.
#Date:                Thu, 19 Jan 2023   Prob (F-statistic):               0.00
#Time:                        17:50:12   Log-Likelihood:                 10141.
#No. Observations:                7200   AIC:                        -2.027e+04
#Df Residuals:                    7194   BIC:                        -2.023e+04
#Df Model:                           5                                         
#Covariance Type:            nonrobust                                         
#===================================================================================
#                      coef    std err          t      P>|t|      [0.025      0.975]
#-----------------------------------------------------------------------------------
#Intercept           0.8041      0.016     51.126      0.000       0.773       0.835
#FD_initial          0.1203      0.002     55.130      0.000       0.116       0.125
#carbon_species     -0.2665      0.005    -53.802      0.000      -0.276      -0.257
#biomass_species     0.0167      0.003      6.537      0.000       0.012       0.022
#Variance           -0.1809      0.005    -37.060      0.000      -0.190      -0.171
#DOC_initial_int     0.0181      0.002     11.264      0.000       0.015       0.021
#==============================================================================
#Omnibus:                     2144.202   Durbin-Watson:                   0.323
#Prob(Omnibus):                  0.000   Jarque-Bera (JB):             9144.918
#Skew:                           1.404   Prob(JB):                         0.00
#Kurtosis:                       7.753   Cond. No.                         175.
#...
#==============================================================================

#Notes:
#[1] Standard Errors assume that the covariance matrix of the errors is correctly specified.
#%%
br_fd_doc_m, br_fd_doc_r = bio_ratio_f_fd_doc (lognona)
#                            OLS Regression Results                            
#==============================================================================
#Dep. Variable:          Biomass_ratio   R-squared:                       0.403
#Model:                            OLS   Adj. R-squared:                  0.403
#Method:                 Least Squares   F-statistic:                     2430.
#Date:                Thu, 19 Jan 2023   Prob (F-statistic):               0.00
#Time:                        17:50:23   Log-Likelihood:                 16722.
#No. Observations:                7200   AIC:                        -3.344e+04
#Df Residuals:                    7197   BIC:                        -3.342e+04
#Df Model:                           2                                         
#Covariance Type:            nonrobust                                         
#===================================================================================
#                      coef    std err          t      P>|t|      [0.025      0.975]
#-----------------------------------------------------------------------------------
#Intercept          -0.0461      0.003    -18.457      0.000      -0.051      -0.041
#FD_initial          0.0054      0.000     36.226      0.000       0.005       0.006
#DOC_initial_int     0.0383      0.001     59.565      0.000       0.037       0.040
#==============================================================================
#Omnibus:                     3392.357   Durbin-Watson:                   0.798
#Prob(Omnibus):                  0.000   Jarque-Bera (JB):            41182.897
#Skew:                          -1.943   Prob(JB):                         0.00
#Kurtosis:                      14.053   Cond. No.                         64.2
#==============================================================================
#
#Notes:
#[1] Standard Errors assume that the covariance matrix of the errors is correctly specified.
#%%
br_cbvdoc_m,br_cbvdoc_r  = bio_ratio_f_cbvdoc (lognona)
#                            OLS Regression Results                            
#==============================================================================
#Dep. Variable:          Biomass_ratio   R-squared:                       0.487
#Model:                            OLS   Adj. R-squared:                  0.487
#Method:                 Least Squares   F-statistic:                     1709.
#Date:                Thu, 19 Jan 2023   Prob (F-statistic):               0.00
#Time:                        17:50:32   Log-Likelihood:                 17268.
#No. Observations:                7200   AIC:                        -3.453e+04
#Df Residuals:                    7195   BIC:                        -3.449e+04
#Df Model:                           4                                         
#Covariance Type:            nonrobust                                         
#===================================================================================
#                      coef    std err          t      P>|t|      [0.025      0.975]
#-----------------------------------------------------------------------------------
#Intercept          -0.1272      0.003    -49.664      0.000      -0.132      -0.122
#carbon_species      0.0172      0.001     19.723      0.000       0.016       0.019
#biomass_species     0.0349      0.001     38.449      0.000       0.033       0.037
#Variance            0.0094      0.000     28.968      0.000       0.009       0.010
#DOC_initial_int     0.0383      0.001     64.253      0.000       0.037       0.039
#==============================================================================
#Omnibus:                     1957.826   Durbin-Watson:                   0.924
#Prob(Omnibus):                  0.000   Jarque-Bera (JB):            34223.819
#Skew:                          -0.845   Prob(JB):                         0.00
#Kurtosis:                      13.546   Cond. No.                         42.2
#==============================================================================
#
#Notes:
#[1] Standard Errors assume that the covariance matrix of the errors is correctly specified.
#%%
br_fcbvdoc_m,br_fcbvdoc_r = bio_ratio_f_fcbvdoc (lognona)
#                            OLS Regression Results                            
#==============================================================================
#Dep. Variable:          Biomass_ratio   R-squared:                       0.523
#Model:                            OLS   Adj. R-squared:                  0.523
#Method:                 Least Squares   F-statistic:                     1577.
#Date:                Thu, 19 Jan 2023   Prob (F-statistic):               0.00
#Time:                        17:50:39   Log-Likelihood:                 17529.
#No. Observations:                7200   AIC:                        -3.505e+04
#Df Residuals:                    7194   BIC:                        -3.500e+04
#Df Model:                           5                                         
#Covariance Type:            nonrobust                                         
#===================================================================================
#                      coef    std err          t      P>|t|      [0.025      0.975]
#-----------------------------------------------------------------------------------
#Intercept          -0.0095      0.006     -1.690      0.091      -0.021       0.002
#FD_initial          0.0182      0.001     23.233      0.000       0.017       0.020
#carbon_species     -0.0191      0.002    -10.744      0.000      -0.023      -0.016
#biomass_species     0.0285      0.001     31.114      0.000       0.027       0.030
#Variance           -0.0306      0.002    -17.464      0.000      -0.034      -0.027
#DOC_initial_int     0.0383      0.001     66.615      0.000       0.037       0.039
#==============================================================================
#Omnibus:                     1824.628   Durbin-Watson:                   0.969
#Prob(Omnibus):                  0.000   Jarque-Bera (JB):            32890.278
#Skew:                          -0.746   Prob(JB):                         0.00
#Kurtosis:                      13.364   Cond. No.                         175.
#...
#==============================================================================

#Notes:
#[1] Standard Errors assume that the covariance matrix of the errors is correctly specified.
#%%
br_cbvdoc_m,br_cbvdoc_r = bio_ratio_f_cbvdoc (lognona)
coefs_fcbvdoc = pd.DataFrame()
coefs_fcbvdoc['Features'] = ['#carbon', '#biomass', 'Variance', 'C availability']
coefs_fcbvdoc['dec_toc'] = dec_toc_cbvdoc_r.params[1:].values
coefs_fcbvdoc['dec_doc'] = dec_doc_cbvdoc_r.params[1:].values
coefs_fcbvdoc['dec_rec'] = dec_rec_cbvdoc_r.params[1:].values
coefs_fcbvdoc['dec_oxc'] = dec_oxc_cbvdoc_r.params[1:].values
coefs_fcbvdoc['FD_ratio'] = fd_ratio_cbvdoc_r.params[1:].values
coefs_fcbvdoc['Biomass_ratio'] = br_cbvdoc_r.params[1:].values
coefs_fcbvdoc['dec_stop_toc'] = dec_stop_toc_cbvdoc_r.params[1:].values
coefs_fcbvdoc['dec_stop_doc'] = dec_stop_doc_cbvdoc_r.params[1:].values
coefs_fcbvdoc['dec_stop_rec'] = dec_stop_rec_cbvdoc_r.params[1:].values
coefs_fcbvdoc['dec_stop_oxc'] = dec_stop_oxc_cbvdoc_r.params[1:].values
coefs_fcbvdoc.set_index('Features', inplace=True)
#%%
sns.heatmap(coefs_fcbvdoc.transpose(), center = 0, annot = True, fmt=".2f", cmap = 'vlag')
#%%
coefs_fd_doc = pd.DataFrame()
coefs_fd_doc['Features'] = ['Initial functional diversity', 'C availability']
coefs_fd_doc['dec_toc'] = dec_toc_fd_doc_r.params[1:].values
coefs_fd_doc['dec_doc'] = dec_doc_fd_doc_r.params[1:].values
coefs_fd_doc['dec_rec'] = dec_rec_fd_doc_r.params[1:].values
coefs_fd_doc['dec_oxc'] = dec_oxc_fd_doc_r.params[1:].values
coefs_fd_doc['FD_ratio'] = fd_ratio_fd_doc_r.params[1:].values
coefs_fd_doc['Biomass_ratio'] = br_fd_doc_r.params[1:].values
coefs_fd_doc['dec_stop_toc'] = dec_stop_toc_fd_doc_r.params[1:].values
coefs_fd_doc['dec_stop_doc'] = dec_stop_doc_fd_doc_r.params[1:].values
coefs_fd_doc['dec_stop_rec'] = dec_stop_rec_fd_doc_r.params[1:].values
coefs_fd_doc['dec_stop_oxc'] = dec_stop_oxc_fd_doc_r.params[1:].values
coefs_fd_doc.set_index('Features', inplace=True)
#%%
sns.heatmap(coefs_fd_doc.transpose(), annot = True, fmt=".2f", center = 0, cmap = 'vlag')
#%%
#TOC plot prediction with scatter plots
x1 = lognona[lognona.C_pool=='TOC'].FD_initial
x2 = lognona[lognona.C_pool=='TOC'].DOC_initial_int
y = lognona[lognona.C_pool=='TOC'].decay_const
x1_sort = x1#np.sort(x1.sort_values())
x2_sort = x2#p.sort(x2.sort_values())
y_pred = dec_toc_fd_doc_r.params[0]+dec_toc_fd_doc_r.params[1]*x1_sort+dec_toc_fd_doc_r.params[2]*x2_sort
fig, axes = plt.subplots(nrows = 1, ncols = 2, sharey = True, figsize = (6,3))
axes[0].scatter(x1,y)
axes[0].scatter(x1,y_pred, c='r')
axes[1].scatter(x2,y)
axes[1].scatter(x2,y_pred,c='r')
axes[0].set_xlabel("log f")
axes[1].set_xlabel("log (initial C\navailability)")
fig.supylabel("log (decay constant)")
fig.suptitle("TOC")
#%%
def dec_const_f_fd (df, cpool_test):
    df_sub  = df[df.C_pool==cpool_test]
    model = smf.ols('decay_const ~ FD_initial', df_sub)
    results = model.fit()
    print(results.summary())

    return model, results
#%%
dec_toc_fd_m, dec_toc_fd_r = dec_const_f_fd(lognona, 'TOC')
#%%
x1 = lognona[lognona.C_pool=='TOC'].FD_initial
y = lognona[lognona.C_pool=='TOC'].decay_const
x1_sort = x1#np.sort(x1.sort_values())
x1_pred = dec_toc_fd_r.params[0]+dec_toc_fd_r.params[1]*x1_sort
fig, axes = plt.subplots(nrows = 1, ncols = 1, sharey = True)
axes.scatter(x1,y)
axes.scatter(x1,x1_pred, c='r')
#%%
#Split data as per similarity
lognona_sim = lognona[lognona.FD_initial<=-7]
lognona_mid = lognona[(lognona.FD_initial>-7)&(lognona.FD_initial<=-5)]
lognona_dis = lognona[lognona.FD_initial>-5]
sim_dec_toc_fd_m, sim_dec_toc_fd_r = dec_const_f_fd(lognona_sim, 'TOC')
mid_dec_toc_fd_m, mid_dec_toc_fd_r = dec_const_f_fd(lognona_mid, 'TOC')
dis_dec_toc_fd_m, dis_dec_toc_fd_r = dec_const_f_fd(lognona_dis, 'TOC')
x1 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial<=-7)].FD_initial
y1 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial<=-7)].decay_const
x2 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-7)&(lognona.FD_initial<=-5)].FD_initial
y2 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-7)&(lognona.FD_initial<=-5)].decay_const
x3 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-5)].FD_initial
y3 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-5)].decay_const
x1_sort = x1#np.sort(x1.sort_values())
x1_pred = sim_dec_toc_fd_r.params[0]+sim_dec_toc_fd_r.params[1]*x1
x2_pred = mid_dec_toc_fd_r.params[0]+mid_dec_toc_fd_r.params[1]*x2
x3_pred = dis_dec_toc_fd_r.params[0]+dis_dec_toc_fd_r.params[1]*x3
fig, axes = plt.subplots(nrows = 1, ncols = 3, sharey = True, sharex=False, figsize = (8,3))
axes[0].scatter(x1,y1)
axes[0].scatter(x1,x1_pred, c='r')
axes[1].scatter(x2,y2)
axes[1].scatter(x2,x2_pred, c='r')
axes[2].scatter(x3,y3)
axes[2].scatter(x3,x3_pred, c='r')
fig.supxlabel("log f")
fig.supylabel("log (decay constant)")
axes[0].set_title("Similar\ncommunities")
axes[1].set_title("Transitioning\ncommunities")
axes[2].set_title("Dissimilar\ncommunities")
fig.subplots_adjust(top=0.75)
fig.suptitle("TOC")
#%%
sim_dec_toc_fd_doc_m, sim_dec_toc_fd_doc_r = dec_cont_f_fd_doc(lognona_sim, 'TOC')
mid_dec_toc_fd_doc_m, mid_dec_toc_fd_doc_r = dec_cont_f_fd_doc(lognona_mid, 'TOC')
dis_dec_toc_fd_doc_m, dis_dec_toc_fd_doc_r = dec_cont_f_fd_doc(lognona_dis, 'TOC')
x11 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial<=-7)].FD_initial
x12 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial<=-7)].DOC_initial_int
y1 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial<=-7)].decay_const
x21 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-7)&(lognona.FD_initial<=-5)].FD_initial
x22 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-7)&(lognona.FD_initial<=-5)].DOC_initial_int
y2 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-7)&(lognona.FD_initial<=-5)].decay_const
x31 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-5)].FD_initial
x32 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-5)].DOC_initial_int
y3 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-5)].decay_const
x1_sort = x1#np.sort(x1.sort_values())
y1_pred = sim_dec_toc_fd_doc_r.params[0]+sim_dec_toc_fd_doc_r.params[1]*x11+sim_dec_toc_fd_doc_r.params[2]*x12
y2_pred = mid_dec_toc_fd_doc_r.params[0]+mid_dec_toc_fd_doc_r.params[1]*x21+mid_dec_toc_fd_doc_r.params[2]*x22
y3_pred = dis_dec_toc_fd_doc_r.params[0]+dis_dec_toc_fd_doc_r.params[1]*x31+dis_dec_toc_fd_doc_r.params[2]*x32


fig, axes = plt.subplots(nrows = 2, ncols = 3, sharey = True, sharex=False, figsize = (10,8))
axes[0][0].scatter(x11,y1)
axes[0][0].scatter(x11,y1_pred, c='r')
axes[0][1].scatter(x21,y2)
axes[0][1].scatter(x21,y2_pred, c='r')
axes[0][2].scatter(x31,y3)
axes[0][2].scatter(x31,y3_pred, c='r')
axes[1][0].scatter(x12,y1)
axes[1][0].scatter(x12,y1_pred, c='r')
axes[1][1].scatter(x22,y2)
axes[1][1].scatter(x22,y2_pred, c='r')
axes[1][2].scatter(x32,y3)
axes[1][2].scatter(x32,y3_pred, c='r')
for a in axes[0,:]:
    a.set_xlabel("log f")
for a in axes[1,:]:
    a.set_xlabel("log (C availability)")
fig.supylabel("log (decay constant)")
axes[0,0].set_title("Similar\ncommunities")
axes[0,1].set_title("Transitioning\ncommunities")
axes[0,2].set_title("Dissimilar\ncommunities")
#fig.subplots_adjust(top=0.8)
fig.suptitle("TOC")
#%%
coefs_fd_toc = pd.DataFrame()
coefs_fd_toc['Features'] = ['Initial functional diversity', 'C availability']
coefs_fd_toc['Similar'] =  sim_dec_toc_fd_doc_r.params[1:].values
coefs_fd_toc['Transitioning'] = mid_dec_toc_fd_doc_r.params[1:].values
coefs_fd_toc['Dissimilar'] = dis_dec_toc_fd_doc_r.params[1:].values
coefs_fd_toc.set_index('Features', inplace=True)
sns.heatmap(coefs_fd_toc.transpose(), annot = True, fmt=".2f", center = 0, cmap = 'vlag')
plt.ylabel("Community type")
#%%
sm.stats.anova_lm(sim_dec_toc_fd_doc_r)#, typ=1)
sm.stats.anova_lm(mid_dec_toc_fd_doc_r)#, typ=1)
sm.stats.anova_lm(dis_dec_toc_fd_doc_r)#, typ=1)
#%%
sim_dec_toc_cbvdoc_m, sim_dec_toc_cbvdoc_r = dec_cont_f_cbvdoc(lognona_sim, 'TOC')
mid_dec_toc_cbvdoc_m, mid_dec_toc_cbvdoc_r = dec_cont_f_cbvdoc(lognona_mid, 'TOC')
dis_dec_toc_cbvdoc_m, dis_dec_toc_cbvdoc_r = dec_cont_f_cbvdoc(lognona_dis, 'TOC')
#x11 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial<=-7)].FD_initial
x11 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial<=-7)].carbon_species
x12 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial<=-7)].biomass_species
x13 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial<=-7)].Variance
x14 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial<=-7)].DOC_initial_int
y1 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial<=-7)].decay_const
#x21 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-7)&(lognona.FD_initial<=-5)].FD_initial
x21 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-7)&(lognona.FD_initial<=-5)].carbon_species
x22 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-7)&(lognona.FD_initial<=-5)].biomass_species
x23 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-7)&(lognona.FD_initial<=-5)].Variance
x24 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-7)&(lognona.FD_initial<=-5)].DOC_initial_int
y2 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-7)&(lognona.FD_initial<=-5)].decay_const
#x31 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-5)].FD_initial
x31 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-5)].carbon_species
x32 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-5)].biomass_species
x33 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-5)].Variance
x34 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-5)].DOC_initial_int
y3 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-5)].decay_const
x1_sort = x1#np.sort(x1.sort_values())
y1_pred = sim_dec_toc_cbvdoc_r.params[0]+sim_dec_toc_cbvdoc_r.params[1]*x11+sim_dec_toc_cbvdoc_r.params[2]*x12+sim_dec_toc_cbvdoc_r.params[3]*x13+sim_dec_toc_cbvdoc_r.params[4]*x14
y2_pred = mid_dec_toc_cbvdoc_r.params[0]+mid_dec_toc_cbvdoc_r.params[1]*x21+mid_dec_toc_cbvdoc_r.params[2]*x22+mid_dec_toc_cbvdoc_r.params[3]*x23+mid_dec_toc_cbvdoc_r.params[4]*x24
y3_pred = dis_dec_toc_cbvdoc_r.params[0]+dis_dec_toc_cbvdoc_r.params[1]*x31+dis_dec_toc_cbvdoc_r.params[2]*x32+dis_dec_toc_cbvdoc_r.params[3]*x33+dis_dec_toc_cbvdoc_r.params[4]*x34

fig, axes = plt.subplots(nrows = 4, ncols = 3, sharey = True, sharex=False, figsize = (7,8))
axes[0][0].scatter(x11,y1)
axes[0][0].scatter(x11,y1_pred, c='r')
axes[0][1].scatter(x21,y2)
axes[0][1].scatter(x21,y2_pred, c='r')
axes[0][2].scatter(x31,y3)
axes[0][2].scatter(x31,y3_pred, c='r')
axes[1][0].scatter(x12,y1)
axes[1][0].scatter(x12,y1_pred, c='r')
axes[1][1].scatter(x22,y2)
axes[1][1].scatter(x22,y2_pred, c='r')
axes[1][2].scatter(x32,y3)
axes[1][2].scatter(x32,y3_pred, c='r')
axes[2][0].scatter(x13,y1)
axes[2][0].scatter(x13,y1_pred, c='r')
axes[2][1].scatter(x23,y2)
axes[2][1].scatter(x23,y2_pred, c='r')
axes[2][2].scatter(x33,y3)
axes[2][2].scatter(x33,y3_pred, c='r')

axes[3][0].scatter(x14,y1)
axes[3][0].scatter(x14,y1_pred, c='r')
axes[3][1].scatter(x24,y2)
axes[3][1].scatter(x24,y2_pred, c='r')
axes[3][2].scatter(x34,y3)
axes[3][2].scatter(x34,y3_pred, c='r')
fig.supylabel("log(decay constant)")
for a in axes[0,:]:
    a.set_xlabel("#carbon species")
for a in axes[1,:]:
    a.set_xlabel("#biomass species")
for a in axes[2,:]:
    a.set_xlabel("Variance")
for a in axes[3,:]:
    a.set_xlabel("C availability")
fig.tight_layout()

#%%
coefs_dec_toc_fcbvdoc = pd.DataFrame()
coefs_dec_toc_fcbvdoc['Features'] = ["#carbon", "#biomass", "Variance", 'C availability']
coefs_dec_toc_fcbvdoc['Similar'] =  sim_dec_toc_cbvdoc_r.params[1:].values
coefs_dec_toc_fcbvdoc['Transitioning'] = mid_dec_toc_cbvdoc_r.params[1:].values
coefs_dec_toc_fcbvdoc['Dissimilar'] = dis_dec_toc_cbvdoc_r.params[1:].values
coefs_dec_toc_fcbvdoc.set_index('Features', inplace=True)
sns.heatmap(coefs_dec_toc_fcbvdoc.transpose(), annot = True, fmt=".2f", center = 0, cmap = 'vlag')
plt.ylabel("Community type")
#%%
sim_fd_ratio_cbvdoc_m, sim_fd_ratio_cbvdoc_r = fd_ratio_f_cbvdoc(lognona_sim)
mid_fd_ratio_cbvdoc_m, mid_fd_ratio_cbvdoc_r = fd_ratio_f_cbvdoc(lognona_mid)
dis_fd_ratio_cbvdoc_m, dis_fd_ratio_cbvdoc_r = fd_ratio_f_cbvdoc(lognona_dis)
y1 = lognona[(lognona.C_pool=='DOC')&(lognona.FD_initial<=-7)].FD_ratio
y2 = lognona[(lognona.C_pool=='DOC')&(lognona.FD_initial>-7)&(lognona.FD_initial<=-5)].FD_ratio
y3 = lognona[(lognona.C_pool=='DOC')&(lognona.FD_initial>-5)].FD_ratio
y1_pred = sim_fd_ratio_cbvdoc_r.params[0]+sim_fd_ratio_cbvdoc_r.params[1]*x11+sim_fd_ratio_cbvdoc_r.params[2]*x12+sim_fd_ratio_cbvdoc_r.params[3]*x13+sim_fd_ratio_cbvdoc_r.params[4]*x14
y2_pred = mid_fd_ratio_cbvdoc_r.params[0]+mid_fd_ratio_cbvdoc_r.params[1]*x21+mid_fd_ratio_cbvdoc_r.params[2]*x22+mid_fd_ratio_cbvdoc_r.params[3]*x23+mid_fd_ratio_cbvdoc_r.params[4]*x24
y3_pred = dis_fd_ratio_cbvdoc_r.params[0]+dis_fd_ratio_cbvdoc_r.params[1]*x31+dis_fd_ratio_cbvdoc_r.params[2]*x32+dis_fd_ratio_cbvdoc_r.params[3]*x33+dis_fd_ratio_cbvdoc_r.params[4]*x34

fig, axes = plt.subplots(nrows = 4, ncols = 3, sharey = True, sharex=False, figsize = (7,8))
axes[0][0].scatter(x11,y1)
axes[0][0].scatter(x11,y1_pred, c='r')
axes[0][1].scatter(x21,y2)
axes[0][1].scatter(x21,y2_pred, c='r')
axes[0][2].scatter(x31,y3)
axes[0][2].scatter(x31,y3_pred, c='r')
axes[1][0].scatter(x12,y1)
axes[1][0].scatter(x12,y1_pred, c='r')
axes[1][1].scatter(x22,y2)
axes[1][1].scatter(x22,y2_pred, c='r')
axes[1][2].scatter(x32,y3)
axes[1][2].scatter(x32,y3_pred, c='r')
axes[2][0].scatter(x13,y1)
axes[2][0].scatter(x13,y1_pred, c='r')
axes[2][1].scatter(x23,y2)
axes[2][1].scatter(x23,y2_pred, c='r')
axes[2][2].scatter(x33,y3)
axes[2][2].scatter(x33,y3_pred, c='r')
axes[3][0].scatter(x14,y1)
axes[3][0].scatter(x14,y1_pred, c='r')
axes[3][1].scatter(x24,y2)
axes[3][1].scatter(x24,y2_pred, c='r')
axes[3][2].scatter(x34,y3)
axes[3][2].scatter(x34,y3_pred, c='r')
fig.supylabel("log(gain in f)")
for a in axes[0,:]:
    a.set_xlabel("#carbon species")
for a in axes[1,:]:
    a.set_xlabel("#biomass species")
for a in axes[2,:]:
    a.set_xlabel("Variance")
for a in axes[3,:]:
    a.set_xlabel("C availability")
fig.tight_layout()
#%%
sim_fd_ratio_fd_doc_m, sim_fd_ratio_fd_doc_r = fd_ratio_f_fd_doc(lognona_sim)
mid_fd_ratio_fd_doc_m, mid_fd_ratio_fd_doc_r = fd_ratio_f_fd_doc(lognona_mid)
dis_fd_ratio_fd_doc_m, dis_fd_ratio_fd_doc_r = fd_ratio_f_fd_doc(lognona_dis)
x11 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial<=-7)].FD_initial
x12 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial<=-7)].DOC_initial_int
y1 = lognona[(lognona.C_pool=='DOC')&(lognona.FD_initial<=-7)].FD_ratio
y2 = lognona[(lognona.C_pool=='DOC')&(lognona.FD_initial>-7)&(lognona.FD_initial<=-5)].FD_ratio
y3 = lognona[(lognona.C_pool=='DOC')&(lognona.FD_initial>-5)].FD_ratio
y1_pred = sim_fd_ratio_fd_doc_r.params[0]+sim_fd_ratio_fd_doc_r.params[1]*x11+sim_fd_ratio_fd_doc_r.params[2]*x12
y2_pred = mid_fd_ratio_fd_doc_r.params[0]+mid_fd_ratio_fd_doc_r.params[1]*x21+mid_fd_ratio_fd_doc_r.params[2]*x22
y3_pred = dis_fd_ratio_fd_doc_r.params[0]+dis_fd_ratio_fd_doc_r.params[1]*x31+dis_fd_ratio_fd_doc_r.params[2]*x32

fig, axes = plt.subplots(nrows = 2, ncols = 3, sharey = True, sharex=False, figsize = (6,4))
axes[0][0].scatter(x11,y1)
axes[0][0].scatter(x11,y1_pred, c='r')
axes[0][1].scatter(x21,y2)
axes[0][1].scatter(x21,y2_pred, c='r')
axes[0][2].scatter(x31,y3)
axes[0][2].scatter(x31,y3_pred, c='r')
axes[1][0].scatter(x12,y1)
axes[1][0].scatter(x12,y1_pred, c='r')
axes[1][1].scatter(x22,y2)
axes[1][1].scatter(x22,y2_pred, c='r')
axes[1][2].scatter(x32,y3)
axes[1][2].scatter(x32,y3_pred, c='r')

fig.supylabel("log(gain in f)")
for a in axes[0,:]:
    a.set_xlabel("log (f)")
for a in axes[1,:]:
    a.set_xlabel("C availability")
fig.tight_layout()

plt.figure()
coefs_fd_ratio_fd_doc = pd.DataFrame()
coefs_fd_ratio_fd_doc['Features'] = ['Initial f','C availability']
coefs_fd_ratio_fd_doc['Similar'] =  sim_fd_ratio_fd_doc_r.params[1:].values
coefs_fd_ratio_fd_doc['Transitioning'] = mid_fd_ratio_fd_doc_r.params[1:].values
coefs_fd_ratio_fd_doc['Dissimilar'] = dis_fd_ratio_fd_doc_r.params[1:].values
coefs_fd_ratio_fd_doc.set_index('Features', inplace=True)
sns.heatmap(coefs_fd_ratio_fd_doc.transpose(), annot = True, fmt=".2f", center = 0, cmap = 'vlag')
plt.ylabel("Community type")
#%%
#%%
sim_br_ratio_cbvdoc_m, sim_br_ratio_cbvdoc_r = bio_ratio_f_cbvdoc(lognona_sim)
mid_br_ratio_cbvdoc_m, mid_br_ratio_cbvdoc_r = bio_ratio_f_cbvdoc(lognona_mid)
dis_br_ratio_cbvdoc_m, dis_br_ratio_cbvdoc_r = bio_ratio_f_cbvdoc(lognona_dis)
y1 = lognona[(lognona.C_pool=='DOC')&(lognona.FD_initial<=-7)].Biomass_ratio
y2 = lognona[(lognona.C_pool=='DOC')&(lognona.FD_initial>-7)&(lognona.FD_initial<=-5)].Biomass_ratio
y3 = lognona[(lognona.C_pool=='DOC')&(lognona.FD_initial>-5)].Biomass_ratio
y1_pred = sim_br_ratio_cbvdoc_r.params[0]+sim_br_ratio_cbvdoc_r.params[1]*x11+sim_br_ratio_cbvdoc_r.params[2]*x12+sim_br_ratio_cbvdoc_r.params[3]*x13+sim_br_ratio_cbvdoc_r.params[4]*x14
y2_pred = mid_br_ratio_cbvdoc_r.params[0]+mid_br_ratio_cbvdoc_r.params[1]*x21+mid_br_ratio_cbvdoc_r.params[2]*x22+mid_br_ratio_cbvdoc_r.params[3]*x23+mid_br_ratio_cbvdoc_r.params[4]*x24
y3_pred = dis_br_ratio_cbvdoc_r.params[0]+dis_br_ratio_cbvdoc_r.params[1]*x31+dis_br_ratio_cbvdoc_r.params[2]*x32+dis_br_ratio_cbvdoc_r.params[3]*x33+dis_br_ratio_cbvdoc_r.params[4]*x34

fig, axes = plt.subplots(nrows = 4, ncols = 3, sharey = True, sharex=False, figsize = (7,8))
axes[0][0].scatter(x11,y1)
axes[0][0].scatter(x11,y1_pred, c='r')
axes[0][1].scatter(x21,y2)
axes[0][1].scatter(x21,y2_pred, c='r')
axes[0][2].scatter(x31,y3)
axes[0][2].scatter(x31,y3_pred, c='r')
axes[1][0].scatter(x12,y1)
axes[1][0].scatter(x12,y1_pred, c='r')
axes[1][1].scatter(x22,y2)
axes[1][1].scatter(x22,y2_pred, c='r')
axes[1][2].scatter(x32,y3)
axes[1][2].scatter(x32,y3_pred, c='r')
axes[2][0].scatter(x13,y1)
axes[2][0].scatter(x13,y1_pred, c='r')
axes[2][1].scatter(x23,y2)
axes[2][1].scatter(x23,y2_pred, c='r')
axes[2][2].scatter(x33,y3)
axes[2][2].scatter(x33,y3_pred, c='r')
axes[3][0].scatter(x14,y1)
axes[3][0].scatter(x14,y1_pred, c='r')
axes[3][1].scatter(x24,y2)
axes[3][1].scatter(x24,y2_pred, c='r')
axes[3][2].scatter(x34,y3)
axes[3][2].scatter(x34,y3_pred, c='r')
fig.supylabel("log(gain in biomass)")
for a in axes[0,:]:
    a.set_xlabel("#carbon species")
for a in axes[1,:]:
    a.set_xlabel("#biomass species")
for a in axes[2,:]:
    a.set_xlabel("Variance")
for a in axes[3,:]:
    a.set_xlabel("C availability")
fig.tight_layout()

#%%
coefs_br_ratio_cbvdoc = pd.DataFrame()
coefs_br_ratio_cbvdoc['Features'] = ["#carbon", "#biomass", "Variance", 'C availability']
coefs_br_ratio_cbvdoc['Similar'] =  sim_br_ratio_cbvdoc_r.params[1:].values
coefs_br_ratio_cbvdoc['Transitioning'] = mid_br_ratio_cbvdoc_r.params[1:].values
coefs_br_ratio_cbvdoc['Dissimilar'] = dis_br_ratio_cbvdoc_r.params[1:].values
coefs_br_ratio_cbvdoc.set_index('Features', inplace=True)
sns.heatmap(coefs_br_ratio_cbvdoc.transpose(), annot = True, fmt=".2f", center = 0, cmap = 'vlag')
plt.ylabel("Community type")

#%%
sim_br_ratio_fd_doc_m, sim_br_ratio_fd_doc_r = bio_ratio_f_fd_doc(lognona_sim)
mid_br_ratio_fd_doc_m, mid_br_ratio_fd_doc_r = bio_ratio_f_fd_doc(lognona_mid)
dis_br_ratio_fd_doc_m, dis_br_ratio_fd_doc_r = bio_ratio_f_fd_doc(lognona_dis)
x11 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial<=-7)].FD_initial
x12 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial<=-7)].DOC_initial_int
x21 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-7)&(lognona.FD_initial<=-5)].FD_initial
x22 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-7)&(lognona.FD_initial<=-5)].DOC_initial_int
x31 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-5)].FD_initial
x32 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-5)].DOC_initial_int
y1 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial<=-7)].Biomass_ratio
y2 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-7)&(lognona.FD_initial<=-5)].Biomass_ratio
y3 = lognona[(lognona.C_pool=='TOC')&(lognona.FD_initial>-5)].Biomass_ratio
y1_pred = sim_br_ratio_fd_doc_r.params[0]+sim_br_ratio_fd_doc_r.params[1]*x11+sim_br_ratio_fd_doc_r.params[2]*x12
y2_pred = mid_br_ratio_fd_doc_r.params[0]+mid_br_ratio_fd_doc_r.params[1]*x21+mid_br_ratio_fd_doc_r.params[2]*x22
y3_pred = dis_br_ratio_fd_doc_r.params[0]+dis_br_ratio_fd_doc_r.params[1]*x31+dis_br_ratio_fd_doc_r.params[2]*x32

fig, axes = plt.subplots(nrows = 2, ncols = 3, sharey = True, sharex=False, figsize = (6,4))
axes[0][0].scatter(x11,y1)
axes[0][0].scatter(x11,y1_pred, c='r')
axes[0][1].scatter(x21,y2)
axes[0][1].scatter(x21,y2_pred, c='r')
axes[0][2].scatter(x31,y3)
axes[0][2].scatter(x31,y3_pred, c='r')
axes[1][0].scatter(x12,y1)
axes[1][0].scatter(x12,y1_pred, c='r')
axes[1][1].scatter(x22,y2)
axes[1][1].scatter(x22,y2_pred, c='r')
axes[1][2].scatter(x31,y3)
axes[1][2].scatter(x31,y3_pred, c='r')
fig.supylabel("log(gain in biomass)")
for a in axes[0,:]:
    a.set_xlabel("log f")
for a in axes[1,:]:
    a.set_xlabel("C availability")
fig.tight_layout()

plt.figure()
coefs_br_ratio_fd_doc = pd.DataFrame()
coefs_br_ratio_fd_doc['Features'] = ['Initial f','C availability']
coefs_br_ratio_fd_doc['Similar'] =  sim_br_ratio_fd_doc_r.params[1:].values
coefs_br_ratio_fd_doc['Transitioning'] = mid_br_ratio_fd_doc_r.params[1:].values
coefs_br_ratio_fd_doc['Dissimilar'] = dis_br_ratio_fd_doc_r.params[1:].values
coefs_br_ratio_fd_doc.set_index('Features', inplace=True)
sns.heatmap(coefs_br_ratio_fd_doc.transpose(), annot = True, fmt=".2f", center = 0, cmap = 'vlag')
plt.ylabel("Community type")

#%%
fd_docfd = smf.ols('decay_const ~ FD_initial + DOC_initial_int*FD_initial', lognona)
fd_docfd = fd_docfd.fit()
print(fd_docfd.summary())
results_dic.update({'fd_docfd_aic': fd_docfd.aic})
#%%
cabvdoc = smf.ols('decay_const ~ carbon_species + biomass_species + activity + Variance + DOC_initial_int', lognona)
cabvdoc = cabvdoc.fit()
print(cabvdoc.summary())
results_dic.update({'cabvdoc_aic': cabvdoc.aic})
#%%
fd_doc_docfd = smf.ols('decay_const ~ FD_initial + DOC_initial_int + DOC_initial_int:FD_initial', lognona)
fd_doc_docfd = fd_doc_docfd.fit()
print(fd_doc_docfd.summary())
results_dic.update({'fd_doc_docfd_aic': fd_doc_docfd.aic})
#%%
cabvdoc_doca = smf.ols('decay_const ~ carbon_species + biomass_species + activity + Variance + DOC_initial_int + DOC_initial_int:activity', lognona)
cabvdoc_doca = cabvdoc_doca.fit()
print(cabvdoc_doca.summary())
results_dic.update({'cabvdoc_doca_aic': cabvdoc_doca.aic})
#%%
cabv_doca = smf.ols('decay_const ~ carbon_species + biomass_species + activity + Variance + DOC_initial_int:activity', lognona)
cabv_doca = cabv_doca.fit()
print(cabv_doca.summary())
results_dic.update({'cabv_doca_aic': cabv_doca.aic})

#%%
# Identify the minimum FD required for the FD_ratio>1.5
act_full_fd_ordered = act_full.sort_values(by=['FD_ratio', 'FD_initial'])
biomass_ratio_1 = act_full_fd_ordered.index[act_full_fd_ordered.FD_ratio>1.100000005]
print(act_full_fd_ordered.loc[biomass_ratio_1]['FD_initial'].min())
biomass_ratio_1_5 = act_full_fd_ordered.index[act_full_fd_ordered.FD_ratio>1.5000005]
print(act_full_fd_ordered.loc[biomass_ratio_1_5]['FD_initial'].min())
#%%
# Identify the minimum FD required for the biomass_ratio>1.5 in fully active communities
act_full_bio_ordered = act_full.sort_values(by=['Biomass_ratio', 'FD_initial'])
biomass_ratio_1 = act_full_bio_ordered.index[act_full_bio_ordered.Biomass_ratio>1.0000005]
print(act_full_bio_ordered.loc[biomass_ratio_1]['FD_initial'].min())
biomass_ratio_1_5 = act_full_bio_ordered.index[act_full_bio_ordered.Biomass_ratio>1.51000005]
print(act_full_bio_ordered.loc[biomass_ratio_1_5]['FD_initial'].min())
#%%
## PREDICTING CHANGE IN FUNCTIONAL DIVERSITY
fd_fd_doca = smf.ols('FD_ratio ~ FD_initial + DOC_initial_int', lognona)
fd_fd_doca = fd_fd_doca.fit()
print(fd_fd_doca.summary())
