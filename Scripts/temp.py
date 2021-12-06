# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

unsat_sub = biomass_ratio_df[biomass_ratio_df.Sat<1]
unsat_sub["eff_sat"] = unsat_sub.Sat/0.6 - 1/3
marklist = ["o","^","s"]
fig, axes = plt.subplots(1,2, figsize = (7,3))
sns.scatterplot(data = unsat_sub, x = "Time", y = "State_Ratio", hue = "Regime",  hue_order = ["Fast","Slow", "Medium"],
                style = "Regime", palette = my_pal, ax = axes.flat[0])
sns.scatterplot(data = unsat_sub, x = "eff_sat", y = "Loc_Ratio", hue = "Regime", hue_order = ["Fast","Slow", "Medium"],
                style = "Regime", palette = my_pal, ax = axes.flat[1], legend=False)
axes.flat[0].set_xscale("log")
axes.flat[0].set_yscale("log")
axes.flat[0].set_title("Ratio of active and\ninactive biomass", fontsize = 12)
axes.flat[1].set_title("Ratio of immobile and\nmobile biomass", fontsize = 12)
axes.flat[0].set_ylabel("Ratio", fontsize = 12)
axes.flat[1].set_ylabel("")
axes.flat[0].set_xlabel("Residence time\nof solutes (days)", fontsize = 12)
axes.flat[1].set_xlabel("Mean saturation", fontsize = 12)
axes.flat[0].legend(title="Flow regime", fontsize = 11, title_fontsize = 11)
axes.flat[0].set_yticks((0.1,1,10),(0.1,1,10))
axes.flat[0].set_xticks((0.1,1,10,100),(0.1,1,10,100))
axes.flat[1].set_xticks((0.45,0.5,0.55,0.6,0.65),(0.45,0.5,0.55,0.6,0.65))
for a in axes[:]:
    a.tick_params(labelsize = 10)
picname = os.path.join(op_dir,"Fig_4_5_Unsaturated_fractions_microbes.png")
plt.savefig(picname, dpi = 300, bbox_inches = 'tight', pad_inches = 0.1)