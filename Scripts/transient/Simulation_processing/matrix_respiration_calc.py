#%%
## Import libraries
import os
import numpy as np
import sys
import h5py 


## LOAD RESULTS
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient", "activity_loss_-02")
results_dir = os.path.join(project_dir, "results")
filestring = sys.argv[1] + '_carbon_' #null
seed_sim_list = [610229235, 983307757, 643338060, 714504443, 277077803, 898393994, 420,13012022,13061989]

cn_list = [3,6,12,18]
bio_n_series = [4,8,12,16,20,24,28,32]
ip = 0
init_dom_list = [1000,2000,5000,10000,15000]
filestring = 'null_carbon_' #null
#%%
for c_n in cn_list[:1]:
    row = []
    c_b_row = []
    results_filename = os.path.join(results_dir, filestring + str(c_n))

    for seed_sim in seed_sim_list[:1]:
        # Load all datasets and save their Shannon and diversity indices in a dataframe
        seed_all = 'seed_'+str(seed_sim)
        
        details_subfolder = filestring + str(c_n) + '_'+str(seed_sim) + '_ip_0'
        simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
        hrw = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'r+')

        for b_n in bio_n_series[:1]:
            for t_dom_initial in init_dom_list[:1]:
                base_case = "b_1_all_"
                c_b = "bio_n_"+ str(b_n)
                dom_init = "dom_initial_" + str(t_dom_initial)
                doc_input = (t_dom_initial)
                #print(c_n, seed_sim, b_n, t_dom_initial)
                if hrw[base_case][c_b][dom_init][seed_all]:
                    sim_data = hrw[base_case][c_b][dom_init][seed_all]
                    sim_seed_details = base_case+"/"+c_b+"/"+dom_init+"/"+seed_all
                    dom_n, bio_n = sim_data['species_number']
                    paras = sim_data['parameters']
                    uptake_coeff = np.asarray(paras['carbon_uptake'])
                    growth_coeff = np.asarray(paras['yield_coefficient'])
                    v_enz = np.asarray(paras['exo_enzyme_rate'])
                    v_params = np.asarray(paras['max_rate'])
                    k_params = np.asarray(paras['half_saturation'])
                    max_cap = 0.15*t_dom_initial
                    x = sim_data['solution']
                    C = np.asarray(x['dom'])[:10,:]
                    B = np.asarray(x['biomass'])[:10,:]
                    params_mult = uptake_coeff*v_params*v_enz
                    #print(np.shape(params_mult))
                    #mm_c = C[0,:][...,None]/(k_params + C[0,:][...,None])
                    B_total = np.sum(B, axis = 1)
                    #R = np.zeros((np.shape(B)[0],))
                    mm_num = np.matmul(C,params_mult)
                    print(np.shape(mm_num))
                    #mm_denom = 
                    freg = growth_coeff*(1 - (np.sum(B,axis=0) - B)/max_cap)
                    B_freg = B*freg
                    G = np.sum(mm_num*B_freg, axis = 1)
                    print(np.shape(G))
                    #for t in list(range(np.shape(B)[0])):
                    #    c_uptake = params_mult*B[t,:]*C[t,:][...,None]/(k_params + C[t,:][...,None])
                    #    B_total = np.sum(B[t,:])
                    #    R_sum = 0
                    #    for i in list(range(bio_n)):
                    #        #print(i)
                    #        R_sum += np.sum(c_uptake[:,i]*(1 - growth_coeff[:,i]*(1 - ((B_total - B[t,i])/max_cap))))
                    #        #print(R_sum)
                    #    R[t] = R_sum
                    #dataset_name = base_case+"/"+c_b+"/"+dom_init+"/"+seed_all+"/solution/respiration"
                    #hrw.create_dataset(dataset_name, data=R)
        hrw.close()