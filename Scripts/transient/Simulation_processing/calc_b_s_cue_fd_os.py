## Import libraries
import os
import pandas as pd
import numpy as np
import sys
import h5py 

def calc_chars(data, timidx):
    x = data['solution']
    C = np.asarray(x['dom'])
    B = np.asarray(x['biomass'])
    NOSC = np.asarray(x['carbon_os'])
    CUE = np.asarray(x['cue'])
    FD = np.asarray(x['func_div'])
    S = np.asarray(x['shannon'])
    Biomass = np.sum(B, axis = 1)
    DOC = np.sum(C, axis = 1)
    si = S[0]
    sf = S[-1]
    doci = DOC[0]
    docf = DOC[-1]
    bi = Biomass[0]
    bf = Biomass[-1]
    #Maximum Biomass
    bbm = np.max(Biomass)
    #Others at max biomass
    timbm = np.argmax(Biomass)
    sbm = S[timbm]
    docbm = DOC[timbm]
    #CUE
    cuei = CUE[0]
    cuef = CUE[-1]
    cuebm = CUE[timbm]
    #Functional diversity
    fdi = FD[0]
    fdf = FD[-1]
    fdbm = FD[timbm]
    #Oxidation state
    osi = NOSC[0]
    osf = NOSC[-1]
    osbm = NOSC[timbm]
    initial = [doci, bi, si, cuei, fdi, osi]
    final = [docf, bf, sf, cuef, fdf, osf]
    maxbio = [docbm, bbm, sbm, cuebm, fdbm, osbm]
    if timidx >0:
        docc, bc, sc, cuec, fdc, osc = DOC[timidx], Biomass[timidx], S[timidx], CUE[timidx], FD[timidx], NOSC[timidx]
    else:
        docc, bc, sc, cuec, fdc, osc = "NA","NA","NA","NA","NA","NA"
    chardata = [docc, bc, sc, cuec, fdc, osc]

    return initial, final, maxbio, chardata

   
## LOAD RESULTS
#project_dir = os.path.join("C:/","Users","swkh9804","Documents","Project_data","HoliSoils", "transient", sys.argv[1])
project_dir = os.path.join("D:/", "Projects", "HoliSoils","data","transient", sys.argv[1])
results_dir = os.path.join(project_dir, "results")
filestring = 'competition_adaptation_carbon_' #null
seed_sim_list = [610229235, 983307757, 643338060, 714504443, 277077803, 898393994, 420,13012022,13061989]

cn_list = [12,18]#[3,6,12,18]
bio_n_series = [4,8,12,16,20,24,28,32]
ip = 0
init_dom_list = [1000,2000,5000,10000,15000]
transient_switch = 0
input_factor = transient_switch*5/365
loss_crit = 0.9#sys.argv[2]#0.37#0.63
result_fstring = "_loss_"+str(loss_crit)
tim_file = os.path.join(results_dir, 'competition_adaptation_carbon__loss_0.9_combined_dataset.pkl')
tim_data = pd.read_pickle(tim_file)

for c_n in cn_list:
    row = []
    c_b_row = []
    results_filename = os.path.join(results_dir, filestring + str(c_n))

    for seed_sim in seed_sim_list:
        # Load all datasets and save their Shannon and diversity indices in a dataframe
        seed_all = 'seed_'+str(seed_sim)
        
        details_subfolder = filestring + str(c_n) + '_'+str(seed_sim) + '_ip_' + str(ip)
        simulations_dir = os.path.join(project_dir, "simulations", details_subfolder)
        hr = h5py.File(os.path.join(simulations_dir,"simulations.h5"), mode = 'r')

        filename = os.path.join(simulations_dir, "seeds_randoms.pkl")
        seed_details = pd.read_pickle(filename)

        for b_n in bio_n_series:
            tim_subset = tim_data[(tim_data['carbon_species']==c_n)&(tim_data['Seed']==seed_sim)&(tim_data['biomass_species']==b_n)].reset_index()
            for t_dom_initial in init_dom_list:
                c_b = "bio_n_"+ str(b_n)
                dom_init = "dom_initial_" + str(t_dom_initial)
                doc_input = (t_dom_initial) * input_factor
                for baseline in ["b_1", "b_2", "b_3", "b_4","b_5"]:
                    if baseline == "b_1":
                        label = "_all_"
                        sim = baseline + "_all_"
                        char_tim = tim_subset[tim_subset.Sim_series==sim]['T_50'].values[0].astype(int)
                        sim_data = hr[sim][c_b][dom_init][seed_all]
                        init, fin, atmaxbio, chartimsc = calc_chars(sim_data,char_tim)
                        doci, docf, docbm, docc = init[0], fin[0], atmaxbio[0], chartimsc[0]
                        bi, bf, bbm, bc = init[1], fin[1], atmaxbio[1], chartimsc[1]
                        si, sf, sbm, sc = init[2], fin[2], atmaxbio[2], chartimsc[2]
                        cuei, cuef, cuebm, cuec = init[3], fin[3], atmaxbio[3], chartimsc[3]
                        fdi, fdf, fdbm, fdc = init[4], fin[4], atmaxbio[4], chartimsc[4]
                        osi, osf, osbm, osc = init[5], fin[5], atmaxbio[5], chartimsc[5]
                        row.append([seed_sim,sim, c_n, b_n, doci, docc, docbm, docf,  si, sc, sbm, sf,  bi, bc, bbm, bf, cuei, cuec, cuebm, cuef, fdi, fdc, fdbm, fdf, osi, osc, osbm, osf])
                    else:
                        for label in ["a", "b", "c","d","e"]:
                            sim = baseline + "_" + label + "_"
                            char_tim = tim_subset[tim_subset.Sim_series==sim]['T_50'].values[0].astype(int)
                            sim_data = hr[sim][c_b][dom_init][seed_all]
                            init, fin, atmaxbio, chartimsc = calc_chars(sim_data,char_tim)
                            bi, bf, bbm, bc = init[1], fin[1], atmaxbio[1], chartimsc[1]
                            si, sf, sbm, sc = init[2], fin[2], atmaxbio[2], chartimsc[2]
                            cuei, cuef, cuebm, cuec = init[3], fin[3], atmaxbio[3], chartimsc[3]
                            fdi, fdf, fdbm, fdc = init[4], fin[4], atmaxbio[4], chartimsc[4]
                            osi, osf, osbm, osc = init[5], fin[5], atmaxbio[5], chartimsc[5]
                            row.append([seed_sim,sim, c_n, b_n, doci, docc, docbm, docf,  si, sc, sbm, sf,  bi, bc, bbm, bf, cuei, cuec, cuebm, cuef, fdi, fdc, fdbm, fdf, osi, osc, osbm, osf])
        hr.close()
    diversity_data = pd.DataFrame.from_records(row, columns = ["Seed", "Sim_series", "carbon_species", "biomass_species", "DOC_initial", "DOC_timscale", "DOC_maxbio", "DOC_final", "S_initial", "S_timscale", "S_maxbio", "S_final","Biomass_initial", "Biomass_timscale", "Biomass_maxbio", "Biomass_final","CUE_initial", "CUE_timscale", "CUE_maxbio", "CUE_final","FD_initial", "FD_timscale", "FD_maxbio", "FD_final","NOSC_initial", "NOSC_timscale", "NOSC_maxbio", "NOSC_final"])
    print("The shape of the dataframe is ", diversity_data.shape)
    print("The dataframe contains the following data types ", diversity_data.dtypes)

    filename = os.path.join(results_dir, results_filename+result_fstring+"_cue_data.pkl")
    diversity_data.to_pickle(filename)
    print ("CUE data is saved here ", filename)