import os
import shutil
import MDAnalysis
import numpy as np
from DEERPREdict.PRE import PREpredict
import pandas as pd
import pickle as pkl
import gzip

def find_dat_filenames(path_to_dir, suffix=".dat"):
    filenames = os.listdir(path_to_dir)
    return [filename for filename in filenames if filename.endswith(suffix)]

def find_traj_filenames(path_to_dir, suffix=".xtc"):
    filenames = os.listdir(path_to_dir)
    return [filename for filename in filenames if filename.endswith(suffix)]

# === CONFIG ===
proteins = ['ACTR_REST2_300K_ensemble']
input_root = './'  
exp_root = 'exp_data_ACTR/'
output_basename = 'PRE'
libname = 'MTSSL MMMx' # PRE rotamer library

# === Load protein info ===

PREproteins_info = {}
for protein in proteins:
    datfiles = find_dat_filenames(exp_root)
    info = pd.read_csv(exp_root+f'info.csv')
    PRE = 'no'
    for file in datfiles:
        if 'PRE' in file:
            PRE='yes'
            continue
    if PRE=='yes':
        

        labelling_sites = []
        temps = []
        for file in datfiles:
            if 'PRE' in file:
                res = int(file.split('.')[0].split('-')[-1])
                dataset = file[:-4]
                temp = np.mean(info[info['Experiment']==dataset]['Temp(K)'])
                labelling_sites.append(res)
                temps.append(temp)

        PREproteins_info[protein] = {}
        PREproteins_info[protein]['Temperatures'] = temps
        PREproteins_info[protein]['Sites'] = labelling_sites
 
print(f"PRE proteins: {proteins}")

# === Iterate over proteins ===

for protein in proteins:

    trajnames = find_traj_filenames(input_root+f"{protein}/")

    for traj in trajnames:

        trajbase = traj[:-4]

        # iterate through all labelling sites
        for site, temp in zip(PREproteins_info[protein]['Sites'], PREproteins_info[protein]['Temperatures']):

            print(f"Processing {protein} trajectory file {traj} PRE-label {site}...")

            if os.path.exists(input_root+f"{protein}/{output_basename}/PRE-{trajbase}-{site}.pkl"):
                print(f'Existing PRE calculations for {protein}/{traj} and skipping...')
                continue
            else:
                os.makedirs(input_root+f"{protein}/{output_basename}", exist_ok=True)
                # load trajectory
                u = MDAnalysis.Universe(input_root+f"{protein}/topology.pdb", input_root+f"{protein}/{trajbase}.xtc")
                
                # check name of amide hydrogen
                C1 = u.select_atoms('name HN').indices.shape[0]
                C2 = u.select_atoms('name H').indices.shape[0]
                if C1>0:
                    amide_name = 'HN'
                else:
                    amide_name = 'H'

                PRE = PREpredict(u, residue = site, libname = libname,
                                    tau_t = .5*1e-9, log_file = input_root+f"{protein}/{output_basename}/log", temperature = temp, z_cutoff = 0.05, atom_selection = amide_name, Cbeta = False)
                PRE.run(output_prefix = input_root+f"{protein}/{output_basename}/PRE-{trajbase}", tau_t = .5*1e-9,delay = 10e-3,
                        tau_c = 5*1e-09,r_2 = 10, wh = 750) # r_2 and wh don't matter yet here
            
print('All done!')
