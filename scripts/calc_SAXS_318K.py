import os
import mdtraj as md
import numpy as np
import pandas as pd
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from io import StringIO
import sys

scripts_dir = '/projects/lindorff_group/people/thd711/software/scripts/'
sys.path.append(scripts_dir)
import pepsi_saxs 

def find_fasta_filenames(path_to_dir, suffix=".fasta"):
    filenames = os.listdir(path_to_dir)
    return [filename for filename in filenames if filename.endswith(suffix)]

def find_traj_filenames(path_to_dir, suffix=".xtc"):
    filenames = os.listdir(path_to_dir)
    return [filename for filename in filenames if filename.endswith(suffix)]

# === CONFIG ===
proteins = ['ACTR_REST2_300K_ensemble']
input_root = './'
exp_root = 'exp_data_ACTR/'
output_basename = 'pepsisaxs_318K'
pepsi_exe = '/projects/lindorff_group/people/thd711/software/pepsi/Pepsi-SAXS'

# === Load protein list ===
SAXSproteins = []
for protein in proteins:
    if os.path.exists(exp_root+f"SAXS_bift_318K.dat"):
        SAXSproteins.append(protein)
print(f"SAXS Proteins: {SAXSproteins}")

# === Load sequences ===
seq_dict = {}
for prot in SAXSproteins:
    with open(exp_root+f'seq.fasta','r') as f:
        seq = str(f.readlines()[1].strip())
        seq_dict[prot]=seq

# === Iterate over proteins ===

for protein in SAXSproteins:

    print(f"Processing {protein}...")

    traj_files = find_traj_filenames(input_root+f"{protein}/")

    for traj in traj_files:
        traj = traj[:-4]

        if os.path.exists(input_root+f"{protein}/"+f"{output_basename}/saxs_curves_{traj}.csv"):
            print(f'Existing saxs calculations for {protein} file {traj} and skipping...')
            continue
        else:
            os.makedirs(input_root+f"{protein}/"+f"{output_basename}", exist_ok = True)
            if not os.path.exists(input_root+f"{protein}/protein.pdb"):
                print(f'Missing {protein} and skipping...')
                continue

            pepsi_saxs.run_pepsi(trajectory=input_root+f"{protein}/{traj}.xtc", topology=input_root+f"{protein}/topology.pdb", saxs=exp_root+f'SAXS_bift_318K.dat', 
            output = input_root+f"{protein}/"+f"{output_basename}", pepsi = pepsi_exe, pH = None, sequence = seq_dict[protein], keep_tmp = False)

            # rename saxs_curves.csv
            # Current file name
            old_name = input_root+f"{protein}/"+f"{output_basename}/saxs_curves.csv"
            # New file name
            new_name = input_root+f"{protein}/"+f"{output_basename}/saxs_curves_{traj}.csv"
            # Rename the file
            os.rename(old_name, new_name)
                    
print('All done!')
