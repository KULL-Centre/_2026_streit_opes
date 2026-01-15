import mdtraj as md
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Process MD simulation input arguments.")

parser.add_argument("--pdb", type=str, help="Path to the PDB file.")
parser.add_argument("--traj", type=str, help="Path to the trajectory file.")
parser.add_argument("--save_dir", type=str, default="outputs", help="Directory to save outputs. Default: outputs")
parser.add_argument("--stride", type=int, default=1, help="Stride for analysing trajectory frames. Default: 1.")

args = parser.parse_args()

pdb = args.pdb
traj = args.traj
save_dir = args.save_dir
stride = args.stride

def compute_radius_of_gyration(traj):
    """Compute radius of gyration for each frame."""
    mass=[]
    for at in traj.topology.atoms:
        mass.append(at.element.mass)
    return md.compute_rg(traj, masses=np.array(mass))

trajectory = md.load(traj, top=pdb)

if stride>1:
    trajectory = trajectory[::stride]

# only select protein atoms
protein_idxs = trajectory.top.select('protein or resname NH2 or resname ACE or resname NME')
subtraj = trajectory.atom_slice(protein_idxs)

rg = compute_radius_of_gyration(subtraj)
np.savetxt(save_dir+'rg_{}.txt'.format(traj[:-4]), rg)

