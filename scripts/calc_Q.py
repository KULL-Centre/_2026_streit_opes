import MDAnalysis as mda
from MDAnalysis.lib.util import convert_aa_code
import numpy as np
import argparse
import mdtraj as md

# Warning - this script assumes that coordinates of the input files are in Angstrom
# Input contacts should be 1 indexed (i.e., the same numbering as in the pdb file)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculate fraction of native contacts using PLUMED Q switching function. Assumes coordinates are in Angstrom.")
    parser.add_argument("-p", "--pdb", required=True, help="Reference PDB file")
    parser.add_argument("-x", "--xtc", required=False, help="Trajectory XTC file")
    parser.add_argument("-c", "--contacts", required=True, help="Atom numbers matching the PDB file")
    parser.add_argument("--r0", required=True, help="Reference distances file")
    parser.add_argument("--lambda_param", type=float, default=1.8, help="Lambda parameter for Q switching function")
    parser.add_argument("--beta", type=float, default=50.0, help="Beta parameter (nm^-1) for Q switching function")
    parser.add_argument("-o", "--output", default="native_contacts.txt", help="Output text file")
    return parser.parse_args()

def q_switching(rij, r0, lambda_param, beta):
    z = beta * (rij - lambda_param * r0)
    # stable sigmoid avoiding overflow completely
    s = np.zeros_like(z)
    pos = z >= 0
    neg = ~pos
    # large positive z -> s ~ 0
    s[pos] = np.exp(-z[pos]) / (1.0 + np.exp(-z[pos]))
    # large negative z -> s ~ 1
    s[neg] = 1.0 / (1.0 + np.exp(z[neg]))
    return s

def main():
    args = parse_arguments()

    r0_values = np.loadtxt(args.r0) / 10 # Å -> nm

    # load contacts information (residue numbers, names, atom names)
    with open(args.contacts, 'r') as f:
        lines = f.readlines()

    # Load Universe
    try:
        u = mda.Universe(args.pdb, args.xtc)
    except:
        u = mda.Universe(args.pdb)

    # find atom indices for contacts
    native_contacts = []
    for line in lines:
        res1 = int(line.strip().split()[0])
        res2 = int(line.strip().split()[3])
        resname1 = str(line.strip().split()[1])
        resname2 = str(line.strip().split()[4])
        atomname1 = str(line.strip().split()[2])
        atomname2 = str(line.strip().split()[5])
        #print(res1, resname1, atomname1, res2, resname2, atomname2)
        # ensure ILE CD vs CD1 handled correctly
        if resname1=='ILE' and 'CD' in atomname1:
            atomname1='CD*'
        elif resname2=='ILE' and 'CD' in atomname2:
            atomname2='CD*'
        # below makes it compatible with different names for e.g. histidine while still checking residue identity
        at1 = u.select_atoms(f"resid {res1} and resname {resname1[:2]}* and name {atomname1}").indices[0]
        at2 = u.select_atoms(f"resid {res2} and resname {resname2[:2]}* and name {atomname2}").indices[0]
        native_contacts.append([at1, at2])
    native_contacts = np.array(native_contacts)

    all_atoms = u.select_atoms("all")

    # Compute fraction of native contacts for each frame
    fractions = []
    for ts in u.trajectory:
        frame_positions = all_atoms.positions / 10.0  # Å -> nm
        rij = np.linalg.norm(
            frame_positions[native_contacts[:, 0]] -
            frame_positions[native_contacts[:, 1]], axis=1
        )
        sij = q_switching(rij, r0_values, args.lambda_param, args.beta)
        fractions.append(np.mean(sij))

    fractions = np.array(fractions)

    # Save
    np.savetxt(args.output, fractions, fmt="%.6f")

    print(f"Fraction of native contacts saved to {args.output}")

if __name__ == "__main__":
    main()

