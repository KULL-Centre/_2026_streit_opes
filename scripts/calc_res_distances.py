import mdtraj as md
import numpy as np

# ==========================
# USER SETTINGS
# ==========================
topology_file = "../topology.pdb"
traj_files = [
    "../traj0.xtc",
    "../traj1.xtc",
    "../traj2.xtc",
    "../traj3.xtc",
    "../traj4.xtc",
]

outdir = './'

# ==========================
# LOAD TOPOLOGY ONCE
# ==========================
top = md.load_topology(topology_file)
n_res = top.n_residues

# generate all residue pairs (i < j)
res_pairs = [(i, j) for i in range(n_res) for j in range(i + 1, n_res)]

# ==========================
# LOOP OVER TRAJECTORIES
# ==========================
for t, trajfile in enumerate(traj_files):
    print(f"Processing trajectory {t}: {trajfile}")

    traj = md.load(trajfile, top=topology_file)

    # compute closest heavy-atom distances (nm)
    distances, pairs = md.compute_contacts(
        traj,
        contacts=res_pairs,
        scheme="closest-heavy"
    )

    # round to 3 decimals to save space
    distances = np.round(distances, 3).astype(np.float16)

    # reshape to (frames, residues, residues)
    n_frames = traj.n_frames
    dist_mat = np.zeros((n_frames, n_res, n_res), dtype=np.float16)

    for k, (i, j) in enumerate(pairs):
        dist_mat[:, i, j] = distances[:, k]
        dist_mat[:, j, i] = distances[:, k]

    # save
    outfile = outdir + f"residue_contacts_traj{t}.npy"
    np.save(outfile, dist_mat)

    print(f"Saved: {outfile}\n")

print("All trajectories processed.")

