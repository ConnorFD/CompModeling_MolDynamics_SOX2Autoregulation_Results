import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns

# =========================
# INPUT FILES
# =========================
TPR = "md_0_100.tpr"
XTC = "md_0_20.xtc"

u = mda.Universe(TPR, XTC)

protein = u.select_atoms("protein")
rna = u.select_atoms("nucleic")

print("Protein atoms:", len(protein))
print("RNA atoms:", len(rna))

# =========================
# H-BOND SETUP
# =========================
donors = protein.select_atoms("name N* or name O* or name S*")
acceptors = rna.select_atoms("name N* or name O* or name P*")

cutoff_dist = 3.5
cutoff_angle = 120.0

n_frames = 10000

hb_frames = {}

print("Running hydrogen bond analysis...")

# =========================
# MAIN LOOP (FIRST 5000 FRAMES)
# =========================
for i, ts in enumerate(tqdm(u.trajectory, desc="Frames")):

    if i >= n_frames:
        break

    d_pos = donors.positions
    a_pos = acceptors.positions

    dist = distance_array(d_pos, a_pos, box=ts.dimensions)

    pairs = np.where(dist < cutoff_dist)

    frame_hits = set()

    for di, aj in zip(*pairs):

        D = donors[di]
        A = acceptors[aj]

        vec = A.position - D.position
        norm = np.linalg.norm(vec)

        if norm == 0:
            continue

        vec /= norm

        H = D.position + vec

        v1 = D.position - H
        v2 = A.position - H

        cosang = np.dot(v1, v2) / (
            np.linalg.norm(v1) * np.linalg.norm(v2)
        )

        angle = np.degrees(np.arccos(np.clip(cosang, -1, 1)))

        if angle >= cutoff_angle:

            key = (
                D.resname, D.resid,
                A.resname, A.resid
            )

            frame_hits.add(key)

    for key in frame_hits:
        if key not in hb_frames:
            hb_frames[key] = set()
        hb_frames[key].add(i)

# =========================
# OCCUPANCY CALCULATION
# =========================
rows = []

for key, frames in hb_frames.items():
    occupancy = len(frames) / n_frames * 100

    rows.append([
        key[0], key[1],
        key[2], key[3],
        occupancy
    ])

df = pd.DataFrame(
    rows,
    columns=["ProteinRes", "ProteinID", "RNARes", "RNAID", "Occupancy"]
)

df = df.sort_values("Occupancy", ascending=False)

df.to_csv("hb_occupancy.csv", index=False)

top15 = df.head(15).copy()

print("\nTop 15 interactions:")
print(top15)

# =========================
# BUILD HEATMAP MATRIX (NO WHITE, FULL 0-FILL)
# =========================

protein_labels = [
    f"{r.ProteinRes}{r.ProteinID}"
    for _, r in top15[["ProteinRes", "ProteinID"]]
    .drop_duplicates()
    .sort_values("ProteinID")
    .iterrows()
]

rna_labels = [
    f"{r.RNARes}{r.RNAID}"
    for _, r in top15[["RNARes", "RNAID"]]
    .drop_duplicates()
    .sort_values("RNAID")
    .iterrows()
]

heatmap = pd.DataFrame(
    0.0,
    index=protein_labels,
    columns=rna_labels
)

for _, row in df.iterrows():

    p = f"{row['ProteinRes']}{row['ProteinID']}"
    r = f"{row['RNARes']}{row['RNAID']}"

    if p in protein_labels and r in rna_labels:
        heatmap.loc[p, r] = row["Occupancy"]

# =========================
# HORIZONTAL HEATMAP PLOT
# =========================
plt.figure(figsize=(10, 6))

ax = sns.heatmap(
    heatmap.T,              # ?? FLIPPED (horizontal layout)
    cmap="viridis",
    linewidths=1.5,
    linecolor="white",
    annot=False,
    vmin=0,
    square=True,
    cbar_kws={"label": "Occupancy (%)"}
)

ax.set_title(
    "Top 15 Protein-RNA Hydrogen Bond Occupancies",
    fontsize=16
)

ax.set_xlabel("Protein residue")
ax.set_ylabel("RNA residue")

plt.xticks(rotation=45, ha="right")
plt.yticks(rotation=0)

plt.tight_layout()

plt.savefig("top15_hbond_heatmap.png", dpi=300)
plt.close()

print("\nDONE")
print("Saved: hb_occupancy.csv")
print("Saved: top15_hbond_heatmap.png")
