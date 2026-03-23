# CompModeling_MolDynamics_SOX2Autoregulation_Results
mmCIF files of AlphaFold3 models, ChimeraX-predicted hydrogen bond tables, GROMACS RMSD/RMSF trajectories, GROMACS MD input files, force field packages, and PDBs of selected GROMACS trajectory frames analyzed in "Computational Modeling and Molecular Dynamics Reveal SOX2 Autoregulation Through Competitive mRNA Binding with HuR". Further data, including but not limited to full JSON files/confidence metric summaries of AlphaFold3 predictions and MD trajectory files, are available from the authors upon request.

The AlphaFold3_Models folder contains twenty mmCIF files of AlphaFold3 models, as five AlphaFold3 models were created for each of the four protein–RNA interactions analyzed in this manuscript: the SOX2–ES2 interaction, the SOX2–S100A14 interaction, the SOX2–SOX2 mRNA interaction, and the HuR–SOX2 mRNA interaction.

Inside the ChimeraX_HBonds folder, there are twenty tables of hydrogen bonds, which were predicted by ChimeraX on all AlphaFold3 models. Importantly, in the first SOX2–SOX2 mRNA interaction and the fourth HuR–SOX2 mRNA interaction, hydrogen bonds were predicted after the structure was modified to ease steric clashes.

The GROMACS_MD_Input_Files folder encompasses twelve sub-folders, each containing the input files for a specific GROMACS MD simulation and the command-line prompts executed, and one sub-folder with the AMBER 14SB OL15 force field package used in the MD simulations (dihedral angle issues have already been fixed). The twelve MD simulation sub-folders come from the five MD simulations performed for both the SOX2–ES2 and SOX2–S100A14 interactions, as well as the one MD simulation performed for both the SOX2–SOX2 mRNA and HuR–SOX2 mRNA interactions.

Within the GROMACS_RMSD_Trajectories folder, there are twelve files that correspond to each of the twelve sub-folders in the GROMACS_MD_Input_Files folder, as the RMSD tables were derived from the GROMACS MD simulations performed on each interaction.

The GROMACS_RMSF_Trajectories folder also corresponds almost exactly to each of the twelve sub-folders in the GROMACS_MD_Input_Files folder. Its only difference from the GROMACS_RMSD_Trajectories folder is that it contains twenty four files, since RMSF tables were calculated both for protein residues and for RNA residues.

There are only two files in the Selected_Frame_PDBs folder; these are both PDB files that are simply individual frames selected from MD simulations to be specifically analyzed in the paper. One file is the frame at 10 ns from the SOX2–SOX2 mRNA MD simulation, and one file is the frame at 0.97 ns from the first SOX2–S100A14 MD simulation.
