# Command lines for GPCR dimer analysis


## Symmetrical Docking
An example with inputs using the crystal structure of inactive-state CXCR4 (PDB ID: 3odu) and code for analysis is in the '3odu_SymDock_Example' folder. The symmetry definition file is present in the input folder, but can otherwise be made with the 'makes_symmetry_file.command' in the input folder

Change 'PATH_TO_ROSETTA' in 'SymDock_native_symmetry.sh' and 'Native_SymDock_flag'. Run it to generate ~10,000 decoys.

Navigate to the 'SymDock_analysis_example' folder and follow the steps in 'SymDock_Analysis.txt' to analyze the symmetric docking results.


## bArr Docking
An example is in the 'bArr_Docking/CXCR4_Arrestin_Example' folder. Inputs can be prepped by aligning receptor models to the 4ZWJA_receptor.pdb complex in the 'prepare' folder and using the 8 arrestin structures (ar1-8.pdb). Change the 'PATH_TO_ROSETTA' in the 'bArr_dock.cmd' and in 'input/dock_flag'. Generate ~1000 decoys for each arrestin model to get ~8000 total decoys. Select the top 1% by interface score and cluster for native-like finger loop position.

## Gi Docking
examples are in the 'Gi_Docking/mOR_Gi_Example' folder and 'Gi_Docking/cxcr4_Gi_Example'. Change the 'PATH_TO_ROSETTA' in the 'Gi_dock.cmd' and generate ~10,000 decoys. Select the top interface scoring model with native-like G protein C-terminal helix position among the top 1% interface scoring decoys.

## AlphaFold Modeling

Sequences in Fasta folder were used as input for ColabFold (v1.4).
