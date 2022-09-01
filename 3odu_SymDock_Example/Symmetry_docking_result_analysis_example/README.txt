Steps to analyze the symmetry docking results.
Please install hdbscan, numpy, matplotlib and seaborn (for plotting) 

1. Select the best 10% from the symmetry docking results by interface energy. In the test case, top 1000 by interface energy from the 10,000 decoys are selected.
   Extract the corresponding PDBs from the silent out file. 
   Keep the 'I_sc' and 'description' information from the score file to 'lowest_score.txt'.

2. Generate cross-angle and distance parameters for the generated decoys. Filter out the dimer structures with too large cross-angles (directly measured by cross-angle) and only extracellular/intracellular-loop-interfaces (featured by large distance). This should be a case-by-case parameters that need to be tuned for each GPCR dimer. In the test case, 744 decoys are selected for further analysis.
   Under the folder of decoys, run:

   mkdir ../selected_structure
   bash ../../Code/cross_angle.sh 

   the output will be stored in cross_angle.infor and the selected ones are stored in selected.infor

   cat selected.infor | awk {'print $1, $4, $6'}  > selected.infor.table

3. Then further select top 10% for hdbscan clustering. Please tune the parameters for loose or strigent cutoff based on the clustering results. 
   ../../Code/cluster_hdbscan.py selected.infor.table 10
   ( 10 means select 1/10 of the total decoys by best interface energy )
