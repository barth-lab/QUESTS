# Command lines for mOR dimer analysis


## Symmetrical Docking
mOR active-state  (6DDF) and inactive-state (4DKL) aligned to 3ODU CXCR4 symmetry for symmetrical docking
W5.34A mutation made with local flexible backbone code

make symdef file:
```
$ROS/source/src/apps/public/symmetry/make_symmdef_file.pl -m NCS -a A -i B  -r 12.0 -p 6ddf_AB.pdb > c2_symmetry
```
run Symdock simulation:
```
for i in {1..50}; do
	printf -v s “%03d” $i
	$ROS/trunk_072512/SymDock.static.linuxgccrelease \
	-database trunk_072512/rosetta_database/ \
	-in:file:s input/6DDF_c2exRpkAA_AB_3oduSym.10Axsep_INPUT.pdb \
	-in:file:native input/6DDF_c2exRpkAA_AB_3oduSymNative.pdb \
	-symmetry:symmetry_definition input/c2_symmetry \
	-out:nstruct 200 \
	-out:file:silent output/"$1"_symdock.out \
	-out:file:scorefile output/score_$1.sc \
	@input/Native_SymDock_flag \
	-out:suffix _$s \
	-seed_offset $s \
	> output/"$1"_$s.log \
	2> output/"$1"_$s.err"
```

Native_SymDock_flag:
	-in:file:spanfile input/uOR.6DDF.span
	-membrane:fixed_membrane
	-membrane:Menv_penalties
	-membrane:Membed_init
	-score:weights membrane_highres_Menv_smooth_dock3.wts
	-symmetry:membrane
	-symmetry:highres_scorefxn docking_membrane_highres_Menv_smooth_dock3.wts
	-symmetry:highres_pack_scorefxn membrane_highres_Menv_smooth_dock3.wts
	-ignore_unrecognized_res
	-symmetry:initialize_rigid_body_dofs
	-symmetry:symmetric_rmsd
	-packing:ex1
	-packing:ex2
	-packing:ex2aro
	-out:file:silent_struct_type binary
	-out:file:fullatom
	-mute core.util
	-unmute core.init

Top 10% by interface score selected for clustering
Top 20% of best 1000 clustered by angle/distance of TM5
Best interface scoring model in each cluster taken as representative model for mutational effects
Major clusters can be further clustered into subclusters for docking structural models of effectors.

Steps to analyze the symmetry docking results.
Install hdbscan, numpy, matplotlib and seaborn (for plotting) 

1. Select the best 10% from the symmetry docking results by interface energy. In the test case, top 1000 by interface energy from the 10,000 decoys are selected.
	Extract the corresponding PDBs from the silent out file. 
	Keep the 'I_sc' and 'description' information from the score file to 'lowest_score.txt'.
```
	grep SCORE: symdock.out |sort -nk20 |head -1000| awk '{print $20"\t"$NF}' > top1k_Isc.sc
	grep SCORE: symdock.out |sort -nk20 |head -1000| awk '{print $NF}' > top1k_Isc.tags
	$ROS/source/bin/extract_pdbs.linuxgccrelease -in:file:silent 6DDFrlx.WT_symdock.out -in:file:silent_struct_type binary -in:file:tagfile top1k_Isc.tags -out:pdb -silent_read_through_errors
	mkdir decoys_top1k_Isc
	mv *pdb decoys_top1k_Isc/
	cp top1k_Isc.sc decoys_top1k_Isc/lowest_score.txt
```

2. Generate cross-angle and distance parameters for the generated decoys. Filter out the dimer structures with either too large cross-angles (directly measured by cross-angle) or only extracellular/intracellular-loop-interfaces (featured by large distance). This should be a case-by-case parameters that need to be tuned for each GPCR dimer. In the test case, 744 decoys are selected for further analysis. The command line inputs of cross_angle.sh define the TM segments that will be analyzed, so adjust for your specific use case.
	Under the folder of decoys, run:
```
	mkdir ../selected_structure
	bash ../../Code/cross_angle.sh
```
	the output will be stored in cross_angle.infor and the selected ones are stored in selected.infor
```
	cat selected.infor | awk {'print $1, $4, $6'}  > selected.infor.table
```

3. Then further select top 10% for hdbscan clustering. Please tune the parameters for loose or strigent cutoff based on the clustering results.

```
	../../Code/cluster_hdbscan.py selected.infor.table 10
```

( 10 means select 1/10 of the total decoys by best interface energy )

## bArr Docking
	1000 decoys for each arrestin model (8)
	8000 total decoys filtered down to top 1% by interface score
	best interface scoring model with native-like finger loop position (4zwj) selected

command:
```
for i in {1..8}; do
$ROS/trunk_072512/rosetta_source/bin/docking_protocol.linuxgccrelease \
	-database $ROS/trunk_072512/rosetta_database/ \
	-nstruct 1000 \
	-ex1 \
	-ex2 \
	-ex2aro \
	-dock_pert 3 8 \
	-s  mORsymdock_bArr_$i.pdb \
	-out:file:fullatom \
	-out:file:silent_struct_type binary \
	> log \
	2> err
```

## Gi Docking
	10000 decoys 
	best interface scoring model with native-like Gi C-terminal position selected (no clustering, just comparison to native 6ddf)
	clustering was performed on top 1% interface-scoring decoys, then selected for native-like C-term helix position

command:
```
for i in {1..40}; do	$ROS/trunk_072512/rosetta_source/bin/docking_protocol.linuxgccrelease
	-database $ROS/trunk_072512/rosetta_database/
	-ex1 \
	-ex2 \
	-ex2aro \
	-dock_pert 3 8 \
	-out:file:fullatom \
	-out:file:silent_struct_type binary \
	-in:file:s input/$1.pdb  \
	-out:file:silent output/"$1"_dock.out \
	-out:file:scorefile output/score_$1.sc \
	-out:suffix _$s \
	-out:nstruct 250 \
	-seed_offset $$ \
	> output.$N/log_$s\
	2> output.$N/err_$s"
```

