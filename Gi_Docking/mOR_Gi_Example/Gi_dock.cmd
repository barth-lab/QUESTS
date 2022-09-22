PATH_TO_ROSETTA/rosetta_source/bin/docking_protocol.linuxgccrelease \
    -database /home/jefferso/sideloads/barth/Source/Minirosetta/trunk_072512/rosetta_database/ \
    -ex1 \
    -ex2 \
    -ex2aro \
    -dock_pert 3 8 \
    -out:file:fullatom \
    -out:file:silent_struct_type binary \
    -in:file:s input/$1.pdb  \
    -out:file:silent output.$N/"$1"_dock.out \
    -out:file:scorefile output.$N/score_$1.sc \
    -out:suffix _$s \
    -out:nstruct 250 \
    -seed_offset $$ \
    > output.$N/log_$s\
    2> output.$N/err_$s
