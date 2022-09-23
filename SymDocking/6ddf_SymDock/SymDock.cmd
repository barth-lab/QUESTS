PATH_TO_ROSETTA/SymDock.linuxgccrelease \
    -database PATH_TO_ROSETTA/trunk_072512/rosetta_database/ \
    -in:file:s input/"$1"_INPUT.pdb \
    -in:file:native input/6DDF_c2exRpkAA_AB_3oduSymNative.pdb \
    -symmetry:symmetry_definition input/$2 \
    -out:nstruct 200 \
    -out:file:silent output/"$1"_symdock.out \
    -out:file:scorefile output/score_$1.sc \
    @input/Native_SymDock_flag \
    -out:suffix _$$ \
    -seed_offset $$ \
    > output/"$1"_$$.log \
    2> output/"$1"_$$.err
