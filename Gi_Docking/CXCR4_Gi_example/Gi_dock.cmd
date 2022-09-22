PATH_TO_ROSETTA/rosetta_source/bin/docking_protocol.static.linuxgccrelease \
	-in:file:s input/cxcr4_Gi.pdb \ 
        @input/Dock_flag \
	-out:suffix _$$ \
        -out:file:fullatom \
        -out:file:silent_struct_type binary \
	-out:file:scorefile score_$$.sc \
	-seed_offset $$ \
	> $$.log \
	2> $$.err
