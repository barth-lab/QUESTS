# both commands can generate C2 symmetry 
# make_symdef_file.pl can generate symmetry based on given structure, making comparing with the X-ray structure easier. 
# rosetta_source/src/apps/public/symmetry/make_symmdef_file_denovo.py -symm_type cn -nsub 2 > c2_symmetry
rosetta_source/src/apps/public/symmetry/make_symmdef_file.pl -m NCS -a A -i B  -r 12.0 -p 3oduAB.pdb > c2_symmetry
