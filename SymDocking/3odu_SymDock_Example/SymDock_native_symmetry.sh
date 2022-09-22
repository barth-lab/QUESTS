#!/bin/bash

ROS="PATH_TO_ROSETTA"

$ROS/trunk_072512/rosetta_source/bin/SymDock.static.linuxgccrelease @input/Native_SymDock_flag -out:suffix _$$ -seed_offset $$ > $$.log
