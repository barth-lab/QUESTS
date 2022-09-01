#! /usr/bin/bash

for i in `tail -n +2 lowest_score.txt | awk {'print $1","$2'}`;
do 
  python ../../Code/analyze_dock_result.py $i 180 201 476 497  >> cross_angle.infor
  wait
done
python ../../Code/summary_result.py cross_angle.infor  >> selected.infor

