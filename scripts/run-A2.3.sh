#!/bin/bash

nprocs=2
files=("pent" "cojack" "drall")
part_types=("dual" "nodal" "classic")

for file_in in "${files[@]}"
do
  for part_type in "${part_types[@]}"
  do
        for((n=2; n<=nprocs; n=$n+1))
        do
              (cd pos-lab/A2.3/code/; mpirun -n $n ./gccg ./data/$file_in.geo.bin $part_type oneread;)
              (cd pos-lab/A2.3/code/; mpirun -n $n ./gccg ./data/$file_in.geo.bin $part_type allread;)
        done
  done
done
