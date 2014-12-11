#!/bin/bash

n=9;
file=pent;
part_type=dual;

##arguments=("l2g" "lcc" "g2l" "bs" "be" "bn" "bw" "bl" "bh" "bp" "su" "partitioning" "recv" "send" "su" "cnorm" "var")
arguments=("cnorm" "var" "direc1" "direc2" "resvec")


cd pos-lab/A2.3/code/output;


ar_length=${#arguments[@]}
for((ar=0; ar<ar_length; ar=$ar+1))
do
	for((i=0; i<n; i=$i+1))
	do
		arg_name=${arguments[ar]}
		echo "$arg_name-$i";
		diff $arg_name.$file.$part_type-oneread.$i.out $arg_name.$file.$part_type-allread.$i.out
	done
done


for((i=0; i<n; i=$i+1))
do
	echo "$file-$i";
	diff $file.$part_type-oneread.$i.out $file.$part_type-allread.$i.out
done
