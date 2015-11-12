#!/usr/bin/env bash

# Combine rij.txt files for each complex
# $1 = number of complexes
# $2 = pool (native or nonnative)

stop=$1
pool=$2

# First complex
if [ "$pool" == "native" ]; then
	echo "native"
	start=0
	mv r12ij.0.txt r12_tmp1
	mv r10ij.0.txt r10_tmp1
	mv r6ij.0.txt r6_tmp1
else	
	echo "nonnative"
	start=1
	mv r12ij.1.txt r12_tmp1
	mv r10ij.1.txt r10_tmp1
	mv r6ij.1.txt r6_tmp1
fi

echo "$start"

# All other complexes
# Increment $start by one for evaluating the next pose
let start=start+1
for ((i=$start;i<=$stop;i++)); do
    
    echo $i
    
    awk '{print $2}' r12ij.$i.txt > r12_tmp2
    paste r12_tmp1 r12_tmp2 > r12_tmp3
    mv r12_tmp3 r12_tmp1
    rm r12ij.$i.txt r12_tmp2 

    awk '{print $2}' r10ij.$i.txt > r10_tmp2
    paste r10_tmp1 r10_tmp2 > r10_tmp3
    mv r10_tmp3 r10_tmp1
    rm r10ij.$i.txt r10_tmp2

    awk '{print $2}' r6ij.$i.txt > r6_tmp2
    paste r6_tmp1 r6_tmp2 > r6_tmp3
    mv r6_tmp3 r6_tmp1
    rm r6ij.$i.txt r6_tmp2

done

# Rename
mv r12_tmp1 r12ij.txt
mv r10_tmp1 r10ij.txt
mv r6_tmp1 r6ij.txt
