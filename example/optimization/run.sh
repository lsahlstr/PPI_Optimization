#!/usr/bin/env bash

EXEC=opt_eps_rmin.R
EXECDIR=/export/nVerde/users/lsahlstr/repos/PPI_Optimization/bin
cp -p $EXECDIR/*.R .

# Define command line arguments for optimization
npar=20
eps_file=eii_lj.txt
rmin_file=rii_wSD.txt
mask_file=mask_20.txt
list=pdblist_lj
potential=lj
fitness=zscore
popsize=10
ncycle=100
min=0.005
max=0.3

# Run optimization
module load R
./opt_eps_rmin.R \
    -e $eps_file \
    -r $rmin_file \
    -m $mask_file \
    -p $potential \
    -l $list \
    -f $fitness \
    -c $ncycle \
    -s $popsize \
    -u $max \
    -o $min


# Plot the data
gnuplot plot_score.gnu

for i in {1..6}; do
    awk '$2 == "'$i'" { print $0 }' rmsdEnerComp.txt > sys${i}.dat
    gnuplot -e "sys='$i'" plot_new_sys.gnu
    gnuplot -e "sys='$i'" plot_old_sys.gnu
done

