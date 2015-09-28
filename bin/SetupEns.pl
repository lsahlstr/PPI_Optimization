#!/usr/bin/env perl

use strict;
use warnings;

my $Nstructs_n = 10;
my $Nstructs_nn = 100;
my $last_a = undef;
my $first_b = undef;
my $n = undef;
my $zdock_dir = "~/CGPPI_FF/decoys_bm4_zd3.0.2_6deg";
my $exec_dir = "~/CGPPI_FF/opt/bin";

my $pdblist = $ARGV[0];
open(PDBS,"$pdblist");
my @pdbs = <PDBS>;

foreach my $pdb (@pdbs) {
    chomp $pdb;
    print "Setting up $pdb\n";
    
    # Make directory for the current system
    `rm -rf $pdb`;
    `mkdir -p $pdb`;
    chdir $pdb;

	# Get the Zdock files for the current system
	my $lig_pdb = "${pdb}_l_u.pdb.ms";
	my $rec_pdb = "${pdb}_r_u.pdb.ms";
	my $zdock_out = "${pdb}.zd3.0.2.fg.out";
	my $zdock_rmsds = "${pdb}.zd3.0.2.fg.out.rmsds";
	`cp -p $zdock_dir/create* .`;
	`cp -p $zdock_dir/input_pdbs/$lig_pdb .`;
	`cp -p $zdock_dir/input_pdbs/$rec_pdb .`;
	`cp -p $zdock_dir/results/$zdock_out .`;
	`cp -p $zdock_dir/results/$zdock_rmsds .`;
	
	# Combine Zdock output and RMSD files
	`head -n 5 $zdock_out > tmp.header`;
	`cat tmp.header $zdock_rmsds > tmp.rmsds`;
	`paste $zdock_out tmp.rmsds | tail -n +6 > tmp.out`;
	
	#### Native-like ensemble ####
	print "Native ensemble\n";
	`mkdir native`;
	`sort -g -k 9 tmp.out | head > tmp.native`;
	`cat tmp.header tmp.native | awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7}' > native.out`;
	`cat tmp.header tmp.native | awk '{print \$8,\$9}' > n_rmsd.dat`;
	# Generate complexes from Zdock output
	`./create.pl native.out`;
	`mv complex*pdb native`;
	
	chdir "native";
	for my $i (1..$Nstructs_n) {
		print "$i ";
		if ($i == 1) {
			# TER line at top of pdb gives the number of ATOM lines in first chain
			# cat -v because ^M formatting
			$n = `cat -v complex.1.pdb | head -n 1 | awk '{print \$2}'`;
			chomp $n;
			$last_a = $n -1;
			$first_b = $n;
			#print "$n\n";
		}
		# Only the ligand PDB is different for each complex. But it is currently unclear if 
		# the first molecule in complex.pdb is always the receptor and the second is always
		# the ligand. So here we construct both molecules from each complex.
		`grep '^ATOM' complex.$i.pdb | head -n $last_a | grep ' CA ' > a.$i.pdb`;
		`grep '^ATOM' complex.$i.pdb | tail -n +$first_b | grep ' CA ' > b.$i.pdb`;
		`$exec_dir/PairDist.pl a.$i.pdb b.$i.pdb $exec_dir/aa_pairs.txt`; # > pairdist.$i.out`;
		`mv r12ij.txt r12ij.$i.txt`;
		`mv r10ij.txt r10ij.$i.txt`;
		`mv r6ij.txt r6ij.$i.txt`;
		`rm a.$i.pdb b.$i.pdb`;
	}
	`$exec_dir/CombineRij.sh $Nstructs_n`;
	`mv r12ij.txt ../n_r12ij.txt`;
	`mv r10ij.txt ../n_r10ij.txt`;
	`mv r6ij.txt ../n_r6ij.txt`;
	
	print "\n";
	chdir "..";
	
	#### Non-native ensemble ####
	print "Non-native ensemble\n";
	`mkdir nonnative`;
	# Near-native (4 A < iRMSD < 8 A)
	`awk ' \$9 > 4 && \$9 < 8 ' tmp.out > tmp.nearn`;
	`sort -R tmp.nearn | head -n 20 > tmp.nearn2`;
	# "Far from"-native (iRMSD >= 8 A)
	`awk ' \$9 >= 8 ' tmp.out > tmp.nn`;
	`sort -R tmp.nn | head -n 80 > tmp.nn2`;
	# Combine near- and far from-native
	`cat tmp.nearn2 tmp.nn2 > tmp.nearn.nn`;
 	`cat tmp.header tmp.nearn.nn | awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7}' > nonnative.out`;
 	`cat tmp.header tmp.nearn.nn | awk '{print \$8,\$9}' > nn_rmsd.dat`;
 	# Generate complexes from Z-dock output
 	`./create.pl nonnative.out`;
 	`mv complex*pdb nonnative`;
	
	chdir "nonnative";
	for my $i (1..$Nstructs_nn) {
		print "$i ";
		`grep '^ATOM' complex.$i.pdb | head -n $last_a | grep ' CA ' > a.$i.pdb`;
		`grep '^ATOM' complex.$i.pdb | tail -n +$first_b | grep ' CA ' > b.$i.pdb`;
		`$exec_dir/PairDist.pl a.$i.pdb b.$i.pdb $exec_dir/aa_pairs.txt`; # > pairdist.$i.out`;
		`mv r12ij.txt r12ij.$i.txt`;
		`mv r10ij.txt r10ij.$i.txt`;
		`mv r6ij.txt r6ij.$i.txt`;
		`rm a.$i.pdb b.$i.pdb`;
	}
	`$exec_dir/CombineRij.sh $Nstructs_nn`;
	`mv r12ij.txt ../nn_r12ij.txt`;
	`mv r10ij.txt ../nn_r10ij.txt`;
	`mv r6ij.txt ../nn_r6ij.txt`;
	
	print "\n";
	chdir "..";
	`rm tmp*`;
	
}