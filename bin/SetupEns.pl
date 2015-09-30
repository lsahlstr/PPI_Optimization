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
	
	`cp -p $exec_dir/my-create.pl .`;
	
	# Combine Zdock output and RMSD files
	`head -n 5 $zdock_out > tmp.header`;
	`cat tmp.header $zdock_rmsds > tmp.rmsds`;
	`paste $zdock_out tmp.rmsds | tail -n +6 > tmp.out`;
	
	#### Native-like ensemble ####
	print "Native ensemble\n";
	`mkdir native`;
	# Sort based upon RMSD values to get native-like poses
	`sort -g -k 9 tmp.out | head > tmp.native`;
	# ZDOCK out file for natives
	`awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7}' tmp.native > tmp.native2`;
	`cat tmp.header tmp.native2 > native.out`;
	# RMSD file for natives
	`awk '{print \$8,\$9}' tmp.native > n_rmsd.dat`;
	# Generate complexes from ZDOCK out file
	`./my-create.pl native.out`;
	`mv complex*pdb lig.*.pdb native`;
	# Get pairwise distance information for the complexes
	chdir "native";
	for my $i (1..$Nstructs_n) {
		print "$i ";
		`$exec_dir/PairDist.pl ../rec.pdb lig.$i.pdb $exec_dir/aa_pairs.txt`; # > pairdist.$i.out`;
		`mv r12ij.txt r12ij.$i.txt`;
		`mv r10ij.txt r10ij.$i.txt`;
		`mv r6ij.txt r6ij.$i.txt`;
	}
	# Combine rij information from all non-native complexes
	`$exec_dir/CombineRij.sh $Nstructs_n`;
	`mv r12ij.txt ../n_r12ij.txt`;
	`mv r10ij.txt ../n_r10ij.txt`;
	`mv r6ij.txt ../n_r6ij.txt`;
	
	`tar czvf native_ligs.tar.gz lig.*.pdb`;
	`rm complex.*.pdb lig.*.pdb`;
	
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
	# ZDOCK out file for non-natives
	`awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7}' tmp.nearn.nn > tmp.nearn.nn2`;
 	`cat tmp.header tmp.nearn.nn2 > nonnative.out`;
 	# RMSD file for non-natives
 	`awk '{print \$8,\$9}' tmp.nearn.nn > nn_rmsd.dat`;
 	# Generate complexes from ZDOCK out file
 	`./my-create.pl nonnative.out`;
 	`mv complex*pdb lig.*.pdb nonnative`;
	# Get pairwise distance information for the complexes
	chdir "nonnative";
	for my $i (1..$Nstructs_nn) {
		print "$i ";
		`$exec_dir/PairDist.pl ../rec.pdb lig.$i.pdb $exec_dir/aa_pairs.txt`; # > pairdist.$i.out`;
		`mv r12ij.txt r12ij.$i.txt`;
		`mv r10ij.txt r10ij.$i.txt`;
		`mv r6ij.txt r6ij.$i.txt`;
		#`rm lig.$i.pdb`;
	}
	# Combine rij information from all non-native complexes
	`$exec_dir/CombineRij.sh $Nstructs_nn`;
	`mv r12ij.txt ../nn_r12ij.txt`;
	`mv r10ij.txt ../nn_r10ij.txt`;
	`mv r6ij.txt ../nn_r6ij.txt`;
	
	`tar czvf nonnative_ligs.tar.gz lig.*.pdb`;
	`rm complex.*.pdb lig.*.pdb`;
	
	chdir "../";
	`rm tmp*`;
	
	chdir "../";
	print "\n\n";
	
}