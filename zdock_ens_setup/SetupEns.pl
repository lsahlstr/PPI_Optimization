#!/usr/bin/env perl

use strict;
use warnings;

my $Nstructs_n = 10;
my $Nstructs_nn = 100;
my $last_a = undef;
my $first_b = undef;
my $nf = undef;
my $zdock_dir = "/export/nVerde/users/lsahlstr/CGPPI_FF/decoys_bm4_zd3.0.2_6deg";
my $exec_dir = "/export/nVerde/users/lsahlstr/CGPPI_FF/opt5/bin";

# List of systems (PDB ID's) for which to construct native and non-native pools
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
	
	####################################################################
	## Native ensemble 
	####################################################################
	print "Native ensemble\n";
	`mkdir native`;
	# 1) Sort based upon RMSD values to get native-like poses
	`sort -g -k 9 tmp.out | head > tmp.native`;
	# 2) ZDOCK out file for natives
	`awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7}' tmp.native > tmp.native2`;
	`cat tmp.header tmp.native2 > native.out`;
	# 3) RMSD file for natives
	`awk 'BEGIN{print "0 0"}{print \$8,\$9}' tmp.native > n_rmsd.dat`;  # first line = "0 0" for native pose
	# 4) Generate complexes from ZDOCK out file
	`./my-create.pl native.out`;
	`mv complex*pdb lig.*.pdb native`;
	`cp -p $lig_pdb native/lig.0.pdb`;  # ligand in the native complex
	# 5) Minimize complexes and get pairwise distance information for the complexes
	for my $i (0..$Nstructs_n) {
		print "$i ";
		# Call external program that builds a Go model of the complex and performs the minimization 
		`$exec_dir/minCA.sh native $i`;		
		`mv r12ij.txt r12ij.$i.txt`;
		`mv r10ij.txt r10ij.$i.txt`;
		`mv r6ij.txt r6ij.$i.txt`;
	}
	# 6) Combine rij information from all native complexes
	`$exec_dir/CombineRij.sh $Nstructs_n`;
	`mv r12ij.txt n_r12ij.txt`;
	`mv r10ij.txt n_r10ij.txt`;
	`mv r6ij.txt n_r6ij.txt`;
	# 7) Check to see if all non-native structures are present
	chdir "native";
	$nf = `ls -1 min.*.log | wc -l`;
	if ($nf == ($Nstructs_n + 1)) {  # Need +1 because 10 + native = 11
	    `echo "n" > ../check`;
	}	
	# 8) Clean up files
	`tar czvf ligs.tar.gz lig.*.pdb`;
	`tar czvf min.tar.gz min.*.pdb`;
	`tar czvf log.tar.gz min.*.log`;
	`rm -f complex.*.pdb lig.*.pdb min.*.pdb min.*.log`; # -f because native pose is write-protected
	
	chdir "..";
	print "\n";
	
	####################################################################
	## Non-native ensemble 
	####################################################################
	print "Non-native ensemble\n";
	`mkdir nonnative`;
	# 1) Near-native (4 A < iRMSD < 8 A)
	`awk ' \$9 > 4 && \$9 < 8 ' tmp.out > tmp.nearn`;
	`sort -R tmp.nearn | head -n 20 > tmp.nearn2`;
	# 2) "Far from"-native (iRMSD >= 8 A)
	`awk ' \$9 >= 8 ' tmp.out > tmp.nn`;
	`sort -R tmp.nn | head -n 80 > tmp.nn2`;
	# 3) Combine near- and far from-native
	`cat tmp.nearn2 tmp.nn2 > tmp.nearn.nn`;
	# 4) ZDOCK out file for non-natives
	`awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7}' tmp.nearn.nn > tmp.nearn.nn2`;
 	`cat tmp.header tmp.nearn.nn2 > nonnative.out`;
 	# 5) RMSD file for non-natives
 	`awk '{print \$8,\$9}' tmp.nearn.nn > nn_rmsd.dat`;
 	# 6) Generate complexes from ZDOCK out file
 	`./my-create.pl nonnative.out`;
 	`mv complex*pdb lig.*.pdb nonnative`;
	# 7) Minimize complexes and get pairwise distance information for the complexes
	for my $i (1..$Nstructs_nn) {
		print "$i ";
		# Call external program that builds a Go model of the complex and performs the minimization
		`$exec_dir/minCA.sh nonnative $i`;
		`mv r12ij.txt r12ij.$i.txt`;
		`mv r10ij.txt r10ij.$i.txt`;
		`mv r6ij.txt r6ij.$i.txt`;
	}
	# 8) Combine rij information from all non-native complexes
	`$exec_dir/CombineRij.sh $Nstructs_nn`;
	`mv r12ij.txt nn_r12ij.txt`;
	`mv r10ij.txt nn_r10ij.txt`;
	`mv r6ij.txt nn_r6ij.txt`;
	# 9) Check to see if all non-native structures are present
	chdir "nonnative";
	$nf = `ls -1 min.*.log | wc -l`;
	if ($nf == $Nstructs_nn) {
	    `echo "nn" >> ../check`;
	}
	# 10) Clean up files
	`tar czvf ligs.tar.gz lig.*.pdb`;
	`tar czvf min.tar.gz min.*.pdb`;
	`tar czvf log.tar.gz min.*.log`;
	`rm complex.*.pdb lig.*.pdb min.*.pdb min.*.log`;
	chdir "../";
	
	`rm tmp*`;
	
	chdir "../";
	print "\n\n";
	
}