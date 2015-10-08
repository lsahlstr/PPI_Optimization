#!/usr/bin/env perl

use strict;
use warnings;

my $pdblist = $ARGV[0];
open(PDBS,"$pdblist");
my @pdbs = <PDBS>;

foreach my $pdb (@pdbs) {
    chomp $pdb;
    print "Submitting $pdb\n";
    #`sed "s/SYSTEM/$pdb/g" template_setup.pbs > $pdb.pbs`;
    `qsub $pdb.pbs`;
    sleep(2);
    print ".\n"; 
}
