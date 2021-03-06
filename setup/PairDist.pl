#!/usr/bin/env perl

use strict;
use warnings;

my $a_mol = $ARGV[0];       # PDB file for mol A
my $b_mol = $ARGV[1];       # PDB file for mol B
my $aa_pairs = $ARGV[2]; 	# File containing all 210 possible amino acid pairs
my $cutoff = 20;            # Cutoff within which to store i,j interactions
my %a_chains = ();          # hash for molecule A info
my %b_chains = ();          # hash for molecule B info
my %r = ();             	# hash for storing rij values for given inter-protein residue pair
#`mkdir -p $outdir`;           


#############################################################################
# Molecule info for molecule A (C-alpha coordinates)
#############################################################################
open (PDBA,"$a_mol");
while ($_ = <PDBA>) {

    if ($_ =~ /^ATOM.*CA/) {
        my $altconf=substr($_,16,1);
        if (($altconf eq " ") || ($altconf eq "A")) { 
            my $chain=substr($_,21,1);
	    	(my $resname=substr($_,17,4))=~s/ //g;
            my $resnum=substr($_,22,6);    
            my $xcoor=substr($_,30,8)+0.0;
	    	my $ycoor=substr($_,38,8)+0.0;
            my $zcoor=substr($_,46,8)+0.0;
            my $t;                         # initialize temporary "hash"

            if (exists $a_chains{$chain}) {  # does hash exist for current chain?
                $t={};                     # initialize $t as empty hash
                $t->{resname}=$resname;    # store resname to hash t
                $t->{resnum}=$resnum;      # store resnum
                $t->{xcoor}=$xcoor;        # store x coordinate
                $t->{ycoor}=$ycoor;        # store y coorindate
                $t->{zcoor}=$zcoor;        # store z coordinate
                push(@{$a_chains{$chain}}, $t);  # store contents of $t to %chains
 
	    	} else {
                $a_chains{$chain} = [];      # chain array with molecule hash
	        	$t={};
	        	$t->{resname}=$resname;
				$t->{resnum}=$resnum;
				$t->{xcoor}=$xcoor;
				$t->{ycoor}=$ycoor;
				$t->{zcoor}=$zcoor;
				
                push(@{$a_chains{$chain}}, $t);
            }
        }
    }
}

my @c_a = ();  # array with a_mol chains
foreach my $key (sort keys %a_chains) {   
  push(@c_a, $key);                  
}
#my $n_spy = $#c_a + 1; # number of spy conformers

#############################################################################
# Molecule info for molecule B (C-alpha coordinates)
#############################################################################
open (PDBB,"$b_mol");
while ($_ = <PDBB>) {

    if ($_ =~ /^ATOM.*CA/) {
        my $altconf=substr($_,16,1);
        if (($altconf eq " ") || ($altconf eq "A")) { 
            my $chain=substr($_,21,1);
	    	(my $resname=substr($_,17,4))=~s/ //g;
            my $resnum=substr($_,22,6);    
            my $xcoor=substr($_,30,8)+0.0;
	    	my $ycoor=substr($_,38,8)+0.0;
            my $zcoor=substr($_,46,8)+0.0;
            my $t;                         # initialize temporary "hash"

            if (exists $b_chains{$chain}) {  # does hash exist for current chain?
                $t={};                     # initialize $t as empty hash
                $t->{resname}=$resname;    # store resname to hash t
                $t->{resnum}=$resnum;      # store resnum
                $t->{xcoor}=$xcoor;        # store x coordinate
                $t->{ycoor}=$ycoor;        # store y coorindate
                $t->{zcoor}=$zcoor;        # store z coordinate
                push(@{$b_chains{$chain}}, $t);  # store contents of $t to %chains
 
	    	} else {
                $b_chains{$chain} = [];      # chain array with molecule hash
	        	$t={};
	        	$t->{resname}=$resname;
				$t->{resnum}=$resnum;
				$t->{xcoor}=$xcoor;
				$t->{ycoor}=$ycoor;
				$t->{zcoor}=$zcoor;
				
                push(@{$b_chains{$chain}}, $t);
            }
        }
    }
}

my @c_b = ();  # array with b_mol chains
foreach my $key (sort keys %b_chains) {   
  push(@c_b, $key);                  
}
#my $n_spy = $#c_a + 1; # number of spy conformers


#############################################################################
# Get inter-molecular Ca-Ca distances
#############################################################################
open(OUT,">rij.txt");
for (my $i=0; $i<=$#c_a; $i++) {  # for each mol_a chain
    for (my $j=0; $j<=$#c_b; $j++){  # for each mol_b chain

		foreach my $ai (@{$a_chains{$c_a[$i]}}) {
			foreach my $aj (@{$b_chains{$c_b[$j]}}) {
		
				# atom i
				my $resname_i = $ai->{resname};
				my $resnum_i = $ai->{resnum};
				my $xcoor_i = $ai->{xcoor};
				my $ycoor_i = $ai->{ycoor};
				my $zcoor_i = $ai->{zcoor};
	
				# atom j
				my $resname_j = $aj->{resname};
				my $resnum_j = $aj->{resnum};
				my $xcoor_j = $aj->{xcoor};
				my $ycoor_j = $aj->{ycoor};
				my $zcoor_j = $aj->{zcoor};
			
				# Identity of residue pair
				my $res_pair;
				if ($resname_i lt $resname_j) {
					$res_pair = $resname_i.$resname_j;
				} else {
					$res_pair = $resname_j.$resname_i;
				}
			
				# Distance information for intermolecular potential
				my $r_ij = abs(sqrt((($xcoor_j-$xcoor_i)**2) 
					+ (($ycoor_j-$ycoor_i)**2)
					+ (($zcoor_j-$zcoor_i)**2)));
				
				if ($r_ij <= $cutoff) {
					printf OUT "%s\t%8.3f\n", $res_pair,$r_ij;
				}
			}
		}
    }
}

# Fill in the gaps in %r12, %r10, %r6:
# Not all 210 possible inter-chain amino acid pairs may be represented
# in the A-B interactions, so just put a zero for them
# open (RMIN,"$aa_pairs");
# while (my $l = <RMIN>) {
#     my @l = split(' ',$l);
#     my $key = $l[0];
#     
#     if (!exists $r{$key}) {
#         $r{$key} = 0;
#     }
# }
# close RMIN;
# 
# # Sort hash keys in alphabetical order
# open(OUT,">rij.txt");
# foreach my $key (sort keys %r) {
#     my $val = $r{$key};
#     print OUT "$key $val\n";
# }
# close OUT;