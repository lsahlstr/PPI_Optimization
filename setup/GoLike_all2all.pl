#!/usr/bin/env perl

####################################################################
#  Script to append all possible pairwise interactions between two 
#  protein chains to the end of a Go model parameter file. 
#
#  c. Logan S. Ahlstrom 2014
####################################################################

use strict;
use warnings;

sub usage {

  printf STDERR "\nUsage:   GoLike_all2all.pl [-options] <aa_pdb>\n";
  printf STDERR "Options: [-sig_file]\n";
  printf STDERR "         [-ene_file]\n";
  printf STDERR "         [-x_inter]";
  printf STDERR "         [-x_intra]";
  printf STDERR "         [-param_file]\n";
  printf STDERR "         [-tag\n]";
  printf STDERR "         [-pcalign]";
  printf STDERR "\n\n";
  exit 1;

}

my $sig_file;             # file with Ca-Ca distances for residue pairs
my $ene_file;             # file with MJ contact energies for residue pairs 
my $pdb_file;             # all-atom PDB (same one that was uploaded to the Go server)
my $param_file;           # original parameter file created by Go server
my $x_inter = 1.0;        # scaling factor for native INTERmolecular contacts
my $x_intra = 1.0;        # scaling factor for native INTRAmolecular contacts
my $x_all2all = 1.0;      # scaling factor for all possible pairwise interactions
my $tag;                  # tag for name of new parameter file (prepended to name of old one)
my $pcalign;              # flag to specify if nativer INTERmolecular interacions are from PCAlign
my $MJ_avg = -3.1757;     # average of all 210 MJ contact energies
my $MJ_scale = $MJ_avg*2; # scaling factor for MJ energies      
my $N = 0;                # number of orientational contacts to include in non-native contacts
my $M = 0;                # number of hydrogen bonds to include in non-native contacts
my $ii_scale = 10000;     # scaling factor for volume exclusion energies 
my %sig = ();             # hash for residue vdW radii
my %ene = ();             # hash for MJ energies
my %chains = ();          # hash for chain, resname, resnum info for each residue
my %natives = ();         # has for native contacts
my @ssbond = ();          # array for chain and resnum for disulfide bond cysteins

while ($#ARGV>=0){
    if ($ARGV[0] eq "-sig_file"){
        shift @ARGV;
        $sig_file = shift @ARGV;
  
    } elsif ($ARGV[0] eq "-ene_file"){
        shift @ARGV;
        $ene_file = shift @ARGV;

    } elsif ($ARGV[0] eq "-param_file"){
        shift @ARGV;
        $param_file = shift @ARGV;

    } elsif ($ARGV[0] eq "-x_inter"){
        shift @ARGV;
        $x_inter = shift @ARGV;

    } elsif ($ARGV[0] eq "-x_intra"){
        shift @ARGV;
        $x_intra = shift @ARGV;
    
    } elsif ($ARGV[0] eq "-x_all2all"){
        shift @ARGV;
        $x_all2all = shift @ARGV;

    } elsif ($ARGV[0] eq "-tag"){
        shift @ARGV;
        $tag = shift @ARGV;

    } elsif ($ARGV[0] eq "-pcalign"){
        shift @ARGV;
        $pcalign = 1;
 
    } elsif ($ARGV[0] eq "-h" || $ARGV[0] eq "-help"){
        &usage();
  
    } else{
        $pdb_file=shift @ARGV;
    }
}

die "\nError: Please provide an all-atom PDB file\n\n" if (!defined $pdb_file || ! -e $pdb_file);


# Create hash for referencing residue van der Waals radii
open (SIG,"$sig_file");
while (my $l = <SIG>) {
    chomp($l);
    my($key,$val) = split("\t", $l);
    $sig{$key} = $val; 
}

# Create hash for referencing MJ contact energies
open (ENE,"$ene_file");
while (my $l = <ENE>) {
    chomp($l);
    my($key,$val) = split("\t", $l);
    $ene{$key} = $val;
}


#############################################################################
##  Construct reference hash (%chains{}) with all chain and residue name/number info
##  Format of this section is MMTSB style
#############################################################################
open (PDB,"$pdb_file");  # all-atom PDB originally processed by Go server 
while (my $_ = <PDB>) {  # go through all-atom PDB file

    if ($_ =~ /^ATOM.*CA/) {  # Find C-alpha lines 

        my $chain=substr($_,21,1);
        (my $resname=substr($_,17,4))=~s/ //g;
        my $resnum=substr($_,22,6);
        my $xcoor=substr($_,30,8)+0.0;
        my $ycoor=substr($_,38,8)+0.0;
        my $zcoor=substr($_,46,8)+0.0; 
        my $t;             # variable for temporary hash below
        
        if (exists $chains{$chain}) {       # check to see if a hash exists for the current chain
            
            $t={};                          # initialize $t as empty hash
            $t->{resname}=$resname;         # store resname to hash t
            $t->{resnum}=$resnum;           # store resnum
            $t->{xcoor}=$xcoor;             # store x coordinate
            $t->{ycoor}=$ycoor;             # store y coorindate
            $t->{zcoor}=$zcoor;             # store z coordinate
            
            push(@{$chains{$chain}}, $t);   # store contents of $t to %chains
            
        } else {
            
            $chains{$chain} = [];           # Create new hash; now a ref to a new empty array
            $t={};                          
            $t->{resname}=$resname;
            $t->{resnum}=$resnum;
            $t->{xcoor}=$xcoor;             
            $t->{ycoor}=$ycoor;             
            $t->{zcoor}=$zcoor;             

	    push(@{$chains{$chain}}, $t);   
            
        }
    }

    # For systems with disulfide bonds:
    # This info is not currently taken into account, but could be used to 
    # construct a stream file that specifies a harmonic restraint for each S-S bond
    if (/^SSBOND/) {
       
        (my $chain1=substr($_,15,1))=~s/ //g;
        (my $resnum1=substr($_,16,5))=~s/ //g;
        (my $chain2=substr($_,29,1))=~s/ //g;
        (my $resnum2=substr($_,30,5))=~s/ //g;

        print "$chain1\t$resnum1\t$chain2\t$resnum2\n";
        
        my $t={}; 
        $t->{chain1}=$chain1;
        $t->{resnum1}=$resnum1;
        $t->{chain2}=$chain2;;
        $t->{resnum2}=$resnum2;
        
        push(@ssbond, $t);
        
    }
}

foreach my $a (@ssbond) {   # foreach chain in the hash %chain
    print "$a->{resnum1}\n";
}


#############################################################################
##  Read through original parameter file and 1) modify volume exclusion energies
##  and 2) make lists of intra- and intermolecular contacts. Print the lines 
##  to the new parameter file. 
############################################################################# 
my $new_param_file = $tag.".".$param_file;  # name for new parameter file
open (OLD,"$param_file");                   # read original parameter file
open (NEW,">$new_param_file");              # open file for writing new parameter file
open (INTRA,">intra_native.cont");          # list of native INTRAmolecular contacts
if ($x_inter > 0) {
    open (INTER,">inter_native.cont");      # list of native INTERmolecular contacts
}
my $count_lines = 1;  # counter for number of lines in the original parameter file
my $num_lines  = `wc -l $param_file | awk '{print \$1}'` - 3; # exclude last three lines from total count

# Go through old parameter file
while (my $l = <OLD>) { 

    if ($count_lines <= $num_lines) {  # If not at the end of the file ("END" statement)

        chomp($l);
        
        # Matches NONBOND section (modify volume exclusion energies)
        if ($l =~ /^(\w)(\d+)\s+0.0\s+-0.000132\s+(\d+.\d+)/) {  # lines w/ vdW nonbond info
            my $c = $1;  # chain id
            my $n = $2;  # residue number
            my $r = $3;  # repulsive radius
      
            my @chn=@{$chains{$c}};   # Grab hash corresponding to the current chain
            for (my $atom=0; $atom<=$#chn; $atom++) { # For each Go atom listed in NONBOND section 
                my $a=$chn[$atom];    # Part of the hash that contains the atom info
    
                if ($a->{resnum} == $n) {  # Compare resnum on current line with reference resnum
                    my $s = $a->{resname}.$a->{resname}; # Lookup string for MJ energy file
                    my $ene_ii = $ene{$s};               # MJ eps_ii self term
                    
                    # Scale the volume exclusion energies so that they are similar in 
                    # magnitude to the original value (-0.000132)
                    my $ene_ii_scaled = $ene_ii/$ii_scale; 
                    
                    printf NEW "$c"."$n"."     0.0  %9.6f  %8.6f\n",$ene_ii_scaled,$r;
                }
            }
        
            $count_lines++;    # increment line counter
        
        # Matches NBFIX section of Go model parameter file
        } elsif ($l =~ /^([A-Z])(\d+)\s+([A-Z])(\d+)\s+(-\d\.\d+)\s+(\d+\.\d+)/) {
            my $ch1  = $1;
            my $num1 = $2;
            my $ch2  = $3;
            my $num2 = $4;
            my $ene  = $5;
            my $dist = $6; 
            my $newene;

            my $string = $ch1.$num1." ".$ch2.$num2." ".$dist; # line for contact lists

            # Print Go-type interactions to new param file and to contact lists
            if ($ch1 eq $ch2) { # INTRAmolecular   
                $newene = $ene * $x_intra;  # scale epsilon 
                $l =~ s/$ene/$newene/g;  # replace old epsilon value with new one
                print NEW "$l\n";   # new param file
                print INTRA "$string\n"; # contact list
            
            } elsif ($ch1 ne $ch2) {  # INTERmolecular
                if ($x_inter > 0) { 
                    my $value = $ch1.$num1.$ch2.$num2;  # Fill in hash of native INTER interactions
                    $natives{$value} = 1;
                    if (!defined $pcalign) {
                        $newene = $ene * $x_inter;  # scale epsilon
                        $l =~ s/$ene/$newene/g;  # replace old epsilon value with new one
                        print NEW "$l\n";   # new param file 
                        print INTER "$string\n";
                    }
                }
            }
            
            $count_lines++;       
            
        } else {
            print NEW "$l\n";  # print all other lines
            $count_lines++;
        }     
    }
}


#############################################################################
##  Append all possible pairwise INTERmolecular interactions to the NBFIX 
##  section of the new parameter file
#############################################################################
my @c=();                          # array for all chains, to loop through below
foreach my $key (sort keys %chains) {   # foreach chain in %chains
  push(@c, $key);                  # store info for each chain to the new array
}

my $NNcount=0; # counter for "non-native" contacts
my $Ncount=0;  # counter for native contacts
# pairwise comparison of residues between chains ("inter")
if ($x_all2all > 0) {
    for (my $i=0; $i<=$#c; $i++){  # foreach "source" chain
        for (my $j=$i+1; $j<=$#c; $j++){  # foreach "target" chain
            foreach my $ai (@{$chains{$c[$i]}}) {  # foreach residue in "source" chain 
                foreach my $aj (@{$chains{$c[$j]}}){  # foreach residue in "target" chain
                    
                    my $num_i = $ai->{resnum};   # resnum from source chain 
                    my $num_j = $aj->{resnum};   # resnum from target chain
                    my $res_i = $ai->{resname};  # resname from source chain
                    my $res_j = $aj->{resname};  # resname from target chain
                    
                    my $res_pair;  # residue pair for ene_ij and sig_ij lookup
                    if ($res_i lt $res_j) {
                        $res_pair = $res_i.$res_j;
                    } else {
                        $res_pair = $res_j.$res_i;
                    }
                    
                    (my $s1="$c[$i]$num_i")=~s/ //g;   # string of chain+resnum (e.g., A85)
                    (my $s2="$c[$j]$num_j")=~s/ //g;   # string of chain+resnum (e.g., B87)
                   
                    # New r_ij (improperly named sigma... will fix this...)
                    my $sig_ij = $sig{$res_pair}; # Ca-Ca distance for the source-target pair

                    # New epsilon_ij  
                    my $eps_ij = $ene{$res_pair}; # MJ contact energy for the source-target pair
                    if ($x_inter > 0) {  # Only include interactions that are not already present in native INTER set
                        my $s12 = $s1.$s2; # Strings to check if residue pair exists in native contact hash
                        my $s21 = $s2.$s1;
                        if ((!exists $natives{$s12}) && (!exists $natives{$s21})) { # Only add if non-native
                            my $eps_ij_scaled = -1.0 * $x_all2all * $eps_ij/$MJ_scale + 0.25*$N + $M; # Scaled MJ energy
                            printf NEW $s1."    ".$s2."       %9.6f   %9.6f\n",$eps_ij_scaled,$sig_ij;
                            $NNcount++;  # increment counter for non-native interactions
                        } elsif ((exists $natives{$s12}) || (exists $natives{$s21})) {
                            my $eps_ij_scaled = -1.0 * $x_inter * $eps_ij/$MJ_scale + 0.25*$N + $M; # Scaled MJ energy
                            printf NEW $s1."    ".$s2."       %9.6f   %9.6f\n",$eps_ij_scaled,$sig_ij;
                            $Ncount++;  # increment counter for non-native interactions 
                        }

                    } elsif ($x_inter == 0) { # Include all pairwise interactions
                        my $eps_ij_scaled = -1.0 * $x_all2all * $eps_ij/$MJ_scale + 0.25*$N + $M; # Scaled MJ energy
                        printf NEW $s1."    ".$s2."       %9.6f   %9.6f\n",$eps_ij_scaled,$sig_ij;
                    }
                }
            }
        }
    }
}

print NEW "\nEND\n\n";
## will only be non-zero if $x_inter > 0
print "$Ncount\n";
print "$NNcount\n";
