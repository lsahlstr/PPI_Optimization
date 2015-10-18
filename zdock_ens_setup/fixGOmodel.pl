#!/usr/bin/env perl

sub usage {
  printf STDERR "\nUsage:   fixGOmodel.pl [-options] <aa_pdb>\n";
  printf STDERR "Options: [-tag name]\n";
  printf STDERR "         [-gotag name]\n";
  printf STDERR "         [-dir name]\n";
  printf STDERR "\n\n";
  exit 1;
}

use vars qw ( $perllibdir );

BEGIN {
    $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
      ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use warnings;
use strict;

use GenUtil;
use Molecule;

my $pdb=undef;
my $gopdb=undef;
my $goprm=undef;
my $gotop=undef;
my $gotag="GO_";
my $tag="new";
my $dir=".";

while ($#ARGV>=0){
  if ($ARGV[0] eq "-tag"){
    shift @ARGV;
    $tag=shift @ARGV;
  }
  elsif ($ARGV[0] eq "-gotag"){
    shift @ARGV;
    $gotag=shift @ARGV;
  }
  elsif ($ARGV[0] eq "-dir"){
    shift @ARGV;
    $dir=shift @ARGV;
  }
  elsif ($ARGV[0] eq "-h" || $ARGV[0] eq "-help"){
    &usage();
  }
  else{
    $pdb=shift @ARGV;
  }
}

die "\nError: Please provide an all-atom PDB file\n\n" if (!defined $pdb || ! -e $pdb);

my $mol=Molecule::new();
$mol->readPDB($pdb);

my ($ref_start,$ref_stop,$ref_prmstart, $ref_prmstop,$ref_inx,$ref_rinx)=&getRanges($mol);

if (&countChains($mol)){
  opendir (DIR, $dir) || die "\nError: Directory \"$dir\" cannot be opened for reading\n\n";

  while (my $file = readdir(DIR)){
    next if ($file !~ /^$gotag/ || $file =~ /^GO_prot/);
    my $out="$tag.$file";
    printf "$file $out\n";
    if ($file =~ /param$/){
      &genNewPrm($file,$out,$ref_prmstart,$ref_prmstop,$ref_inx,$mol);
    }
    elsif ($file =~ /top$/){
      &genNewTop($file,$out,$ref_inx,$mol);
    }
    elsif ($file =~ /pdb$/){
      my $gomol=Molecule::new();
      $gomol->readPDB($file,ignoreseg=>1);
      &genNewPDB($gomol,$out,$ref_start,$ref_stop,$ref_inx,$ref_rinx,$mol)
    }
    else{
      next;
    }
  }

  closedir (DIR);
}

###################
#                 #
#   SUBROUTINES   #
#                 #
###################

sub countChains {
  my $self=shift;
  
  my $nchains=0;

  foreach my $c (@{$self->activeChains()}){
    #printf "$c->{id}\n";
    $nchains++;
  }
  return $nchains;
}

sub getRanges {
  my $self=shift;
  my @start=();
  my @stop=();
	my @prmstart=();
	my @prmstop=();
  my %inx=();
  my %rinx=(); #Reverse lookup

  my $chainid="A";
  my $segid="PROA";

	my %chainlu=(); #Lookup
  #Can conversion for each chain to original PDB ids
  foreach my $c (@{$self->activeChains()}){
    $chainlu{$chainid}=$c->{id};
    $chainid++;
    $segid++;
  }

  $chainid="A";
  $segid="PROA";

  my $natom=0;
  foreach my $c (@{$self->activeChains()}){
    my $inx=$#start+1;
    foreach my $r (@{$c->{res}}){
      $natom++;
#      $start[$inx]=$r->{num} if (!defined $start[$inx]);
#      $stop[$inx]=$r->{num}; #No restriction
			$prmstart[$inx]=$r->{num} if (!defined $prmstart[$inx]);
			$prmstop[$inx]=$r->{num};
			$start[$inx]=$natom if (!defined $start[$inx]);
			$stop[$inx]=$natom;
			$inx{"G$natom"}="$chainlu{$chainid}$r->{num}";
      $rinx{"$chainid$r->{num}"}=$natom;
    }
    $chainid++;
    $segid++;
  }

  for (my $i=0; $i<=$#start; $i++){
    #print "$start[$i] $stop[$i]\n";
  }

  return (\@start,\@stop,\@prmstart,\@prmstop,\%inx,\%rinx);
}

sub genNewPDB {
  my $self=shift;
  my $out=shift;
  my $ref_start=shift;
  my $ref_stop=shift;
  my $ref_inx=shift;
  my $ref_rinx=shift;
	my $ref_mol=shift;

  my @start=@{$ref_start};
  my @stop=@{$ref_stop};
  my %inx=%{$ref_inx};
  my %rinx=%{$ref_rinx};

  my $chainid="A";
  my $segid="PROA";

  my $finalmol=Molecule::new();

	my %chainlu=(); #Lookup
	#Conversion for each chain to original PDB ids
	foreach my $c (@{$ref_mol->activeChains()}){
		$chainlu{$chainid}=$c->{id};
		$chainid++;
		$segid++;
	}

	$chainid="A";
	$segid="PROA";

  for (my $i=0; $i<=$#start; $i++){
    #print "$start[$i] $stop[$i]\n";
    my $tstart=$rinx{"$chainid$start[$i]"};
    my $tstop=$rinx{"$chainid$stop[$i]"};

    #print "$tstart $tstop\n";
		my $natom=0;
		my $nres=0;
    #Only first chain is renumbered starting with residue 1
    #$self->setValidSelection("$tstart-$tstop");
		$self->resetValidResidues(0);
		foreach my $c ( @{$self->activeChains()} ) {
    	foreach my $r ( @{$c->{res}} ) {
				$nres++;
				if ($nres >= $start[$i] && $nres <= $stop[$i]){
      		$r->{valid}=1;
				}
			}
     	foreach my $a ( @{$c->{atom}} ) {
				$natom++;
				if ($natom >= $start[$i] && $natom <= $stop[$i]){
  				$a->{valid}=1;
				}
     	}
   	}

    my $tmpmol=Molecule::new();
    $tmpmol=$self->clone(1);

    $tmpmol->setChain($chainid);
    $tmpmol->generateSegNames($segid);

    foreach my $c (@{$tmpmol->activeChains()}){
			$c->{id}=$chainlu{$chainid};	
      foreach my $a (@{$c->{atom}}){	
        die "\nError: Total number of residues (999) exceeded\n\n" if ($a->{resnum} > 999 || $a->{atominx} > 999);
#        $a->{resname}="G$a->{resnum}";
				$a->{chain}=$chainlu{$chainid};
        if ($chainid eq "A"){
          $a->{resname}=$inx{"G$a->{resnum}"};
					$a->{resname}=~ s/A/$chainlu{$chainid}/;
        }
        else{
#          $a->{resname}="$chainid$a->{resnum}";
					$a->{resname}="$chainlu{$chainid}$a->{resnum}";
        }
      }
    }

    $finalmol->merge($tmpmol);
    $chainid++;
    $segid++;
  }

  $finalmol->writePDB($out);

  return;
}

sub genNewPrm {
  my $inp=shift;
  my $out=shift;
  my $ref_start=shift;
  my $ref_stop=shift;
  my $ref_inx=shift;
	my $ref_mol=shift;

  my @start=@{$ref_start};
  my @stop=@{$ref_stop};
  my %inx=%{$ref_inx};

  open (INP, "$inp") || die "\nError: Cannot open file \"$inp\"\n\n";
  open (OUT, ">$out");

  my @s=();
  my $i;

  my $param=undef;
  my $flag=1;
  while (<INP>){
    if ($flag){
      $flag=0 if ($_ =~ /read param card/);
      next;
    }
    elsif ($_ =~ /(BOND|ANGLE|DIHEDRAL|NONBOND|NBFIX)/){
      $param=$1;
      printf OUT $_;
    }
    else{
      if ($_ =~ /^\n/){
        printf OUT $_;
      }
      elsif (defined $param){
        @s=split/\s+/,$_;
				my %chains=();
        if (defined $s[0] && $s[0] =~ /^G/){
          $s[0]=~s/$s[0]/$inx{$s[0]}/g;
					$chains{substr($s[0],0,1)}=1;
          $s[0]=~s/([A-Z])(\d+)/$2/g;
        }
        if (defined $s[1] && $s[1] =~ /^G/){
          $s[1]=~s/$s[1]/$inx{$s[1]}/g;
					$chains{substr($s[1],0,1)}=1;
          $s[1]=~s/([A-Z])(\d+)/$2/g;
        }
        if (defined $s[2] && $s[2] =~ /^G/){
          $s[2]=~s/$s[2]/$inx{$s[2]}/g;
					$chains{substr($s[2],0,1)}=1;
          $s[2]=~s/([A-Z])(\d+)/$2/g;
        }
        if (defined $s[3] && $s[3] =~ /^G/){
          $s[3]=~s/$s[3]/$inx{$s[3]}/g;
					$chains{substr($s[3],0,1)}=1;
          $s[3]=~s/([A-Z])(\d+)/$2/g;
        }
        $_ =~ s/^(G\d+)(\s+)/$inx{$1}$2/g;
        $_ =~ s/(\s+)(G\d+)/$1$inx{$2}/g;
        if ($param =~ /^BOND/){
          for ($i=0; $i<=$#start; $i++){
            if (scalar keys %chains == 1 && $s[0] >= $start[$i] && $s[0] <= $stop[$i] 
                && $s[1] >= $start[$i] && $s[1] <= $stop[$i]){
              printf OUT $_;
              last;
            }
          }
        }
        elsif ($param =~ /^ANGLE/){
          for ($i=0; $i<=$#start; $i++){
            if (scalar keys %chains == 1 && $s[0] >= $start[$i] && $s[0] <= $stop[$i]
                && $s[1] >= $start[$i] && $s[1] <= $stop[$i]
                && $s[2] >= $start[$i] && $s[2] <= $stop[$i]){
              printf OUT $_;
              last;
            }
          }
        }
        elsif ($param =~ /^DIHEDRAL/){
          for ($i=0; $i<=$#start; $i++){
            if (scalar keys %chains == 1 && $s[0] >= $start[$i] && $s[0] <= $stop[$i] 
                && $s[1] >= $start[$i] && $s[1] <= $stop[$i]
                && $s[2] >= $start[$i] && $s[2] <= $stop[$i]
                && $s[3] >= $start[$i] && $s[3] <= $stop[$i]){
              printf OUT $_;
              last;
            }
          }
        }
				elsif ($param =~ /^NBFIX/){
					if (defined $s[3] && $s[3] > 13){
						printf STDERR "\nWarning: NBFIX distance is larger than 13 Angstroms!\n";
						printf STDERR "Warning: Please verify that this contact is correct!\n";
						printf STDERR "$_\n";
					}
					printf OUT $_;
				}
        else{
          printf OUT $_;
        }
      }
      else{
        printf OUT $_;
      }
    }
  }

  close(OUT);
  close(INP);

  return;
}

sub genNewTop {
  my $inp=shift;
  my $out=shift;
  my $ref_inx=shift;

  my %inx=%{$ref_inx};

  open (INP,$inp) || die "\nError: Could not open file \"$inp\"\n\n";
  open (OUT,">$out");

  my $flag=1;
  while (<INP>){
    if ($flag){
      $flag=0 if ($_ =~ /read rtf card/);
      next;
    }
    else{
      $_ =~ s/(G\d+)/$inx{$1}/g;
      printf OUT $_;
    }
  }

  close(OUT);
  close(INP);

  return;
}
