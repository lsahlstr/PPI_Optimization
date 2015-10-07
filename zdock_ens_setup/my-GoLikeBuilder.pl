#!/usr/bin/perl

# Build a Go model from a PDB structure

# Requires the MMTSB toolset (including charmm)

# Requires as command-line arguments:
#			1) tag ($ARGV[0])
#			2) pdb filename ($ARGV[1])

# J. Karanicolas, 2003

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
  #open( outfile,">>go_error.log" );
  #system( "env>>go_error.log");
  #print outfile "\n mmtsbdir ",$ENV{MMTSBDIR}," charmmexec ",$ENV{CHARMMEXEC}," charmmdata\n",$ENV{CHARMMDATA};
}

use lib $perllibdir;
use lib ".";

use strict;

use IO::File;
use IO::Handle;

use GenUtil;
use Molecule;
use CHARMM;

my $ARGC = @ARGV;
if (($ARGC < 2) || ($ARGV[0] eq "-h") || ($ARGV[0] eq "-help")) {
  printf "\nProgram requires a tag and a pdb filename.\n";
  die "Please try again with appropriate command-line arguments.\n\n";
}

my $prefix=$ARGV[0];
my $pdbfile=$ARGV[1];

&echoPrefix($prefix);

my %par = ( 

           # For defining Hbond contacts
           Hmaxdist_sq    =>   27.2484,              # 5.22 ^ 2
	   HBcut          =>   -0.5,
           KSfac          =>   27.888,               # 0.42 * 0.2 * 332
           H_factor       =>   1.0,

           # For defining sidechain contacts
           SC_dist_sq     =>   20.25,                # 4.5 ^ 2
           SC_CA_dist_sq  =>   196,                  # 14 ^ 2
           MJ_ave         =>   -3.1757,
           MJ_factor      =>   0.5,
           MJfile         =>   "/export/nVerde/users/lsahlstr/go_interact/bin/MJ_values.list",
           #MJfile         =>   "/export/www/mmtsb/data/MJ_values.list",

           # For charmm parameter file
           TfScaling      =>   0.0054,
           BondFactor     =>   200,
           AngleFactor    =>   40,
           DiheFactor     =>   0.4,
           RepFactor      =>   0.00007,
           RepSigFactor   =>   0.561231,
           DiheFile       =>   "/export/nVerde/users/lsahlstr/go_interact/bin/pseudodihedral.param"
           #DiheFile       =>   "/export/www/mmtsb/data/pseudodihedral.param"

	   );

# Setting Tf here scales everything such that the folding temperature will
# be approx. Tf. Setting TfScaling instead uses a fixed scaling factor.
my $Tf = 350.0;
my $TfScaling = undef;
$Tf = undef if ((defined $TfScaling) && (defined $Tf));

die "Unable to read input file" if (! -e $pdbfile);
my $inpdb=&GenUtil::getInputFile($pdbfile);
my $models=0;
while (<$inpdb>) {
  $models++ if (/^MODEL +([0-9]+)/);
}

my $mol;
my $qlist;


my $charmm=&CHARMM::new("cha.log");

$charmm->loadParameters(param=>19,dcel=>0.5,hsd=>"",hse=>"");
if ($models == 0) {
  $mol=&loadMolecule($charmm,$prefix,$pdbfile);
  $qlist=&generateQlist($mol);
} else {
  my $qarr=();
  $mol=&loadMolecule($charmm,$prefix,$pdbfile,1);
  push(@{$qarr},&generateQlist($mol));
  for (my $m=2; $m <= $models; $m++) {
    push(@{$qarr},
     &generateQlist(&loadMolecule($charmm,$prefix,$pdbfile,$m)));
  }
  $qlist=&aveQlist($qarr);
}
$charmm->finish();

# Write files for charmm
&writeGoPDB($prefix,$mol);
&writeGoSeq($prefix,$mol);
&writeGoTop($prefix,$mol);
&writeGoParam($prefix,$mol,$qlist,$Tf,$TfScaling);

# Write the list of native contacts
&writeQdetails($mol,$qlist,"GO_".$prefix.".Qdetails");
&writeQlist($mol,$qlist,"GO_".$prefix.".Qlist");

#printf "All done\n";
exit(0);



sub echoPrefix {
  my $prefix=shift;

# printf "Processing %s\n", $prefix;
  return 1;
}

sub loadMolecule {
  my $charmm=shift;
  my $prefix=shift;
  my $fname=shift;
  my $model=shift;

  my $Hfname="tmp$$.pdb";

#  printf "Creating clean pdb file";
#  printf " for model %d",$model if (defined $model);
#  printf "\n";
  my $tmol=Molecule::new();
  $tmol->readPDB($fname,model=>$model);
  $tmol->removeHetero();
  $tmol->fixHistidine();
  $tmol->translate("CHARMM19");
  $tmol->generateSegNames();
  $tmol->renumber(1);
  $charmm->setupFromMolecule($tmol);
  $charmm->verbose("print coor");
  $charmm->orient();
  $charmm->verbose("cons fix sele .not. type H end");
  $charmm->minimizeSD;
  $charmm->writePDB($Hfname);
  $charmm->verbose("delete atom sele all end");

  my $mol=Molecule::new($Hfname);
  $mol->fixHistidine();
  $mol->translate("CHARMM19");
  $mol->generateSegNames();
  $mol->selectChain($mol->getChain());
  system("/bin/rm ".$Hfname);

  return $mol;
}

sub generateQlist {
  my $mol=shift;

  my $ra=$mol->{defchain}->{res};
  my $numres=$#{$ra}+1;

  my $Hlist=&_findHbonds($mol);
  my $Slist=&_findSidechainContacts($mol);

  my $qlist;
  foreach my $contact ( @{$Slist} ) {
    my $i=$contact->{i};
    my $j=$contact->{j};
    $qlist->{$i}->{$j}->{E}=$contact->{E};
    $qlist->{$i}->{$j}->{type}="S";
  }

  my $not_placed=();
  foreach my $contact ( @{$Hlist} ) {
    my $i=$contact->{i};
    my $j=$contact->{j};
    if (! defined $qlist->{$i}->{$j}) {
      $qlist->{$i}->{$j}->{E}=-1.0*$par{H_factor};
      $qlist->{$i}->{$j}->{type}="H";
      $contact->{n}--;
      push(@{$not_placed},$contact) if ($contact->{n} > 0);
    } else {
      push(@{$not_placed},$contact);
    }
  }

  foreach my $contact ( @{$not_placed} ) {
    my $i=$contact->{i};
    my $j=$contact->{j};

    $qlist->{$i}->{$j}->{type}.="H";
    $qlist->{$i}->{$j}->{type}.="H" if ($contact->{n} == 2);

    my $numorie=0;
    $numorie++ if ($i > 0);
    $numorie++ if ($j < $#{$ra});
    $numorie+=2 if (($j-$i) >= 4);
    my $Eorie=-1.0*$par{H_factor}*$contact->{n}/$numorie;
    if ($i > 0) {
      $qlist->{$i-1}->{$j}->{type}.="O";
      $qlist->{$i-1}->{$j}->{type}.="O" if ($contact->{n} == 2);
      $qlist->{$i-1}->{$j}->{E}+=$Eorie;
    }
    if ($j < $#{$ra}) {
      $qlist->{$i}->{$j+1}->{type}.="O";
      $qlist->{$i}->{$j+1}->{type}.="O" if ($contact->{n} == 2);
      $qlist->{$i}->{$j+1}->{E}+=$Eorie;
    }
    if (($j-$i) >= 4) {
      $qlist->{$i}->{$j-1}->{type}.="O";
      $qlist->{$i}->{$j-1}->{type}.="O" if ($contact->{n} == 2);
      $qlist->{$i}->{$j-1}->{E}+=$Eorie;
      $qlist->{$i+1}->{$j}->{type}.="O";
      $qlist->{$i+1}->{$j}->{type}.="O" if ($contact->{n} == 2);
      $qlist->{$i+1}->{$j}->{E}+=$Eorie;
    }
  }

  return $qlist;
}

sub aveQlist {
  my $qarr=shift;

  my $models=$#{$qarr}+1;
  my $qlist={};
  while ($#{$qarr} > -1) {
    my $l=shift(@{$qarr});
    foreach my $i ( keys %{$l} ) {
      foreach my $j ( keys %{$l->{$i}} ) {
        $qlist->{$i}->{$j}->{type}.=$l->{$i}->{$j}->{type};
        $qlist->{$i}->{$j}->{E}+=$l->{$i}->{$j}->{E};
      }
    }
  }

  foreach my $i ( keys %{$qlist} ) {
    foreach my $j ( keys %{$qlist->{$i}} ) {
      $qlist->{$i}->{$j}->{E}/=$models;
    }
  }

  return $qlist;
}

sub writeGoPDB {
  my $prefix=shift;
  my $mol=shift;

  my $outcoor=&GenUtil::getOutputFile("GO_".$prefix.".pdb");

  my $line="HEADER  This pdb file contains the coordinates for a Go model of ".$prefix;
  printf $outcoor "%s\n", $line;

  my $c=$mol->{defchain};
  my $aa=$c->{atom};

  for (my $i=0; $i<=$#{$aa}; $i++) {
    if ($aa->[$i]->{atomname} eq "CA") {
      my $a=$aa->[$i];
      printf $outcoor "ATOM %6d  CA  ALA %5d", $a->{resnum}, $a->{resnum};
      printf $outcoor "%12.3f %7.3f ", $a->{xcoor}, $a->{ycoor};
      printf $outcoor "%7.3f  1.00  1.00      PROT\n", $a->{zcoor};
    }
  }

  printf $outcoor "END\n";
  undef $outcoor;
  return 1;

}

sub writeGoSeq {
  my $prefix=shift;
  my $mol=shift;

  my $ra=$mol->{defchain}->{res};
  my $numres=$#{$ra}+1;

  my $outseq=&GenUtil::getOutputFile("GO_".$prefix.".seq");
  printf $outseq "* This CHARMM .seq file describes a Go model of %s\n", $prefix;
  printf $outseq "*\n\n";
  printf $outseq "read sequence card\n";
  printf $outseq "* Sequence for Go model of %s\n", $prefix;
  printf $outseq "*\n";
  printf $outseq "%d\n", $numres;

  my $curr=1;

  while ($curr <= $numres) {
    my @a=();
    for (my $i=0; $i < 16; $i++) {
      push(@a,$curr) if ($curr <= $numres);
      $curr++;
    }
    printf $outseq "G%s\n", join(" G",@a);
  }

  printf $outseq "\ngenerate prot setup\n";
  printf $outseq "rename resname ala sele all end\n";
  undef $outseq;
  return 1;
}

sub writeGoTop {
  my $prefix=shift;
  my $mol=shift;

  my $ra=$mol->{defchain}->{res};
  my $numres=$#{$ra}+1;

  my $outtop=&GenUtil::getOutputFile("GO_".$prefix.".top");
  printf $outtop "* This CHARMM .top file describes a Go model of %s\n", $prefix;
  printf $outtop "*\n\n";
  printf $outtop "read rtf card\n";
  printf $outtop "* Topology for Go model of %s\n", $prefix;
  printf $outtop "*\n";
  printf $outtop "   20   1\n";
  my $curr=1;
  foreach my $r ( @{$ra} ) {
    printf $outtop "MASS %-3d G%-3d     %f\n", $curr, $curr, &_resMass($r->{name});
    $curr++;
  }
  printf $outtop "\nDECL +CA\n\n";
  printf $outtop "AUTOGENERATE ANGLES DIHEDRAL\n\n";
  for (my $i=1; $i<=$numres; $i++) {
    printf $outtop "RESI G%-3d       0.0\n", $i;
    printf $outtop "GROU\n";
    printf $outtop "Atom  CA  G%-3d     0.0\n", $i;
    printf $outtop "Bond CA +CA\n\n";
  }
  printf $outtop "END\n\n";

  undef $outtop;
  return 1;
}

sub writeGoParam {
  my $prefix=shift;
  my $mol=shift;
  my $qlist=shift;
  my $Tf=shift;
  my $TfScaling=shift;

  my $c=$mol->{defchain};
  my $aa=$c->{atom};
  my $ra=$c->{res};
  my $numres=$#{$ra}+1;

  my $outparam=&GenUtil::getOutputFile("GO_".$prefix.".param");
  printf $outparam "* This CHARMM .param file describes a Go model of %s\n", $prefix;

  my $Etot_curr = 0.0;
  foreach my $i ( keys %{$qlist} ) {
    foreach my $j ( keys %{$qlist->{$i}} ) {
      $Etot_curr += $qlist->{$i}->{$j}->{E};
    }
  }

  if (! defined $TfScaling) {
    die "Either Tf or TfScaling must be defined for writeGoParam" if (! defined $Tf);
    printf $outparam "* Parameters scaled to obtain a Tf of: %f\n", $Tf;
    $TfScaling = -1.0 * $par{TfScaling} * $Tf * $numres / $Etot_curr;
    printf $outparam "* Scaling all native contacts by a factor of %f\n*\n\n", $TfScaling;
  } else {
    printf $outparam "* Parameters scaled based on a selection of: %f\n*\n\n", $TfScaling;
  }
#  printf "Scaling all native contacts by a factor of %f\n", $TfScaling;
  my $Etot = $TfScaling * $Etot_curr;
  my $EperRes = $Etot/$numres;

  my $Kb = -1.0 * $par{BondFactor} * $EperRes;
  my $Ka = -1.0 * $par{AngleFactor} * $EperRes;
  my $Kd = -1.0 * $par{DiheFactor} * $EperRes;
  my $Erep = $par{RepFactor} * $EperRes;

  printf $outparam "read param card\n";
  printf $outparam "* Parameters for Go model of %s\n*\n", $prefix;

  printf $outparam "\nBOND\n";
  my $currcoor=&_findCA($c,0);
  for (my $i=1; $i < $numres; $i++) {
    my $nextcoor=&_findCA($c,$i);
    printf $outparam "G%-3d    G%-3d      %f  %f\n",$i,$i+1,$Kb,sqrt(&_sqdist($currcoor,$nextcoor));
    $currcoor=$nextcoor;
  }

  printf $outparam "\nANGLE\n";
  my $A=&_findCA($c,0);
  my $B=&_findCA($c,1);
  for (my $i=2; $i < $numres; $i++) {
    my $C=&_findCA($c,$i);
    printf $outparam "G%-3d    G%-3d    G%-3d      %f  %f\n",$i-1,$i,$i+1,$Ka,&_getTheta($A,$B,$C);
    $A=$B;
    $B=$C;
  }

  printf $outparam "\nDIHEDRAL\n";
  my $all_dihe_params=&_loadDiheParams($par{DiheFile});
  my $prevname=$ra->[1]->{name};
  for (my $i = 2; $i < ($numres-1); $i++) {
    my $newname=$ra->[$i]->{name};
    my $dihe_params=$all_dihe_params->{&_oneFromThree($prevname)}->{&_oneFromThree($newname)};
    foreach my $period ( sort {$a <=> $b} keys %{$dihe_params} ) {
      printf $outparam "G%-3d G%-3d G%-3d G%-3d   %f  %d  %f\n",$i-1,$i,$i+1,$i+2,
        $dihe_params->{$period}->{K}*$Kd,$period,$dihe_params->{$period}->{min};
    }
    $prevname=$newname;
  }

  printf $outparam "\nNONBONDED NBXMOD 3 ATOM CDIEL SHIFT VATOM VDISTANCE VSWITCH -\n";
  printf $outparam "  CUTNB 399.0 CTOFNB 398.5 CTONNB 395.5 EPS 1.0 WMIN 1.5\n\n";
  for (my $i=0; $i < $numres; $i++) {
    my $atom=&_findCA($c,$i);
    my $closest = 999.9;
    for (my $j=0; $j < $numres; $j++) {
      if ((abs($j - $i) > 3) && (!defined $qlist->{$i}->{$j}) && (!defined $qlist->{$j}->{$i})) {
        my $comp = &_sqdist($atom,&_findCA($c,$j));
        $closest = $comp if ($comp < $closest);
      }
    }
    printf $outparam "G%-3d     0.0  %f  %f\n",$i+1,$Erep,sqrt($closest)*$par{RepSigFactor};
  }

  printf $outparam "\nNBFIX\n";
  foreach my $i ( sort {$a <=> $b} keys %{$qlist} ) {
    my $atom=&_findCA($c,$i);
    foreach my $j ( sort {$a <=> $b} keys %{$qlist->{$i}} ) {
      printf $outparam "G%-3d    G%-3d       %f    %f\n",$i+1,$j+1,
        $TfScaling*$qlist->{$i}->{$j}->{E},sqrt(&_sqdist($atom,&_findCA($c,$j)));
    }
  }

  printf $outparam "\nEND\n\n";
  undef $outparam;
  return 1;
}

sub writeQdetails {
  my $mol=shift;
  my $list=shift;
  my $fname=shift;

  my $ra=$mol->{defchain}->{res};
  my $numres=$#{$ra}+1;

  my $outdet=&GenUtil::getOutputFile($fname);
  foreach my $i ( sort {$a <=> $b} keys %{$list} ) {
    my $resA=$ra->[$i];
    foreach my $j ( sort {$a <=> $b} keys %{$list->{$i}} ) {
      my $resB=$ra->[$j];
      printf $outdet "%s %-3d %s %-3d  %5f %s\n",$resA->{name},$resA->{num},
        $resB->{name},$resB->{num},
        $list->{$i}->{$j}->{E},$list->{$i}->{$j}->{type};
    }
  }

  undef $outdet;
  return 1;
}

sub writeQlist {
  my $mol=shift;
  my $list=shift;
  my $fname=shift;

  my $ra=$mol->{defchain}->{res};

  my $outlist=&GenUtil::getOutputFile($fname);
  foreach my $i ( sort {$a <=> $b} keys %{$list} ) {
    my $resA=$ra->[$i];
    foreach my $j ( sort {$a <=> $b} keys %{$list->{$i}} ) {
      if (($list->{$i}->{$j}->{type} =~ "H") || 
                        ($list->{$i}->{$j}->{type} =~ "S")) {
	printf $outlist "%d %d\n", $resA->{num}, $ra->[$j]->{num};
      }
    }
  }

  undef $outlist;
  return 1;
}

sub _findHbonds {
  my $mol=shift;

  my $Hlist=();

  my $c=$mol->{defchain};
  for (my $i=0; $i < $#{$c->{res}}; $i++) {
    for (my $j=$i+3; $j <= $#{$c->{res}}; $j++) {
      my $hbs=&_testHbond($c,$i,$j);
      if ($hbs > 0) {
        my $Hbond={};
	$Hbond->{i}=$i;
	$Hbond->{j}=$j;
	$Hbond->{n}=$hbs;
	push(@{$Hlist}, $Hbond);
      }
    }
  }
  return $Hlist;
}

sub _testHbond {
  my $c=shift;
  my $i=shift;
  my $j=shift;

  return 0 if ((! defined $c->{res}->[$i]) || (! defined $c->{res}->[$j]));

  my $aa=$c->{atom};

  my $icoor={};
  my $ri=$c->{res}->[$i];
  for (my $an=$ri->{start}; $an<=$ri->{end}; $an++) {
    my $atomname=$aa->[$an]->{atomname};
    if (($atomname eq "N") || ($atomname eq "C") || ($atomname eq "O") || ($atomname eq "H")) {
      $icoor->{$atomname}=$aa->[$an];
    }
  }

  my $jcoor={};
  my $rj=$c->{res}->[$j];
  for (my $an=$rj->{start}; $an<=$rj->{end}; $an++) {
    my $atomname=$aa->[$an]->{atomname};
    if (($atomname eq "N") || ($atomname eq "C") || ($atomname eq "O") || ($atomname eq "H")) {
      $jcoor->{$atomname}=$aa->[$an];
    }
  }

  my $Hbonds=0;

  if ((defined $icoor->{N}) && (defined $icoor->{H}) && 
      (defined $jcoor->{C}) && (defined $jcoor->{O})) {
    my $rON=&_sqdist($icoor->{N},$jcoor->{O});
    if ($rON < $par{Hmaxdist_sq}) {
      $rON=sqrt($rON);
      my $rCH = sqrt(&_sqdist($icoor->{H},$jcoor->{C}));
      my $rOH = sqrt(&_sqdist($icoor->{H},$jcoor->{O}));
      my $rCN = sqrt(&_sqdist($icoor->{N},$jcoor->{C}));
      my $E = (1 / $rON) + (1 / $rCH) - (1 / $rOH) - (1 / $rCN);
      $E *= $par{KSfac};
      if ($E < $par{HBcut}) {
	$Hbonds++;
#	printf "%d %d\n", $c->{res}->[$i]->{num}, $c->{res}->[$j]->{num};
      }
    }
  }

  if ((defined $jcoor->{N}) && (defined $jcoor->{H}) && 
      (defined $icoor->{C}) && (defined $icoor->{O})) {
    my $rON=&_sqdist($jcoor->{N},$icoor->{O});
    if ($rON < $par{Hmaxdist_sq}) {
      $rON=sqrt($rON);
      my $rCH = sqrt(&_sqdist($jcoor->{H},$icoor->{C}));
      my $rOH = sqrt(&_sqdist($jcoor->{H},$icoor->{O}));
      my $rCN = sqrt(&_sqdist($jcoor->{N},$icoor->{C}));
      my $E = (1 / $rON) + (1 / $rCH) - (1 / $rOH) - (1 / $rCN);
      $E *= $par{KSfac};
      if ($E < $par{HBcut}) {
          $Hbonds++;
#	  printf "%d %d\n", $c->{res}->[$j]->{num}, $c->{res}->[$i]->{num};
      }
    }
  }

  return $Hbonds;
}

sub _findSidechainContacts {
  my $mol=shift;

  my $A=0;
  my $B=0;
  my $MJparams;
  my $fname=&GenUtil::getInputFile($par{MJfile});
  while (<$fname>) {
    my $inline = $_;
    chomp($inline);
    $MJparams->{$A}->{$B} = $inline;
    $MJparams->{$B}->{$A} = $inline;
    $B++;
    if ($B == 20) {
      $A++;
      $B=$A;
    }
  }
  undef $fname;

  my $c=$mol->{defchain};
  my $Slist=();
  for (my $i=0; $i < $#{$c->{res}}; $i++) {
    for (my $j=$i+3; $j <= $#{$c->{res}}; $j++) {
      my $sc=&_testScontact($c,$i,$j);
      if ($sc > 0) {
        my $SC={};
	$SC->{i}=$i;
	$SC->{j}=$j;
	$SC->{E}=-1.0*$MJparams->{&_MJcode($c->{res}->[$i]->{name})}->{&_MJcode($c->{res}->[$j]->{name})}
                   *$par{MJ_factor}/$par{MJ_ave};
	push(@{$Slist}, $SC);
      }
    }
  }
  return $Slist;
}

sub _testScontact {
  my $c=shift;
  my $i=shift;
  my $j=shift;

  return 0 if ((! defined $c->{res}->[$i]) || (! defined $c->{res}->[$j]));

  my $aa=$c->{atom};
  my $ir=$c->{res}->[$i];
  my $jr=$c->{res}->[$j];

  return 0 if (&_sqdist(&_findCA($c,$i),&_findCA($c,$j)) > $par{SC_CA_dist_sq});

  for (my $ai=$ir->{start}; $ai<=$ir->{end}; $ai++) {
    if (!$aa->[$ai]->{hyd}) {
      my $iname=$aa->[$ai]->{atomname};
      if (($iname ne "N") && ($iname ne "C") && ($iname ne "O") && ($iname ne "CA")) {
	my $icoor=$aa->[$ai];
        for (my $aj=$jr->{start}; $aj<=$jr->{end}; $aj++) {
          if (!$aa->[$aj]->{hyd}) {
            my $jname=$aa->[$aj]->{atomname};
            if (($jname ne "N") && ($jname ne "C") && ($jname ne "O") && ($jname ne "CA")) {
	      my $jcoor=$aa->[$aj];
	      return 1 if (&_sqdist($icoor,$jcoor) < $par{SC_dist_sq});
	    }
	  }
        }
      }
    }
  }
  return 0;

}

sub _findCA {
  my $c=shift;
  my $i=shift;

  my $atom={};
  my $aa=$c->{atom};
  my $ri=$c->{res}->[$i];
  for (my $an=$ri->{start}; $an<=$ri->{end}; $an++) {
    if ($aa->[$an]->{atomname} eq "CA") {
      $atom=$aa->[$an];
    }
  }

  return $atom;
}

sub _sqdist {
  my $i=shift;
  my $j=shift;

  my $dist = ($i->{xcoor} - $j->{xcoor}) * ($i->{xcoor} - $j->{xcoor});
  $dist += ($i->{ycoor} - $j->{ycoor}) * ($i->{ycoor} - $j->{ycoor});
  $dist += ($i->{zcoor} - $j->{zcoor}) * ($i->{zcoor} - $j->{zcoor});
  return $dist;
}

sub _getTheta {
  my $i=shift;
  my $j=shift;
  my $k=shift;

  my $ax = $i->{xcoor}-$j->{xcoor};
  my $bx = $k->{xcoor}-$j->{xcoor};
  my $ay = $i->{ycoor}-$j->{ycoor};
  my $by = $k->{ycoor}-$j->{ycoor};
  my $az = $i->{zcoor}-$j->{zcoor};
  my $bz = $k->{zcoor}-$j->{zcoor};
  my $dot_prod = ($ax * $bx) + ($ay * $by) + ($az * $bz);
  my $adist = sqrt(($ax * $ax) + ($ay * $ay) + ($az * $az));
  my $bdist = sqrt(($bx * $bx) + ($by * $by) + ($bz * $bz));
  my $cos_thet = $dot_prod / ($adist * $bdist);
  my $theta = atan2(sqrt(1 - $cos_thet * $cos_thet), $cos_thet) * 180 / 3.1415;

  return $theta;
}

sub _getRg {
  my $mol=shift;

  my $aa=$mol->{defchain}->{atom};

  my $cx=0.0;
  my $cy=0.0;
  my $cz=0.0;

  my $numCA=0;
  for (my $i=0; $i<=$#{$aa}; $i++) {
    if ($aa->[$i]->{atomname} eq "CA") {
      my $atom=$aa->[$i];
      $cx+=$atom->{xcoor};
      $cy+=$atom->{ycoor};
      $cz+=$atom->{zcoor};
      $numCA++;
    }
  }

  my $rg=0.0;
  my $cen={};
  $cen->{xcoor}=$cx/$numCA;
  $cen->{ycoor}=$cy/$numCA;
  $cen->{zcoor}=$cz/$numCA;

  for (my $i=0; $i<=$#{$aa}; $i++) {
    if ($aa->[$i]->{atomname} eq "CA") {
      $rg+=&_sqdist($aa->[$i],$cen);
    }
  }
  return sqrt($rg/$numCA);
}

sub _oneFromThree {
  my $three=shift;

  my %code = ('GLY', 'G', 'PRO', 'P', 'ALA', 'A', 'VAL', 'V',
	      'LEU', 'L', 'ILE', 'I', 'MET', 'M', 'CYS', 'C',
	      'PHE', 'F', 'TYR', 'Y', 'TRP', 'W', 'HIS', 'H',
	      'HSD', 'H', 'HSE', 'H', 'LYS', 'K', 'ARG', 'R',
	      'GLN', 'Q', 'ASN', 'N', 'GLU', 'E', 'ASP', 'D',
	      'SER', 'S', 'THR', 'T');

  return $code{uc($three)};

}

sub _resMass {
  my $three=shift;

  my %mass = ('GLY', 57, 'PRO', 97, 'ALA', 71, 'VAL', 99,
	      'LEU', 113, 'ILE', 113, 'MET', 131, 'CYS', 103,
	      'PHE', 147, 'TYR', 163, 'TRP', 186, 'HIS', 138,
	      'HSD', 138, 'HSE', 138, 'LYS', 128, 'ARG', 157,
	      'GLN', 128, 'ASN', 114, 'GLU', 128, 'ASP', 114,
	      'SER', 87, 'THR', 101);

  return $mass{uc($three)};

}

sub _MJcode {
  my $three=shift;

  my %code = ('GLY', 9, 'PRO', 19, 'ALA', 8, 'VAL', 5,
	      'LEU', 4, 'ILE', 3, 'MET', 1, 'CYS', 0,
	      'PHE', 2, 'TYR', 7, 'TRP', 6, 'HIS', 16,
	      'HSD', 16, 'HSE', 16, 'LYS', 18, 'ARG', 17,
	      'GLN', 13, 'ASN', 12, 'GLU', 15, 'ASP', 14,
	      'SER', 11, 'THR', 10);

  return $code{uc($three)};

}


sub _loadDiheParams {
  my $paramfile=shift;

  my $params={};

  my $fname=&GenUtil::getInputFile($paramfile);
  while (<$fname>) {
    my $inline = $_;
    chomp($inline);
    my @readline = split(' ', $inline);
    $params->{$readline[0]}->{$readline[1]}->{$readline[3]}->{K}=$readline[2];
    $params->{$readline[0]}->{$readline[1]}->{$readline[3]}->{min}=$readline[4];
  }

  undef $fname;
  return $params;
}


