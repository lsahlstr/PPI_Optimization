#!/usr/bin/env bash

# Script variables
pool=$1
member=$2

# Environment variables
EXEC_DIR=/export/nVerde/users/lsahlstr/CGPPI_FF/opt6/bin
#CHARMMEXEC=/export/nVerde/users/lsahlstr/apps/charmm/exec/em64t_M/charmm
CHARMMEXEC=/export/nVerde/users/lsahlstr/apps/charmm_bigMAXATC/exec/em64t_M/charmm

# Make all-atom PDBs in CHARMM-friendly format
# Receptor from bound form (for native pose)
if ([ "$pool" == "native" ] && [ $member -eq 0 ]); then
	convpdb.pl -nsel protein -setchain A -segnames rec.0.pdb | sed 's/PRO0/PROA/g' | sed 's/HSD/HIS/g' | grep -v END > a_AA.pdb
	cp -p a_AA.pdb a_AA.0.pdb
fi
# Receptor from unbound form
if ([ "$pool" == "native" ] && [ $member -eq 1 ]); then
	convpdb.pl -nsel protein -setchain A -segnames rec.pdb | sed 's/PRO0/PROA/g' | sed 's/HSD/HIS/g' | grep -v END > a_AA.pdb
fi
# Ligand from unbound form
convpdb.pl -nsel protein -setchain B -segnames ${pool}/lig.${member}.pdb | sed 's/PRO0/PROB/g' | sed 's/HSD/HIS/g' > b_AA.pdb
cat a_AA.pdb b_AA.pdb > ab_AA.pdb

# Build the Go-like model
$EXEC_DIR/my-GoLikeBuilder.pl complex ab_AA.pdb

# Rename the Go-like model and remove interactions between chains
$EXEC_DIR/fixGOmodel.pl -gotag GO_complex -tag rename ab_AA.pdb

# Make Calpha-only PDB files, with Go-like model naming
grep 'PROA' rename.GO_complex.pdb | awk '{print}END{print "TER\nEND"}' > a.pdb
grep 'PROB' rename.GO_complex.pdb | awk '{print}END{print "TER\nEND"}' > b.pdb

# Append the generic inter-protein potential to the Go-like model
$EXEC_DIR/GoLike_all2all.pl -sig_file $EXEC_DIR/dists_respair_modGLY.txt -ene_file $EXEC_DIR/MJ_respair.txt -param_file rename.GO_complex.param -tag all2all -x_intra 2.0 -x_inter 0 -x_all2all 0.25 ab_AA.pdb

# Rename the parameter and topology files
mv all2all.rename.GO_complex.param go.param
mv rename.GO_complex.top go.top

# Minimize the complex
cp -p $EXEC_DIR/min_go.inp .
# minsteps=1000
$CHARMMEXEC snbf=200000 top=go.top param=go.param < min_go.inp | egrep '^MINI|^INTE>' > $pool/min.$member.log

# Split minimized structure into two components
grep PROA min.pdb | awk '{print}END{print "TER\nEND"}' > tmp.a.pdb
grep PROB min.pdb | awk '{print}END{print "TER\nEND"}' > tmp.b.pdb
mv init.pdb $pool/init.$member.pdb
mv min.pdb $pool/min.$member.pdb

# Template Calpha-only PDB files for Go model --> "real" residue naming
grep CA a_AA.pdb | awk '{print}END{print "TER\nEND"}' > tmp.template_a.pdb
grep CA b_AA.pdb | awk '{print}END{print "TER\nEND"}' > tmp.template_b.pdb

# PDB -> Traj -> PDB to get PDBs with "real" residue naming
pdb2traj.pl -out tmp.a.dcd tmp.a.pdb
traj2PDB -pdb tmp.template_a.pdb tmp.a.dcd
mv file.1.pdb tmp.a_resnames.pdb
 
pdb2traj.pl -out tmp.b.dcd tmp.b.pdb
traj2PDB -pdb tmp.template_b.pdb tmp.b.dcd
mv file.1.pdb tmp.b_resnames.pdb

$EXEC_DIR/PairDist.pl tmp.a_resnames.pdb tmp.b_resnames.pdb $EXEC_DIR/aa_pairs.txt

rm *GO_complex*
