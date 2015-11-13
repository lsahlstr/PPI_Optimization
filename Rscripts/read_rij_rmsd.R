#!/export/apps/R/3.1.2/bin/Rscript
#######################################################################
# R script to read distance and RMSD information as part of a Genetic 
# Algorithm-based force-field optimization for coarse-grained protein-
# protein binding simulations. The distance and RMSD information for each
# system considered in the optimization needs to be stored in its own 
# subdirectory (e.g., with a name corresonding to the 4-letter PDB ID)
#
# For each binding pose included in the optimization procedure, the 
# inter-molecular pairwise distance information is stored in a single file with
# the following format:
#	#ResPair  # distance (r_ij)
#	LYSPRO	  19.895
#	LYSVAL	  19.917
#	GLYLYS	  16.950
#	LYSPHE	  19.580
#	...
# These files need to be named according to the following format:
# 	[pool]_rij.[pose].txt
#			pool => native ("n") or non-native ("nn")
#			pose => an integeger to delineate poses within a given pool 
#
# RMSD information is stored as "n_rmsd.dat" and "nn_rmsd.dat" 
# within each of the system subdirectories.
#
# Logan S. Ahlstrom and Aaron T. Frank, U. Michigan c. 2015
#######################################################################

#######################################################################
# Load required libraries
#######################################################################
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ptw"))

#######################################################################
# Command line arguments
#######################################################################
option_list <- list(
	make_option(c("-l", "--list"), type="character", default="pdblist",
		help="One-column file listing the sub-directories for each system to be included in the optimization [default %default]"),
	make_option(c("-p", "--pairs"), type="character", default="respairs.txt",
		help="File containing all 210 possible amino acid pairs in alphabetical order [default %default]"),
	make_option(c("-o", "--outname"), type="character", default="rij_data",
		help="Name for data structure in which rij and RMSD information is stored [default %default]")
)
parser <- OptionParser(usage = "%prog [options] pdblist respairs", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

#######################################################################
# Read distance and RMSD information
#######################################################################
source('~/repos/PPI_Optimization/Rscripts/format_rij_rmsd.R')
source('~/repos/PPI_Optimization/Rscripts/combine_rij_rmsd.R')
pdbs <- read.table(opt$list)$V1
respairs <- read.table(opt$pairs)$V1
rij_rmsd_data <- combineRijRMSD(pdbs,1000)

image_name <- paste(opt$outname,".RData",sep="")
save.image(image_name)
