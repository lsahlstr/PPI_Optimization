#!/export/apps/R/3.1.2/bin/Rscript
#######################################################################
# R script for implementing the Genetic Algorithm in the force-field 
# optimization of coarse-grained protein-protein binding simulations.
# Energies and distances are optimized simultaneously. 
#
# 	Initial version of the code written to optimize eij parameters 
#	with different fitness functions
#		Aaron Frank, U. Michigan, c. 2015
#	Code further developed to optimize eij and rmin simultaneously 
#	and with different potential and fitness functions
#		Logan S. Ahlstrom, U. Michigan c. 2015
#######################################################################


#######################################################################
# Load required libraries. 
# 	GA = genetic algorithm
# 	plyr = package for splitting, applying, and combining data.
#######################################################################
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("GA"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("doParallel"))
suppressPackageStartupMessages(library("ptw"))


#######################################################################
# Command line arguments
#######################################################################
option_list <- list( 
	make_option(c("-e", "--eps"), type="character", default="eii.txt",
		help="File containing initial guesses of interaction strengths, epsilon_ij [default %default]"),
	make_option(c("-r", "--rmin"), type="character", default="rii.txt",
		help="File containing initial guesses of interaction distances, rmin_ij [default %default]"),
	make_option(c("-m", "--mask"), type="character", default="mask_20.txt",
		help="File containing 1's and 0's to specify which parameters to optimize and which to ignore [default %default]"),
	make_option(c("-w", "--weights"), type="character", default="weights.txt",
		help="File containing weights specifying the contribution of each system to the fitness function [default %default]"),
	make_option(c("-p", "--potential"), type="character", default="lj",
		help="String defining the potential function to use for optimization (lj or eten) [default %default]"),
	make_option(c("-f", "--fitness"), type="character", default="zscore",
		help="Fitness function (zscore or slr) for GA optimization [default %default]"),
	make_option(c("-s", "--popsize"), type="integer", default=10,
        help="Population size for GA optimization [default %default]"),
    make_option(c("-c", "--ncycles"), type="integer", default=10,
        help="Number of iterations for GA optimization [default %default]"),
    make_option(c("-u", "--maxval"), type="double", default=2,
        help="Upper boundary/maximum value for epsilon parameter [default %default]"),
    make_option(c("-o", "--minval"), type="double", default=0,
        help="Lower boundary/minimum value for epsilon parameter [default %default]"),
    make_option(c("-l", "--list"), type="character", default="pdblist",
		help="One-column file listing the sub-directories for each system to be included in the optimization [default %default]"),
	make_option(c("-t", "--opttype"), type="character", default="gentle",
		help="Flag for optimization type: gentle (refinement) or stringent [default %default]"),
	make_option(c("-d", "--rdatafile"), type="character", default="rij.RData",
		help="R data file containing data structure with rij and RMSD information [default %default]")
)
parser <- OptionParser(usage = "%prog [options] epsilonFile rminFile maskFile potentialType", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

cat("\nProject: Energy Landscape Theory optimization of a coarse-grained protein-protein interaction force field.\n")
cat("Purpose: To determine a set of parameters (eps_ij and rmin_ij) that yield native-like interfaces as the lowest energy states.\n")
cat("Authors: Logan S. Ahlstrom and Aaron T. Frank\n")
cat(sprintf("%s\n\n",date()))


#######################################################################
# Load R image file with distance and RMSD information; rij_data; rij_rmsd_data
#######################################################################
load(opt$rdatafile)

# number of instances of each i,j interaction in rij_rmsd_data 
lsize <- 1000

# Training set:
# Get system info: flag, system, rmsd, and conformer
sysinfo <- data.frame(flag=rij_rmsd_data[,1],system=rij_rmsd_data[,2],rmsd=rij_rmsd_data[,3],conf=rij_rmsd_data[,4])
# Just the distance information (rij)
rij <- rij_rmsd_data[,c(-1,-2,-3,-4)]


#######################################################################
# Source external routines
#######################################################################
# Directory with R scripts
workdir <- '~/repos/PPI_Optimization/Rscripts/'
# Combine energy and distance parameters
source(paste(workdir,'combine_pars.R',sep=""))
# Fitness functions
source(paste(workdir,'fitness.R',sep=""))
# Energy routine
source(paste(workdir,'ener.R',sep=""))
# Epsilon and Rmin data structures
source(paste(workdir,'big_pars.R',sep=""))
# Energy check
source(paste(workdir,'ener_check.R',sep=""))
# Fitness functions
source(paste(workdir,'fitness.R',sep=""))
# Evaluate fitness functions
source(paste(workdir,'fitness_eval.R',sep=""))
# Evaluate fitness functions for each individual system
source(paste(workdir,'fitness_eval_indiv.R',sep=""))


#######################################################################
# Read data
#######################################################################    
# List of PDB ID's for each system
pdbs <- read.table(opt$list)$V1

# Initial eps values
eps_file <- read.table(opt$eps)
eps <- eps_file$V2

# Initial rmin values
rmin_file <- read.table(opt$rmin)
rmin <- rmin_file$V2
rminSD <- rmin_file$V3
  
# Combine eps and rmin into a single data structure
ipars <- c(eps,rmin)

# Get residue pair names
respairs <- eps_file$V1
		
# Mask file
mask_file <- read.table(opt$mask)
mask <- mask_file$V2

# Weights
weights_file <- read.table(opt$weights)
weights <- weights_file$V2 

# Potential
potFlag <- opt$potential

# Fitness function
fitFlag <- opt$fitness

# Gentle or stringent optimization
opttypeFlag <- opt$opttype

# Energy test
if (potFlag == "lj") {
	check <- ener_lj(ipars)
} else if (potFlag == "eten"){
	check <- ener_eten(ipars)
} else if (potFlag == "etsr"){
	check <- ener_etsr(ipars)
}

ener_check(check,0.1)

#######################################################################
# Genetic Algorithm optimization
#######################################################################
# Setup parameters
len <- length(ipars)
popSize <- opt$popsize
iters <- opt$ncycles
eps_min_val <- opt$minval
eps_max_val <- opt$maxval

# Initial solution
initialSolution <- matrix(ipars,ncol=length(ipars),nrow=popSize,byrow=T)

# Run GA
source(paste(workdir,'ga_opt.R',sep=""))
GAReal <- ga_opt()

# Best solution from GA optimization
bestPars <- as.vector(GAReal@bestSol[[iters]][1,])


#######################################################################
# Analysis
#######################################################################
# Get the mean Z-score and SLR at each cycle
source(paste(workdir,'fitness_cycle.R',sep=""))
fitness_cycle()

# Z-score and SLR for each system before and after optimization
source(paste(workdir,'fitness_before_after.R',sep=""))
fitness_before_after()

# Compute energies with initial and optimized parameters; combine with RMSD and system information
source(paste(workdir,'rmsd_ener.R',sep=""))
rmsd_ener <- rmsd_ener()

# Enrichment score
source(paste(workdir,'enrich.R',sep=""))
test <- ddply(.dat=rmsd_ener,.var=c("system"),.fun=enrichment)

# Save R environment variables and image
#save(list=c("initialSolution","ipars","bestPars","min","max","list","rmsd_ener","old_new"),file="check.RData")
save.image("opt.RData")

cat("Done!\n")
cat(sprintf("%s\n\n",date()))