#!/export/apps/R/3.1.2/bin/Rscript
#######################################################################
# R script for running a Genetic Algorithm-based optimization of a 
# coarse-grained protein-protein interaction force field. 
#		Logan S. Ahlstrom, Blair Whittington, and Aaron T. Frank, and Blair Whittington, U. Michigan c. 2015
#######################################################################


#######################################################################
# Load libraries
#######################################################################
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("GA"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("doParallel"))
suppressPackageStartupMessages(library("ptw"))

cat("Analyze GA optimization.")
cat("Authors: Logan S. Ahlstrom, Blair Whittington, and Aaron T. Frank\n")
cat(sprintf("%s\n\n",date()))


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
# Genetic Algorithm optimization
source(paste(workdir,'ga_opt.R',sep=""))
# Mean Z-score and SLR at each cycle
source(paste(workdir,'fitness_cycle.R',sep=""))
# Z-score and SLR for each system before and after optimization
source(paste(workdir,'fitness_before_after.R',sep=""))
# Compute energies with initial and optimized parameters; combine with RMSD and system information
source(paste(workdir,'rmsd_ener.R',sep=""))
# Enrichment score
source(paste(workdir,'enrich.R',sep=""))


#######################################################################
# Setup and GA data from optimization run
#######################################################################
load('init.RData')
load('GA.RData')

#######################################################################
# Analysis
#######################################################################
# Best solution from GA optimization
# GA parameters
len <- length(ipars)
popSize <- opt$popsize
iters <- opt$ncycles
eps_min_val <- opt$minval
eps_max_val <- opt$maxval
bestPars <- as.vector(GAReal@bestSol[[iters]][1,])

# Mean Z-score and SLR at each cycle
#fitness_cycle()

# Z-score and SLR for each system before and after optimization
fitness_before_after()

# Compute energies with initial and optimized parameters; combine with RMSD and system information
rmsd_ener_out <- rmsd_ener()

# Enrichment score
enrich <- ddply(.data=rmsd_ener_out,.var=c("system"),.fun=enrichment)

# Save R environment variables and image
save.image("analysis.RData")

cat("Done!\n")
cat(sprintf("%s\n\n",date()))