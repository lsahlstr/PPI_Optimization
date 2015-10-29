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
		help="Flag for optimization type: gentle (refinement) or stringent [default %default]")
)
parser <- OptionParser(usage = "%prog [options] epsilonFile rminFile maskFile potentialType", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

cat("\nProject: Energy Landscape Theory optimization of a coarse-grained protein-protein interaction force field.\n")
cat("Purpose: To determine a set of parameters (eps_ij and rmin_ij) that yield native-like interfaces as the lowest energy states.\n")
cat("Authors: Logan S. Ahlstrom and Aaron T. Frank\n")
cat(sprintf("%s\n\n",date()))


#######################################################################
# Source external routines
#######################################################################
# Combine energy and distance parameters
source('combine_pars.R')
# Fitness and energy functions
source('fitness_energy.R')
# Read disance and RMSD information
source('read_distance_rmsd.R')


#######################################################################
# Read data for each system
#######################################################################
pdbs <- read.table(opt$list)$V1
for (p in seq_along(pdbs)){
    pdb <- pdbs[p]
    
    cat(sprintf("%s\n",pdbs[p]))
    
    if (p==1){
        
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
        iparsRes <- eps_file$V1
                
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
        
        list <- readDistRMSDinfo(pdb,p,potFlag)
       
    } else {
		list <- readDistRMSDinfo(pdb,p,potFlag)    
    }
    
    r12 <- list$r12
    r12_test <- list$r12_test
    r6 <- list$r6
    r6_test <- list$r6_test
    
    if (potFlag == "eten") {
		r10 <- list$r10
		r10_test <- list$r10_test
    }    
}

# Get flag, system, and rmsd info
# Training set: first three columns of r12
tmp_info <- data.frame(flag=r12[,1],system=r12[,2],rmsd=r12[,3])

# Testing set: first three columns of r12_test
tmp_info_test <- data.frame(flag=r12_test[,1],system=r12_test[,2],rmsd=r12_test[,3])


# Get size (m x n) for defining eps and rmin matrices in energy and fitness function subroutines
# For training set
tmp_r12 <- r12[,c(-1,-2,-3)]
m <- nrow(tmp_r12)
n <- ncol(tmp_r12)

# For testing set
tmp_r12_test <- r12_test[,c(-1,-2,-3)]
m_test <- nrow(tmp_r12_test)
n_test <- ncol(tmp_r12_test)


#######################################################################
# Genetic Algorithm optimization
#######################################################################
# Setup parameters
popSize <- opt$popsize
iters <- opt$ncycles

# Define min and max
min_eps <- NULL
max_eps <- NULL
min_rmin <- NULL
max_rmin <- NULL

len <- length(ipars)
for (i in 1:len) {
	# eps
	if (i <= (len/2)) {	
		if (opttypeFlag == "gentle") {
			min_eps[i] <- (ipars[i] - (0.5*ipars[i])) # opt$minval
			max_eps[i] <- (ipars[i] + (0.5*ipars[i])) # opt$maxval
		} else {
			min_eps[i] <- opt$minval
			max_eps[i] <- opt$maxval
		}
	# rmin		
	} else {
		min_rmin[(i-(len/2))] <- (ipars[i] - (1.0*rminSD[(i-(len/2))]))
		max_rmin[(i-(len/2))] <- (ipars[i] + (1.0*rminSD[(i-(len/2))]))
	}
}
min <- c(min_eps,min_rmin) 
max <- c(max_eps,max_rmin)
	
# Initial solution
initialSolution <- matrix(ipars,ncol=length(ipars),nrow=popSize,byrow=T)  # ipars_mask

# GA to assign weights; fitness function = Z-score for single system
ffunc <- NULL
if (potFlag == "eten") {
	if (fitFlag == "zscore") {
		ffunc <- fitnessZ_eten
    } else {
    	ffunc <- fitnessSLR_eten
	}
} else {
	if (fitFlag == "zscore") {
		ffunc <- fitnessZ_lj
    } else {
    	ffunc <- fitnessSLR_lj
	}
}

# Run GA
cat(sprintf("%s\n\n",date()))
GAReal <- ga(type = "real-valued", fitness=ffunc, min=min, max=max, popSize=popSize, maxiter=iters, suggestions=jitter(initialSolution), keepBest=T, parallel=TRUE) 
cat(sprintf("%s\n\n",date()))

# Save data from all GA cycles to R data structure
save(GAReal,file="GA.RData")

# Best solution from GA optimization
bestPars <- as.vector(GAReal@bestSol[[iters]][1,])


#######################################################################
# Get the mean Z-score and SLR at each iteration
#######################################################################
zscore <- NULL
slr <- NULL
if (potFlag == "eten") {
	for (i in 1:iters) {
		zscore <- c(zscore,calZ_eten(as.vector(GAReal@bestSol[[i]][1,])))
		slr <- c(slr,calSLR_eten(as.vector(GAReal@bestSol[[i]][1,])))
	}
} else {
	for (i in 1:iters) {
		zscore <- c(zscore,calZ_lj(as.vector(GAReal@bestSol[[i]][1,])))
		slr <- c(slr,calSLR_lj(as.vector(GAReal@bestSol[[i]][1,])))
	}
}

zscore <- data.frame(iters=1:iters,score=zscore)
slr <- data.frame(iters=1:iters,score=slr)


#######################################################################
# Z-score and SLR for each system after final optimization cycle
#######################################################################
zscore_indiv_old <- NULL
zscore_indiv_new <- NULL
slr_indiv_old <- NULL
slr_indiv_new <- NULL

if (potFlag == "eten") {

	zscore_indiv_old <- calZ_eten_indiv(ipars)
	zscore_indiv_new <- calZ_eten_indiv(bestPars)
	zscore_indiv_write <- data.frame(old=zscore_indiv_old,new=zscore_indiv_new)
	zscore_indiv_write$pdb <- pdbs
	
	slr_indiv_old <- calSLR_eten_indiv(ipars)
	slr_indiv_new <- calSLR_eten_indiv(bestPars)
	slr_indiv_write <- data.frame(old=slr_indiv_old,new=slr_indiv_new)
	slr_indiv_write$pdb <- pdbs

} else {

	zscore_indiv_old <- calZ_lj_indiv(ipars)
	zscore_indiv_new <- calZ_lj_indiv(bestPars)
	zscore_indiv_write <- data.frame(old=zscore_indiv_old,new=zscore_indiv_new)
	zscore_indiv_write$pdb <- pdbs

	slr_indiv_old <- calSLR_lj_indiv(ipars)
	slr_indiv_new <- calSLR_lj_indiv(bestPars)
	slr_indiv_write <- data.frame(old=slr_indiv_old,new=slr_indiv_new)
	slr_indiv_write$pdb <- pdbs
	
}


#######################################################################
# Compute energies with initial ("old") and optimized ("new") parameters
# and combine with RMSD and system information. Organize energy and RMSD 
# information and old and new parameters, and then write to output files.
#######################################################################
# RMSD and energy
if (potFlag == "eten") {
	
	rmsd_ener <- calE_eten(bestPars,r12,r10,r6,m,n,tmp_info)
	rmsd_ener$oldenergy <- calE_eten(ipars,r12,r10,r6,m,n,tmp_info)$ener
	
	rmsd_ener_test <- calE_eten(bestPars,r12_test,r10_test,r6_test,m_test,n_test,tmp_info_test)
	rmsd_ener_test$oldenergy <- calE_eten(ipars,r12_test,r10_test,r6_test,m_test,n_test,tmp_info_test)$ener
	
} else {

	rmsd_ener <- calE_lj(bestPars,r12,r6,m,n,tmp_info)
	rmsd_ener$oldenergy <- calE_lj(ipars,r12,r6,m,n,tmp_info)$ener
	
	rmsd_ener_test <- calE_lj(bestPars,r12_test,r6_test,m_test,n_test,tmp_info_test)
	rmsd_ener_test$oldenergy <- calE_lj(ipars,r12_test,r6_test,m_test,n_test,tmp_info_test)$ener
	
}

rmsd_ener <- rmsd_ener[order(rmsd_ener$system,rmsd_ener$ener),]
rmsd_ener_test <- rmsd_ener_test[order(rmsd_ener_test$system,rmsd_ener_test$ener),]

# Old and new parameters
old_new <- data.frame(respair=iparsRes, 
	eij_initial=ipars[1:(len/2)], 
	eij_optimized=bestPars[1:(len/2)],
	rij_initial=ipars[((len/2)+1):len],
	rij_optimized=bestPars[((len/2)+1):len])


#######################################################################
# Compute enrichment score based upon sorted energies and RMSD values
#######################################################################
enrich_write <- ddply(.dat=rmsd_ener,.var=c("system"),.fun=enrichment)


#######################################################################
# Write data to output
#######################################################################
write.table(old_new,file="OldNewComp.txt",row.names=F,quote=F,col.names=T)
write.table(rmsd_ener,file="rmsdEnerComp.txt",row.names=F,quote=F,col.names=F)
write.table(rmsd_ener_test,file="rmsdEnerComp_test.txt",row.names=F,quote=F,col.names=F)
write.table(zscore,file="zscore.dat",row.names=F,quote=F,col.names=F)
write.table(zscore_indiv_write,file="zscore_indiv.dat",row.names=F,quote=F,col.names=F)
write.table(slr,file="slr.dat",row.names=F,quote=F,col.names=F)
write.table(slr_indiv_write,file="slr_indiv.dat",row.names=F,quote=F,col.names=F)
write.table(enrich_write,file="enrich.dat",row.names=F,quote=F,col.names=F)


#######################################################################
# Save R environment variables and image
#######################################################################
save(list=c("initialSolution","ipars","bestPars","min","max","list","rmsd_ener","old_new"),file="check.RData")
save.image("image.RData")

cat("Done!\n")
cat(sprintf("%s\n\n",date()))