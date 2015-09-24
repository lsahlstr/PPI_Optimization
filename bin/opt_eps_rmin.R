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
		help="One-column file listing the sub-directories for each system to be included in the optimization [default %default]")
	#make_option(c("-g", "--gentle"), type="character", default="T",
	#	help="True/False flag for gentle optimization [default %default]")
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
        
        # Potential
		potFlag <- opt$potential
		
		# Fitness function
		fitFlag <- opt$fitness
        
        list <- readDistRMSDinfo(pdb,p,potFlag)
       
    } else {
		list <- readDistRMSDinfo(pdb,p,potFlag)    
    }
    
    if (potFlag == "eten") {
		#r112 <- list$r112
		r12 <- list$r12
		r10 <- list$r10
		r6 <- list$r6
    } else {
    	r12 <- list$r12
		r6 <- list$r6
    }
    
}

# Get flag, system, and rmsd info
tmp_r12 <- r12[,c(-1,-2,-3)]
tmp_info <- data.frame(flag=r12[,1],system=r12[,2],rmsd=r12[,3])

# Get size (m x n) for defining eps and rmin matrices in energy and fitness function subroutines
m <- nrow(tmp_r12)
n <- ncol(tmp_r12)


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
	if (i <= (len/2)) {
		min_eps[i] <- (ipars[i] - (0.5*ipars[i])) # opt$minval
		max_eps[i] <- (ipars[i] + (0.5*ipars[i])) # opt$maxval
	} else {
		min_rmin[(i-(len/2))] <- (ipars[i] - (0.5*rminSD[(i-(len/2))]))
		max_rmin[(i-(len/2))] <- (ipars[i] + (0.5*rminSD[(i-(len/2))]))
	}
}
min <- c(min_eps,min_rmin) 
max <- c(max_eps,max_rmin)
	
# Initial solution
initialSolution <- matrix(ipars,ncol=length(ipars),nrow=popSize,byrow=T)  # ipars_mask

# GA to assign weights; fitness function = Z-score for single system
if (potFlag == "eten") {
	if (fitFlag == "zscore") {
    	GAReal <- ga(type = "real-valued", fitness=fitnessZ_eten, min=min, max=max, popSize=popSize, maxiter=iters, suggestions=initialSolution, keepBest=T)
    } else {
		GAReal <- ga(type = "real-valued", fitness=fitnessSLR_eten, min=min, max=max, popSize=popSize, maxiter=iters, suggestions=initialSolution, keepBest=T)
	}
    	
} else {
	if (fitFlag == "zscore") {
    	GAReal <- ga(type = "real-valued", fitness=fitnessZ_lj, min=min, max=max, popSize=popSize, maxiter=iters, suggestions=initialSolution, keepBest=T)
    } else {
		GAReal <- ga(type = "real-valued", fitness=fitnessSLR_lj, min=min, max=max, popSize=popSize, maxiter=iters, suggestions=initialSolution, keepBest=T)
	}
}

# Save data from 
save(GAReal,file="GA.RData")

#######################################################################
# Summarize results
#######################################################################
# Get value of Z-score at each iteration
zscore <- NULL
slr <- NULL
if (potFlag == "eten") {
	#if (fitFlag == "zscore") {
		for (i in 1:iters) {
			zscore <- c(zscore,calZ_eten(as.vector(GAReal@bestSol[[i]][1,])))
		}
	#} else {
		for (i in 1:iters) {
			slr <- c(slr,calSLR_eten(as.vector(GAReal@bestSol[[i]][1,])))
		}
	#}
} else {
	#if (fitFlag == "zscore") {
		for (i in 1:iters) {
			zscore <- c(zscore,calZ_lj(as.vector(GAReal@bestSol[[i]][1,])))
		}
	#} else {
		for (i in 1:iters) {
			slr <- c(slr,calSLR_lj(as.vector(GAReal@bestSol[[i]][1,])))
		}
	#}
}

zscore <- data.frame(iters=1:iters,score=zscore)
slr <- data.frame(iters=1:iters,score=slr)


#######################################################################
# Compare initial ("old") and optimized ("best") parameters and save results
#######################################################################
bestPars <- as.vector(GAReal@bestSol[[iters]][1,])

if (potFlag == "eten") {
	check <- calE_eten(bestPars)
	check$oldenergy <- calE_eten(ipars)$ener
} else {
	check <- calE_lj(bestPars)
	check$oldenergy <- calE_lj(ipars)$ener
}

check <- check[order(check$system,check$ener),]

results <- data.frame(respair=iparsRes, 
	eij_initial=ipars[1:(len/2)], 
	eij_optimized=bestPars[1:(len/2)],
	rij_initial=ipars[((len/2)+1):len],
	rij_optimized=bestPars[((len/2)+1):len])

save(list=c("ipars","bestPars","min","max","list","check","results"),file="check.RData")
write.table(results,file="InitOptComp.txt",row.names=F,quote=F,col.names=T)
write.table(check,file="rmsdEnerComp.txt",row.names=F,quote=F,col.names=F)
write.table(zscore,file="zscore.dat",row.names=F,quote=F,col.names=F)
write.table(slr,file="slr.dat",row.names=F,quote=F,col.names=F)

cat("Done!\n")
cat(sprintf("%s\n\n",date()))