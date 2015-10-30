#######################################################################
# R code to read distance and RMSD information as part of a Genetic 
# Algorithm-based force-field optimization for coarse-grained protein-
# protein binding simulations. The distance and RMSD information for each
# system considered in the optimization needs to be stored in its own 
# subdirectory (e.g., "hdea_lj" for simulations of HdeA run using a Lennard
# Jones potential for the inter-protein interactions).
#
# Within these subdirectories, the distance information is the 1/r^12, 1/r^10, 
# and 1/r^6 terms for solving the standar 12-6 Lennard Jones ("LJ") or 
# 12-10-6 ("ETEN") potentials expressed in matrix form. They are named "n_r12ij.txt"
# and "nn_r12ij.txt" for the 1/r^12 term for the native and non-native pools,
# respectively. RMSD informatoin is stored as "n_rmsd.dat" and "nn_rmsd.dat" 
# within each of the system subdirectories.
#
#       The initial version of the routine was written for the optimization
#		of energy parameters with a Lennard Jones function 
#               Aaron Frank, U. Michigan, c. 2015
#       The subroutine was further developed to consider energies and distances
#		simultaneously and with the 12-10-6 potential function
#               Logan S. Ahlstrom, U. Michigan c. 2015
#######################################################################


#######################################################################
# Read distance and RMSD information
#######################################################################
readPool <- function(pool="n",lengthm) {
	
	superM <- matrix(ncol=lengthm)
		
	poolRMSD <- read.table(paste(pdb,"/",pool,"_rmsd.dat",sep=""))$V2

	for (conf in 1:length(poolRMSD)) {

		# read in data
		R <- read.table(paste(pdb,"/",pool,"_rij.",conf,".txt",sep=""),col.names=c("pair","rij"))

		# Get system info
		#if (conf == 1) {		
			info <- NULL
			info[1] <- 0
			if (pool == "n") {
				info[1] <- 1
			}
			info[2] <- p
			info[3] <- poolRMSD[conf]
			info[4] <- conf
		#}

		# order based on pair (resnames of the ij pair)
		R <- R[order(R$pair),]

		R$pair <- as.character(R$pair)
		names <- as.character(respairs)

		# Loop over ij keys and create list (and matrix)
		# initial list l and matrix m
		l <- list(NULL)
		m <- matrix(nrow=length(names),ncol=lsize)

		for (i in seq_along(names)) {
			name <- names[i]
			tmp <- c(subset(R,pair==name)$rij)
			nzeros <- lsize-length(tmp)

			if(!is.null(tmp)){
				tmp <- as.vector(padzeros(tmp, nzeros, side="right"))
			} else {
				tmp <- rep(0,lsize)
			}						
			l[[name]] <- tmp
			m[i,] <- tmp					
		}
		m <- c(info,as.vector(t(m)))
		superM <- rbind(superM,m)
	}
		
	return(superM[-1,])
	
}

#######################################################################
# Read distance and RMSD information
#######################################################################
readDistRMSDinfo <- function(pdbs,lsize=1000) {		

	lengthm <- 4 + (210 * lsize)
	superM2 <- matrix(ncol=lengthm)
	
	for (p in seq_along(pdbs)){
		pdb <- pdbs[p]
		superM2 <- rbind(superM2,readPool("n",lengthm))
		superM2 <- rbind(superM2,readPool("nn",lengthm))
		
	}
	
	return(superM2[-1,])

}