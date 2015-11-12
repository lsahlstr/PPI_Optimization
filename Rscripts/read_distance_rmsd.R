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
# Read distance and RMSD information
#######################################################################
readPool <- function(pool="n",lengthm) {
	
	superM <- matrix(ncol=lengthm)
		
	poolRMSD <- read.table(paste(pdb,"/",pool,"_rmsd.dat",sep=""))$V2

	for (conf in 1:length(poolRMSD)) {

		# read in data
		R <- read.table(paste(pdb,"/",pool,"_rij.",conf,".txt",sep=""),col.names=c("pair","rij"))

		# Get system info		
		info <- NULL
		info[1] <- 0
		if (pool == "n") {
			info[1] <- 1
		}
		info[2] <- p
		info[3] <- poolRMSD[conf]
		info[4] <- conf

		# order based on pair (resnames of the ij pair)
		R <- R[order(R$pair),]
		R$pair <- as.character(R$pair)
		
		# Reference with all 210 aimino acid pairs
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