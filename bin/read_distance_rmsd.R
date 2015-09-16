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
readDistRMSDinfo <- function(pdb,p,potFlag) {		

	# RMSD files
	nativeRMSD <- read.table(paste(pdb,"/n_rmsd.dat",sep=""))$V2
	nonnativeRMSD <- read.table(paste(pdb,"/nn_rmsd.dat",sep=""))$V2
        
	# rij^12
	native_r12 <- t(read.table(paste(pdb,"/n_r12ij.txt",sep="")))
	native_r12 <- data.matrix(native_r12[2:nrow(native_r12),])
	native_r12 <- cbind(rep(1,nrow(native_r12)),cbind(rep(p,nrow(native_r12)),cbind(nativeRMSD,native_r12)))
		   
	nonnative_r12 <- t(read.table(paste(pdb,"/nn_r12ij.txt",sep="")))
	nonnative_r12 <- data.matrix(nonnative_r12[2:nrow(nonnative_r12),])
	nonnative_r12 <- cbind(rep(0,nrow(nonnative_r12)),cbind(rep(p,nrow(nonnative_r12)),cbind(nonnativeRMSD,nonnative_r12)))
	
	if (p==1){
		r12 <- as.matrix(rbind(native_r12,nonnative_r12))
		r12 <- matrix(as.numeric(r12),ncol=ncol(r12),nrow=nrow(r12),byrow=F)
	} else {
		tmp <- as.matrix(rbind(native_r12,nonnative_r12))
        tmp <- matrix(as.numeric(tmp),ncol=ncol(tmp),nrow=nrow(tmp),byrow=F)
        r12 <- rbind(r12,tmp)
	}
	
	# rij^6
	native_r6 <- t(read.table(paste(pdb,"/n_r6ij.txt",sep="")))
	native_r6 <- data.matrix(native_r6[2:nrow(native_r6),])
	native_r6 <- cbind(rep(1,nrow(native_r6)),cbind(rep(p,nrow(native_r6)),cbind(nativeRMSD,native_r6)))      
	
	nonnative_r6 <- t(read.table(paste(pdb,"/nn_r6ij.txt",sep="")))
	nonnative_r6 <- data.matrix(nonnative_r6[2:nrow(nonnative_r6),])
	nonnative_r6 <- cbind(rep(0,nrow(nonnative_r6)),cbind(rep(p,nrow(nonnative_r6)),cbind(nonnativeRMSD,nonnative_r6)))
	
	if (p==1){
		r6 <- as.matrix(rbind(native_r6,nonnative_r6))
		r6 <- matrix(as.numeric(r6),ncol=ncol(r6),nrow=nrow(r6),byrow=F)
	} else {
		tmp <- as.matrix(rbind(native_r6,nonnative_r6))
        tmp <- matrix(as.numeric(tmp),ncol=ncol(tmp),nrow=nrow(tmp),byrow=F)
        r6 <- rbind(r6,tmp)
	}
	
	if (potFlag == "eten") {
		# rij^10
		native_r10 <- t(read.table(paste(pdb,"/n_r10ij.txt",sep="")))
		native_r10 <- data.matrix(native_r10[2:nrow(native_r10),])
		native_r10 <- cbind(rep(1,nrow(native_r10)),cbind(rep(p,nrow(native_r10)),cbind(nativeRMSD,native_r10)))
		   
		nonnative_r10 <- t(read.table(paste(pdb,"/nn_r10ij.txt",sep="")))
		nonnative_r10 <- data.matrix(nonnative_r10[2:nrow(nonnative_r10),])
		nonnative_r10 <- cbind(rep(0,nrow(nonnative_r10)),cbind(rep(p,nrow(nonnative_r10)),cbind(nonnativeRMSD,nonnative_r10)))
	
		if (p==1){
			r10 <- as.matrix(rbind(native_r10,nonnative_r10))
			r10 <- matrix(as.numeric(r10),ncol=ncol(r10),nrow=nrow(r10),byrow=F)
		} else {
			tmp <- as.matrix(rbind(native_r10,nonnative_r10))
			tmp <- matrix(as.numeric(tmp),ncol=ncol(tmp),nrow=nrow(tmp),byrow=F)
			r10 <- rbind(r10,tmp)
		}
	
		list <- list(r12=r12,r10=r10,r6=r6)
	
	} else {
		list <- list(r12=r12,r6=r6)
    }
    
    return(list)

}