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
	
	testidx_native <- sample(1:nrow(native_r12),2)
	testidx_nonnative <- sample(1:nrow(nonnative_r12),5)
	
	# Separate r^12 into training and testing sets
	native_r12_test <- native_r12[testidx_native,]
	native_r12 <- native_r12[-testidx_native,]
	
	nonnative_r12_test <- nonnative_r12[testidx_nonnative,]
	nonnative_r12 <- nonnative_r12[-testidx_nonnative,]
	
	if (p==1){
		r12 <- as.matrix(rbind(native_r12,nonnative_r12))
		r12 <- matrix(as.numeric(r12),ncol=ncol(r12),nrow=nrow(r12),byrow=F)
		
		r12_test <- as.matrix(rbind(native_r12_test,nonnative_r12_test))
		r12_test <- matrix(as.numeric(r12_test),ncol=ncol(r12_test),nrow=nrow(r12_test),byrow=F)
		
	} else {
		tmp <- as.matrix(rbind(native_r12,nonnative_r12))
        tmp <- matrix(as.numeric(tmp),ncol=ncol(tmp),nrow=nrow(tmp),byrow=F)
        r12 <- rbind(r12,tmp)
        
        tmp <- as.matrix(rbind(native_r12_test,nonnative_r12_test))
        tmp <- matrix(as.numeric(tmp),ncol=ncol(tmp),nrow=nrow(tmp),byrow=F)
        r12_test <- rbind(r12_test,tmp)
	}
	
	# rij^6
	native_r6 <- t(read.table(paste(pdb,"/n_r6ij.txt",sep="")))
	native_r6 <- data.matrix(native_r6[2:nrow(native_r6),])
	native_r6 <- cbind(rep(1,nrow(native_r6)),cbind(rep(p,nrow(native_r6)),cbind(nativeRMSD,native_r6)))      
	
	nonnative_r6 <- t(read.table(paste(pdb,"/nn_r6ij.txt",sep="")))
	nonnative_r6 <- data.matrix(nonnative_r6[2:nrow(nonnative_r6),])
	nonnative_r6 <- cbind(rep(0,nrow(nonnative_r6)),cbind(rep(p,nrow(nonnative_r6)),cbind(nonnativeRMSD,nonnative_r6)))
	
	# Separate r^6 into training and testing sets
	native_r6_test <- native_r6[testidx_native,]
	native_r6 <- native_r6[-testidx_native,]
	
	nonnative_r6_test <- nonnative_r6[testidx_nonnative,]
	nonnative_r6 <- nonnative_r6[-testidx_nonnative,]
	
	if (p==1){
		r6 <- as.matrix(rbind(native_r6,nonnative_r6))
		r6 <- matrix(as.numeric(r6),ncol=ncol(r6),nrow=nrow(r6),byrow=F)
		
		r6_test <- as.matrix(rbind(native_r6_test,nonnative_r6_test))
		r6_test <- matrix(as.numeric(r6_test),ncol=ncol(r6_test),nrow=nrow(r6_test),byrow=F)
		
	} else {
		tmp <- as.matrix(rbind(native_r6,nonnative_r6))
        tmp <- matrix(as.numeric(tmp),ncol=ncol(tmp),nrow=nrow(tmp),byrow=F)
        r6 <- rbind(r6,tmp)
        
        tmp <- as.matrix(rbind(native_r6_test,nonnative_r6_test))
        tmp <- matrix(as.numeric(tmp),ncol=ncol(tmp),nrow=nrow(tmp),byrow=F)
        r6_test <- rbind(r6_test,tmp)
	}
	
	if (potFlag == "eten") {
		# rij^10
		native_r10 <- t(read.table(paste(pdb,"/n_r10ij.txt",sep="")))
		native_r10 <- data.matrix(native_r10[2:nrow(native_r10),])
		native_r10 <- cbind(rep(1,nrow(native_r10)),cbind(rep(p,nrow(native_r10)),cbind(nativeRMSD,native_r10)))
		   
		nonnative_r10 <- t(read.table(paste(pdb,"/nn_r10ij.txt",sep="")))
		nonnative_r10 <- data.matrix(nonnative_r10[2:nrow(nonnative_r10),])
		nonnative_r10 <- cbind(rep(0,nrow(nonnative_r10)),cbind(rep(p,nrow(nonnative_r10)),cbind(nonnativeRMSD,nonnative_r10)))
	
		# Separate r^10 into training and testing sets
		native_r10_test <- native_r10[testidx_native,]
		native_r10 <- native_r10[-testidx_native,]
	
		nonnative_r10_test <- nonnative_r10[testidx_nonnative,]
		nonnative_r10 <- nonnative_r10[-testidx_nonnative,]
	
		if (p==1){
			r10 <- as.matrix(rbind(native_r10,nonnative_r10))
			r10 <- matrix(as.numeric(r10),ncol=ncol(r10),nrow=nrow(r10),byrow=F)
		
			r10_test <- as.matrix(rbind(native_r10_test,nonnative_r10_test))
			r10_test <- matrix(as.numeric(r10_test),ncol=ncol(r10_test),nrow=nrow(r10_test),byrow=F)
		
		} else {
			tmp <- as.matrix(rbind(native_r10,nonnative_r10))
			tmp <- matrix(as.numeric(tmp),ncol=ncol(tmp),nrow=nrow(tmp),byrow=F)
			r10 <- rbind(r10,tmp)
		
			tmp <- as.matrix(rbind(native_r10_test,nonnative_r10_test))
			tmp <- matrix(as.numeric(tmp),ncol=ncol(tmp),nrow=nrow(tmp),byrow=F)
			r10_test <- rbind(r10_test,tmp)
		}
	
		list <- list(r12=r12,r12_test=r12_test,r10=r10,r10_test=r10_test,r6=r6,r6_test=r6_test)
	
	} else {
		list <- list(r12=r12,r12_test=r12_test,r6=r6,r6_test=r6_test)
    }
    
    return(list)

}
