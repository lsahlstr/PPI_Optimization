#######################################################################
# R script for implementing the Genetic Algorithm in the force-field 
# optimization of HA-antibody systems.
# Code written initially by Aaron Frank, U. Michigan, c. 2014
#######################################################################

#######################################################################
# Load required libraries. GA = genetic algorithm. plyr = package for
# splitting, applying, and combining data.
#######################################################################
library(GA)
library(plyr)

#######################################################################
# Determine i,j interaction energies: eps_ij = (eps_i * eps_j)**0.5
#######################################################################
combineEPS <- function(eps,mask) {
    mat <- NULL
    for (i in 1:20) {
        for (j in 1:20) {
            if (i <= j) {
                mat <- c(mat,sqrt(eps[i]*eps[j]))
            }
        }
    }
    return(mat)
}

combineEPSOutput <- function(eps) {
    mat <- res <- NULL
    for (i in 1:20){
        for (j in 1:20){
            if (i <= j){
                mat <- c(mat,sqrt(eps[i]*eps[j]))
                res <- c(res,paste(iparsRes[i],sep=""))
            }
        }
    }
    return(data.frame(res,mat))
}

# Make second mask - apply to "ij" interactions NOT being optimized
# when computing total energy
MakeMask <- function(mask) {
	mask2 <- NULL
	k <- 1
	for (i in 1:20){
        for (j in 1:20){
            if (i <= j){
				val <- mask[k]
				if (mask[k] == 0) {
					mask2[k] <- 1;
				} else {
					mask2[k] <- 0;
				}
				#print (c(k,mask[k],mask2[k]))
				k <- k + 1 
			}
        }
    }
    return(mask2)
}


#######################################################################
# Define fitness functions
#######################################################################
# Sum of Log Ranks (SLR) early recognition metric
# Ref: Venkatraman, J Chem. Info. Model. 2010; 50:2079
SLR <- function(result) {   
    ri <- which(result$flag==1)
    N <- nrow(result)
    i <-  1:length(ri)
    SLRmax <- -sum(log(i/N))
    return(-sum(log(ri/N))/SLRmax)
}

#fitness <- function(eps) {
#  tmp <- comb[,c(-1,-2,-3)]
#  eps <- combineEPS(eps)
#  eps <- eps*mask
#  eps <- matrix(eps,nrow=nrow(tmp),ncol=ncol(tmp),byrow=T)
#  tmp <- data.frame(flag=comb[,1],system=comb[,2],rmsd=comb[,3],ener=rowSums(eps*tmp))
#  tmp <- tmp[order(tmp$ener),]
#  mean(ddply(.data=tmp,.var=c("system"),.fun=SLR)$V1)
#  #+(cor(tmp$ener,tmp$rmsd))
#}

# Z-score
Zscore <- function(tmp) {
    ((mean(tmp$ener[tmp$flag==1])-mean(tmp$ener[tmp$flag==0]))/sd(tmp$ener[tmp$flag==0]))
}

calZ <- function(eps) {
  tmp <- comb[,c(-1,-2,-3)]
  eps <- combineEPS(eps)
  eps <- eps*mask
  eps <- matrix(eps,nrow=nrow(tmp),ncol=ncol(tmp),byrow=T)
  tmp <- data.frame(flag=comb[,1],system=comb[,2],rmsd=comb[,3],ener=rowSums(eps*tmp))
  mean(ddply(.data=tmp,.var=c("system"),.fun=Zscore)$V1)
}

fitnessZ <- function(eps) {
  tmp <- comb[,c(-1,-2,-3)]
  eps <- combineEPS(eps)
  eps <- eps*mask
  eps <- matrix(eps,nrow=nrow(tmp),ncol=ncol(tmp),byrow=T)
  tmp <- data.frame(flag=comb[,1],system=comb[,2],rmsd=comb[,3],ener=rowSums(eps*tmp))
  -1.0*mean(ddply(.data=tmp,.var=c("system"),.fun=Zscore)$V1)+(cor(tmp$ener,tmp$rmsd))+mean(ddply(.data=tmp,.var=c("system"),.fun=SLR)$V1)
}

calE <- function(eps) {
  tmp <- comb[,c(-1,-2,-3)]
  eps <- combineEPS(eps)
  #eps <- eps*mask
  eps <- matrix(eps,nrow=nrow(tmp),ncol=ncol(tmp),byrow=T)
  tmp <- data.frame(flag=comb[,1],system=comb[,2],rmsd=comb[,3],ener=rowSums(eps*tmp))
}

calE2 <- function(eps1,eps2) {
	# 1 = new; 2 = old params
  	tmp <- comb[,c(-1,-2,-3)]
  	eps1 <- combineEPS(eps1)
  	eps2 <- combineEPS(eps2)
  	eps1 <- eps1*mask
  	eps2 <- eps2*mask2
  	eps1 <- matrix(eps1,nrow=nrow(tmp),ncol=ncol(tmp),byrow=T)
  	eps2 <- matrix(eps2,nrow=nrow(tmp),ncol=ncol(tmp),byrow=T)
  	tmp <- data.frame(flag=comb[,1],system=comb[,2],rmsd=comb[,3],ener=(rowSums(eps1*tmp)+rowSums(eps2*tmp)))
}

#######################################################################
# Read data
#######################################################################
pdbs <- c("optfiles")
for (p in seq_along(pdbs)) {
    pdb <- pdbs[p]
    
    if (p == 1) {
    
    	# Residues and initial parameters
        ipars <- read.table(paste(pdb,"/eii.txt",sep=""))
        iparsRes <- ipars$V1
        ipars <- ipars$V2
        
        # Mask
        maskfile <- read.table(paste(pdb,"/mask.dat",sep=""))
        mask <- maskfile$V2
        
        mask2 <- MakeMask(mask)
        # mask2 <- as.numeric(mask == 0)
        # mask2 <- as.numeric(!(0 %in% mask))
        #stop("Exiting")
                
        # fij files
        native <- t(read.table(paste(pdb,"/n_fij.txt",sep="")))
        nonnative <- t(read.table(paste(pdb,"/nn_fij.txt",sep="")))
        
        # RMSD files
        nativeRMSD <- read.table(paste(pdb,"/n_rmsd.dat",sep=""))$V2
        nonnativeRMSD <- read.table(paste(pdb,"/nn_rmsd.dat",sep=""))$V2
            
        # Convert to matrices
        native <- data.matrix(native[2:nrow(native),])
        nonnative <- data.matrix(nonnative[2:nrow(nonnative),])
		native <- cbind(rep(1,nrow(native)),cbind(rep(p,nrow(native)),cbind(nativeRMSD,native)))
        nonnative <- cbind(rep(0,nrow(nonnative)),cbind(rep(p,nrow(nonnative)),cbind(nonnativeRMSD,nonnative)))
        comb <- as.matrix(rbind(native,nonnative))
        comb <- matrix(as.numeric(comb),ncol=ncol(comb),nrow=nrow(comb),byrow=F)
    	
    } else {
    
    	# fij files
        native <- t(read.table(paste(pdb,"/n_fij.txt",sep="")))
        nonnative <- t(read.table(paste(pdb,"/nn_fij.txt",sep="")))
        
        # RMSD files
        nativeRMSD <- read.table(paste(pdb,"/n_rmsd.dat",sep=""))$V2
        nonnativeRMSD <- read.table(paste(pdb,"/nn_rmsd.dat",sep=""))$V2
        
        # Convert to matrices
        native <- data.matrix(native[2:nrow(native),])
        nonnative <- data.matrix(nonnative[2:nrow(nonnative),])
        native <- cbind(rep(1,nrow(native)),cbind(rep(p,nrow(native)),cbind(nativeRMSD,native)))
        nonnative <- cbind(rep(0,nrow(nonnative)),cbind(rep(p,nrow(nonnative)),cbind(nonnativeRMSD,nonnative)))
        tmp <- as.matrix(rbind(native,nonnative))
        tmp <- matrix(as.numeric(tmp),ncol=ncol(tmp),nrow=nrow(tmp),byrow=F)
        comb <- rbind(comb,tmp)
        
    }
}

# Compute initial energy
check <- calE(ipars)
check$oldenergy <- calE(ipars)$ener
check <- check[order(check$system,check$ener),]
write.table(check,file="check1.txt",row.names=F,quote=F,col.names=T)

check <- calE2(ipars,ipars)
check$oldenergy <- calE2(ipars,ipars)$ener
check <- check[order(check$system,check$ener),]
write.table(check,file="check2.txt",row.names=F,quote=F,col.names=T)


#stop("Exiting")

# Create a test set
test <- comb[comb[,2]==3,]
test <- comb
comb <- comb[comb[,2]!=3,]


#######################################################################
# Genetic Algorithm
#######################################################################
# Setup parameters
popSize <- 10
iters <- 1000
min <- rep(-1,length(ipars))
max <- rep(0,length(ipars))
initialSolution <- matrix(ipars,ncol=length(ipars),nrow=popSize,byrow=F)

# GA to assign weights; fitness function = Z-score
GAReal <- ga(type = "real-valued", fitness=fitnessZ, min=min, max=max, popSize=popSize, maxiter=iters,suggestions=initialSolution,keepBest=T)

# Get value of cost function at each iteration
Z <- NULL
for (i in 1:iters){
    Z <- c(Z,calZ(as.vector(GAReal@bestSol[[i]][1,])))
}

# Summarize Results
Z <- data.frame(iters=1:iters,Z=Z)

train <- comb
comb <- test

bestPars <- as.vector(GAReal@bestSol[[iters]])
check <- calE2(bestPars,ipars)
check$oldenergy <- calE2(ipars,ipars)$ener

check <- check[order(check$system,check$ener),]
results <- data.frame(resname=iparsRes, eii_initial=ipars,eii_optimized=as.matrix(as.vector(GAReal@bestSol[[iters]])))

write.table(combineEPSOutput(as.vector(GAReal@bestSol[[iters]])),file="parametersGA.txt",row.names=F,quote=F,col.names=F)
write.table(results,file="oldNewComparison.txt",row.names=F,quote=F,col.names=T)
write.table(check,file="rmsdEnergyComparison.txt",row.names=F,quote=F,col.names=T)
write.table(Z,file="Zscore.txt",row.names=F,quote=F,col.names=T)