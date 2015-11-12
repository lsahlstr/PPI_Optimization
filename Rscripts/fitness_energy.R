#######################################################################
# R code containing subroutines to evaluate the fitness and energy functions
# as part of a Genetic Algorithm-based force-field optimization for coarse-grained 
# protein-protein binding simulations.
#
# 	Initial versions of the routines were written to optimize energy parameters 
#	with different fitness functions
#		Aaron Frank, U. Michigan, c. 2015
#	Subroutines were further developed for the simultaneous optimization of energies 
#	and distances, and for optimization with different potential and fitness functions
#		Logan S. Ahlstrom, U. Michigan c. 2015
#######################################################################


#######################################################################
# Define and evaluate the fitness functions
#######################################################################
# Z-score
Zscore <- function(tmp) {
    ((mean(tmp$ener[tmp$flag==1])-mean(tmp$ener[tmp$flag==0]))/sd(tmp$ener[tmp$flag==0]))
    #((mean(tmp$ener[tmp$flag==1])-mean(tmp$ener[tmp$flag==0]))/sd(tmp$ener[tmp$flag==0])) + 1/sd(tmp$ener[tmp$flag==1])
}

# Sum of Logarithms of Ranks (SLR) early recognition metric
SLR <- function(result){
    ri <- which(result$flag==1)
    N <- nrow(result)
    i <-  1:length(ri)
    SLRmax <- -sum(log(i/N))
    return(-sum(log(ri/N))/SLRmax)
    #return(-sum(log(ri/N))/SLRmax + 1/sd(results$ener[tmp$flag==1]))
}


#######################################################################
# Routines for computing the fitness function during optimization
#######################################################################
## 12-6 potential ##
# Z-score
fitnessZ_lj <- function(pars) {

	list <- combinePARS(pars)
	eps <- list$eps
	eps <- matrix(eps,nrow=m,ncol=n,byrow=T)
	rmin <- list$rmin
	rmin <- matrix(rmin,nrow=m,ncol=n,byrow=T)
	
	r12min <- rmin**12
	tmp_r12 <- r12[,c(-1,-2,-3)]
	tmp1 <- r12min*tmp_r12
	
	r6min <- rmin**6
	tmp_r6 <- r6[,c(-1,-2,-3)]
	tmp2 <- r6min*tmp_r6
	
	tmp3 <- eps*(tmp1 - 2*tmp2)
		
	tmp <- tmp_info	
	tmp$ener <- rowSums(tmp3)
    
    tmp4 <- -1.0*weighted.mean(ddply(.data=tmp,.var=c("system"),.fun=Zscore)$V1,weights)
    
    print(tmp4)
    cat(sprintf("%s\n\n",tmp4))
    
    return(tmp4)
}

# SLR
fitnessSLR_lj <- function(pars) {

	list <- combinePARS(pars)
	eps <- list$eps
	eps <- matrix(eps,nrow=m,ncol=n,byrow=T)
	rmin <- list$rmin
	rmin <- matrix(rmin,nrow=m,ncol=n,byrow=T)
	
	r12min <- rmin**12
	tmp_r12 <- r12[,c(-1,-2,-3)]
	tmp1 <- r12min*tmp_r12
	
	r6min <- rmin**6
	tmp_r6 <- r6[,c(-1,-2,-3)]
	tmp2 <- r6min*tmp_r6
	
	tmp3 <- eps*(tmp1 - 2*tmp2)
		
	tmp <- tmp_info	
	tmp$ener <- rowSums(tmp3)
    tmp <- tmp[order(tmp$ener),]
    
    weighted.mean(ddply(.data=tmp,.var=c("system"),.fun=SLR)$V1,weights)
}

## 12-10-6 potential ##
# Z-score
fitnessZ_eten <- function(pars) {

	list <- combinePARS(pars)
	eps <- list$eps
	eps <- matrix(eps,nrow=m,ncol=n,byrow=T)
	rmin <- list$rmin
	rmin <- matrix(rmin,nrow=m,ncol=n,byrow=T)
	
	r12min <- rmin**12
	tmp_r12 <- r12[,c(-1,-2,-3)]
	tmp1 <- r12min*tmp_r12
	
	r10min <- rmin**10
	tmp_r10 <- r10[,c(-1,-2,-3)]
	tmp2 <- r10min*tmp_r10
	
	r6min <- rmin**6
	tmp_r6 <- r6[,c(-1,-2,-3)]
	tmp3 <- r6min*tmp_r6
	
	tmp4 <- eps*((13*tmp1) - (18*tmp2) + (4*tmp3))
		
	tmp <- tmp_info	
	tmp$ener <- rowSums(tmp4)
    
    -1.0*weighted.mean(ddply(.data=tmp,.var=c("system"),.fun=Zscore)$V1,weights)
}

# SLR
fitnessSLR_eten <- function(pars) {

	list <- combinePARS(pars)
	eps <- list$eps
	eps <- matrix(eps,nrow=m,ncol=n,byrow=T)
	rmin <- list$rmin
	rmin <- matrix(rmin,nrow=m,ncol=n,byrow=T)
	
	r12min <- rmin**12
	tmp_r12 <- r12[,c(-1,-2,-3)]
	tmp1 <- r12min*tmp_r12
	
	r10min <- rmin**10
	tmp_r10 <- r10[,c(-1,-2,-3)]
	tmp2 <- r10min*tmp_r10
	
	r6min <- rmin**6
	tmp_r6 <- r6[,c(-1,-2,-3)]
	tmp3 <- r6min*tmp_r6
	
	tmp4 <- eps*((13*tmp1) - (18*tmp2) + (4*tmp3))
		
	tmp <- tmp_info	
	tmp$ener <- rowSums(tmp4)
	tmp <- tmp[order(tmp$ener),]
    
    weighted.mean(ddply(.data=tmp,.var=c("system"),.fun=SLR)$V1,weights)
}


#######################################################################
# Routines for computing the mean value of the fitness function
#######################################################################
## 12-6 potential ##
calZ_lj <- function(pars){

  	list <- combinePARS(pars)
	eps <- list$eps
	eps <- matrix(eps,nrow=m,ncol=n,byrow=T)
	rmin <- list$rmin
	rmin <- matrix(rmin,nrow=m,ncol=n,byrow=T)
	
	r12min <- rmin**12
	tmp_r12 <- r12[,c(-1,-2,-3)]
	tmp1 <- r12min*tmp_r12
	
	r6min <- rmin**6
	tmp_r6 <- r6[,c(-1,-2,-3)]
	tmp2 <- r6min*tmp_r6
	
	tmp3 <- eps*(tmp1 - 2*tmp2)
	
	tmp <- tmp_info	
	tmp$ener <- rowSums(tmp3)
  	
  	weighted.mean(ddply(.data=tmp,.var=c("system"),.fun=Zscore)$V1,weights)
  	  
}

calSLR_lj <- function(pars){

  	list <- combinePARS(pars)
	eps <- list$eps
	eps <- matrix(eps,nrow=m,ncol=n,byrow=T)
	rmin <- list$rmin
	rmin <- matrix(rmin,nrow=m,ncol=n,byrow=T)
	
	r12min <- rmin**12
	tmp_r12 <- r12[,c(-1,-2,-3)]
	tmp1 <- r12min*tmp_r12
	
	r6min <- rmin**6
	tmp_r6 <- r6[,c(-1,-2,-3)]
	tmp2 <- r6min*tmp_r6
	
	tmp3 <- eps*(tmp1 - 2*tmp2)
	
	tmp <- tmp_info	
	tmp$ener <- rowSums(tmp3)
    tmp <- tmp[order(tmp$ener),]	
  	
  	weighted.mean(ddply(.data=tmp,.var=c("system"),.fun=SLR)$V1,weights)
  	  
}

## 12-10-6 potential ##
calZ_eten <- function(pars){

  	list <- combinePARS(pars)
	eps <- list$eps
	eps <- matrix(eps,nrow=m,ncol=n,byrow=T)
	rmin <- list$rmin
	rmin <- matrix(rmin,nrow=m,ncol=n,byrow=T)
	
	r12min <- rmin**12
	tmp_r12 <- r12[,c(-1,-2,-3)]
	tmp1 <- r12min*tmp_r12
	
	r10min <- rmin**10
	tmp_r10 <- r10[,c(-1,-2,-3)]
	tmp2 <- r10min*tmp_r10
	
	r6min <- rmin**6
	tmp_r6 <- r6[,c(-1,-2,-3)]
	tmp3 <- r6min*tmp_r6
	
	tmp4 <- eps*((13*tmp1) - (18*tmp2) + (4*tmp3))
	
	tmp <- tmp_info	
	tmp$ener <- rowSums(tmp4)
  	
  	weighted.mean(ddply(.data=tmp,.var=c("system"),.fun=Zscore)$V1,weights)
  	  
}

calSLR_eten <- function(pars){

  	list <- combinePARS(pars)
	eps <- list$eps
	eps <- matrix(eps,nrow=m,ncol=n,byrow=T)
	rmin <- list$rmin
	rmin <- matrix(rmin,nrow=m,ncol=n,byrow=T)
	
	r12min <- rmin**12
	tmp_r12 <- r12[,c(-1,-2,-3)]
	tmp1 <- r12min*tmp_r12
	
	r10min <- rmin**10
	tmp_r10 <- r10[,c(-1,-2,-3)]
	tmp2 <- r10min*tmp_r10
	
	r6min <- rmin**6
	tmp_r6 <- r6[,c(-1,-2,-3)]
	tmp3 <- r6min*tmp_r6
	
	tmp4 <- eps*((13*tmp1) - (18*tmp2) + (4*tmp3))
	
	tmp <- tmp_info	
	tmp$ener <- rowSums(tmp4)
	tmp <- tmp[order(tmp$ener),]
  	
  	weighted.mean(ddply(.data=tmp,.var=c("system"),.fun=SLR)$V1,weights)
  	  
}


#######################################################################
# Routines for computing the fitness function for each individual system
# included in the opitimization
#######################################################################
## 12-6 potential ##
calZ_lj_indiv <- function(pars){

  	list <- combinePARS(pars)
	eps <- list$eps
	eps <- matrix(eps,nrow=m,ncol=n,byrow=T)
	rmin <- list$rmin
	rmin <- matrix(rmin,nrow=m,ncol=n,byrow=T)
	
	r12min <- rmin**12
	tmp_r12 <- r12[,c(-1,-2,-3)]
	tmp1 <- r12min*tmp_r12
	
	r6min <- rmin**6
	tmp_r6 <- r6[,c(-1,-2,-3)]
	tmp2 <- r6min*tmp_r6
	
	tmp3 <- eps*(tmp1 - 2*tmp2)
	
	tmp <- tmp_info	
	tmp$ener <- rowSums(tmp3)

  	ddply(.data=tmp,.var=c("system"),.fun=Zscore)$V1
  	  
}

calSLR_lj_indiv <- function(pars){

  	list <- combinePARS(pars)
	eps <- list$eps
	eps <- matrix(eps,nrow=m,ncol=n,byrow=T)
	rmin <- list$rmin
	rmin <- matrix(rmin,nrow=m,ncol=n,byrow=T)
	
	r12min <- rmin**12
	tmp_r12 <- r12[,c(-1,-2,-3)]
	tmp1 <- r12min*tmp_r12
	
	r6min <- rmin**6
	tmp_r6 <- r6[,c(-1,-2,-3)]
	tmp2 <- r6min*tmp_r6
	
	tmp3 <- eps*(tmp1 - 2*tmp2)
	
	tmp <- tmp_info	
	tmp$ener <- rowSums(tmp3)
    tmp <- tmp[order(tmp$ener),]	
  	
  	ddply(.data=tmp,.var=c("system"),.fun=SLR)$V1
  	  
}

## 12-10-6 potential ##
calZ_eten_indiv <- function(pars){

  	list <- combinePARS(pars)
	eps <- list$eps
	eps <- matrix(eps,nrow=m,ncol=n,byrow=T)
	rmin <- list$rmin
	rmin <- matrix(rmin,nrow=m,ncol=n,byrow=T)
	
	r12min <- rmin**12
	tmp_r12 <- r12[,c(-1,-2,-3)]
	tmp1 <- r12min*tmp_r12
	
	r10min <- rmin**10
	tmp_r10 <- r10[,c(-1,-2,-3)]
	tmp2 <- r10min*tmp_r10
	
	r6min <- rmin**6
	tmp_r6 <- r6[,c(-1,-2,-3)]
	tmp3 <- r6min*tmp_r6
	
	tmp4 <- eps*((13*tmp1) - (18*tmp2) + (4*tmp3))
	
	tmp <- tmp_info	
	tmp$ener <- rowSums(tmp4)
  	
  	ddply(.data=tmp,.var=c("system"),.fun=Zscore)$V1
  	  
}

calSLR_eten_indiv <- function(pars){

  	list <- combinePARS(pars)
	eps <- list$eps
	eps <- matrix(eps,nrow=m,ncol=n,byrow=T)
	rmin <- list$rmin
	rmin <- matrix(rmin,nrow=m,ncol=n,byrow=T)
	
	r12min <- rmin**12
	tmp_r12 <- r12[,c(-1,-2,-3)]
	tmp1 <- r12min*tmp_r12
	
	r10min <- rmin**10
	tmp_r10 <- r10[,c(-1,-2,-3)]
	tmp2 <- r10min*tmp_r10
	
	r6min <- rmin**6
	tmp_r6 <- r6[,c(-1,-2,-3)]
	tmp3 <- r6min*tmp_r6
	
	tmp4 <- eps*((13*tmp1) - (18*tmp2) + (4*tmp3))
	
	tmp <- tmp_info	
	tmp$ener <- rowSums(tmp4)
	tmp <- tmp[order(tmp$ener),]
  	
  	ddply(.data=tmp,.var=c("system"),.fun=SLR)$V1
  	  
}


#######################################################################
# Routines to evaluate energies after optimization
#######################################################################
## 12-6 potential ##
calE_lj <- function(pars,r12,r6,m,n,tmp_info) {

	list <- combinePARS(pars)
	eps <- list$eps
	eps <- matrix(eps,nrow=m,ncol=n,byrow=T)
	rmin <- list$rmin
	rmin <- matrix(rmin,nrow=m,ncol=n,byrow=T)
	
	r12min <- rmin**12
	tmp_r12 <- r12[,c(-1,-2,-3)]
	tmp1 <- r12min*tmp_r12
	
	r6min <- rmin**6
	tmp_r6 <- r6[,c(-1,-2,-3)]
	tmp2 <- r6min*tmp_r6
	
	tmp3 <- eps*(tmp1 - 2*tmp2)
	
	tmp <- data.frame(flag=tmp_info[,1],system=tmp_info[,2],rmsd=tmp_info[,3],ener=rowSums(tmp3))
	
}

## 12-10-6 potential ##
calE_eten <- function(pars,r12,r10,r6,m,n,tmp_info) {

	list <- combinePARS(pars)
	eps <- list$eps
	eps <- matrix(eps,nrow=m,ncol=n,byrow=T)
	rmin <- list$rmin
	rmin <- matrix(rmin,nrow=m,ncol=n,byrow=T)
	
	r12min <- rmin**12
	tmp_r12 <- r12[,c(-1,-2,-3)]
	tmp1 <- r12min*tmp_r12
	
	r10min <- rmin**10
	tmp_r10 <- r10[,c(-1,-2,-3)]
	tmp2 <- r10min*tmp_r10
	
	r6min <- rmin**6
	tmp_r6 <- r6[,c(-1,-2,-3)]
	tmp3 <- r6min*tmp_r6
	
	tmp4 <- eps*((13*tmp1) - (18*tmp2) + (4*tmp3))
	
	tmp <- data.frame(flag=tmp_info[,1],system=tmp_info[,2],rmsd=tmp_info[,3],ener=rowSums(tmp4))

}


#######################################################################
# Routine to compute enrichment score
#######################################################################
enrichment <- function(tt,thresholds=seq(0.05,0.20,0.05)){
    enrich <- NULL
    for (threshold in thresholds) {
		nhits <- floor(threshold*nrow(tt))
		tt$index_energy <- rownames(tt)[order(tt$ener)]
		tt$index_rmsd <- rownames(tt)[order(tt$rmsd)]
		enrich <- c(enrich,sum(tt$index_energy[1:nhits] %in% tt$index_rmsd[1:nhits])/((threshold^2)*nrow(tt)))
	}
	enrich <- as.data.frame(t(enrich))
	colnames(enrich) <- paste("ES_",thresholds,sep="")
	enrich$pdb <- as.character(pdbs[unique(tt$system)])
	return(enrich)
}