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
    
    -1.0*mean(ddply(.data=tmp,.var=c("system"),.fun=Zscore)$V1)
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
    
    mean(ddply(.data=tmp,.var=c("system"),.fun=SLR)$V1)
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
    
    -1.0*mean(ddply(.data=tmp,.var=c("system"),.fun=Zscore)$V1)
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
    
    mean(ddply(.data=tmp,.var=c("system"),.fun=SLR)$V1)
}


#######################################################################
# Routines for computing the fitness function at each iteration after optimization
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
  	
  	mean(ddply(.data=tmp,.var=c("system"),.fun=Zscore)$V1)
  	  
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
  	
  	mean(ddply(.data=tmp,.var=c("system"),.fun=SLR)$V1)
  	  
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
  	
  	mean(ddply(.data=tmp,.var=c("system"),.fun=Zscore)$V1)
  	  
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
  	
  	mean(ddply(.data=tmp,.var=c("system"),.fun=SLR)$V1)
  	  
}


#######################################################################
# Routines to evaluate energies after optimization
#######################################################################
## 12-6 potential ##
calE_lj <- function(pars) {

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
calE_eten <- function(pars) {

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