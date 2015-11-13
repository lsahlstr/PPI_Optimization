#######################################################################
# Inter-protein potential functions
#######################################################################
## 12-6 potential ##
ener_lj <- function(pars){

	dim(rij_rmsd_data)
  	list <- combinePARS(pars)
  	
	#eps <- list$eps
	#eps <- matrix(eps,nrow=m,ncol=n,byrow=T)
	#rmin <- list$rmin
	#rmin <- matrix(rmin,nrow=m,ncol=n,byrow=T)
	
	
	
	#r12min <- rmin**12
	#tmp_r12 <- r12[,c(-1,-2,-3)]
	#tmp1 <- r12min*tmp_r12
	
	#r6min <- rmin**6
	#tmp_r6 <- r6[,c(-1,-2,-3)]
	#tmp2 <- r6min*tmp_r6
	
	#tmp3 <- eps*(tmp1 - 2*tmp2)
	
	#tmp <- tmp_info	
	#tmp$ener <- rowSums(tmp3)
  	
  	#weighted.mean(ddply(.data=tmp,.var=c("system"),.fun=Zscore)$V1,weights)
  	  
}