#######################################################################
# Inter-protein potential functions
#######################################################################
## 12-6 potential ##
ener_lj <- function(pars){
	
	lsize <- 2
	#lenrow <- 210*lsize
	
  	list <- combinePARS(pars)
  	  	
	eps <- list$eps
	rmin <- list$rmin
	
	big_eps <- matrix(nrow=210,ncol=lsize)
	big_rmin <- matrix(nrow=210,ncol=lsize)
	for (i in 1:210) {
		tmp_eps <- rep(eps[i],lsize)
		tmp_rmin <- rep(rmin[i],lsize)
		big_eps[i,] <- tmp_eps
		big_rmin[i,] <- tmp_rmin
	}	
	big_eps <- as.vector(t(big_eps))
	big_rmin <- as.vector(t(big_rmin))
	
	list <- list(big_eps=big_eps,big_rmin=big_rmin)
	
	return(list)
	
	#eps <- matrix(eps,nrow=m,ncol=n,byrow=T)
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