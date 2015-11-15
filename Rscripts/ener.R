#######################################################################
# Inter-protein potential functions
#######################################################################
## 12-6 potential ##
ener_lj <- function(pars){
	
	#lsize <- 1000 (Defined globally in main.R now)
	
  	list <- combinePARS(pars)
  	  	
	eps <- list$eps
	rmin <- list$rmin
	
	big_eps_mat <- matrix(nrow=nrow(rij_rmsd_data),ncol=(ncol(rij_rmsd_data) - 4))
	big_rmin_mat <- matrix(nrow=nrow(rij_rmsd_data),ncol=(ncol(rij_rmsd_data) - 4))
	
	for (j in 1:nrow(rij_rmsd_data)) {
			
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
		
		#big_eps_mat <- rbind(big_eps_mat,big_eps)
		big_eps_mat[j,] <- c(big_eps)
		big_rmin_mat[j,] <- c(big_rmin)
		
	}
	
	rij <- rij_rmsd_data[,c(-1,-2,-3,-4)]
	
	rmin_12 <- big_rmin_mat**12
	rij_12 <- rij**12
	tmp1 <- rmin_12/rij_12
	
	rmin_6 <- big_rmin_mat**6
	rij_6 <- rij**6
	tmp2 <- rmin_6/rij_6
	
	tmp3 <- big_eps_mat*(tmp1 - 2*tmp2)
	
	tmp <- sysinfo	
	tmp$ener <- rowSums(tmp3)
  	
  	return(tmp)
  	
}