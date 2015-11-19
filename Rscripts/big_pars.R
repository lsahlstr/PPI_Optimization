#######################################################################
# Make big eps and rij data structures, i.e. the same size as rij_rmsd_data
#######################################################################
big_pars <- function(pars){

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
	
	list <- list(big_eps_mat=big_eps_mat,big_rmin_mat=big_rmin_mat)
    return(list)
    	
}