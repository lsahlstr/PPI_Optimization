#######################################################################
# Make big eps and rmin data structures, i.e. the same size as rij
#######################################################################
big_pars <- function(pars){

	list <- combinePARS(pars)
  	  	
	eps <- unlist(lapply(list$eps,rep,lsize))
	rmin <- unlist(lapply(list$rmin,rep,lsize))

	big_eps_mat <- matrix(eps,nrow=nrow(rij),ncol=ncol(rij),byrow=TRUE)
	big_rmin_mat <- matrix(rmin,nrow=nrow(rij),ncol=ncol(rij),byrow=TRUE)
	
	list <- list(big_eps_mat=big_eps_mat,big_rmin_mat=big_rmin_mat)
    return(list)
    	
}