#######################################################################
# Make big eps and rmin data structures, i.e. the same size as rij
#######################################################################
big_rmin <- function(pars,n=12){

	list <- combinePARS(pars)
	
	rmin <- unlist(lapply((list$rmin**n),rep,lsize))
	big_rmin_mat <- matrix(rmin,nrow=nrow(rij),ncol=ncol(rij),byrow=TRUE)
	
    return(big_rmin_mat)
    	
}

big_eps <- function(pars,n=12){

	list <- combinePARS(pars)
	
	eps <- unlist(lapply((list$eps),rep,lsize))
	big_eps_mat <- matrix(eps,nrow=nrow(rij),ncol=ncol(rij),byrow=TRUE)
	
    return(big_eps_mat)
    	
}