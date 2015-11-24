#######################################################################
# Combine interaction energy (eps) and radii (rmin) parameters
#######################################################################
combinePARS <- function(pars) {
	
	matEPS <- NULL
	matRMIN <- NULL
	
	len <- length(ipars)
	eps <- pars[1:(len/2)]
	rmin <- pars[((len/2)+1):len]
	
	eps_mask <- mask*eps
	
	if ((len/2) == 20) {
		for (i in 1:(len/2)){
			for (j in 1:(len/2)){
				if (i <= j){
					# eps_ij = (eps_i * eps_j)**0.5
					matEPS <- c(matEPS,sqrt(eps_mask[i]*eps_mask[j]))
					# rmin_ij = (rmin_i + rmin_j)/2
					matRMIN <- c(matRMIN,(rmin[i]+rmin[j])/2)
				}
			}
		}
		
    } else if ((len/2) == 210) {
    	for (i in 1:(len/2)){
    		#j <- i+length(eps)
			matEPS <- c(matEPS,eps_mask[i])
			matRMIN <- c(matRMIN,rmin[i])
		}
    
    } else {
    	stop("wrong number of parameters, i.e., not 20 or 210 eps_ij (rmin_ij)")
    }
    
    list <- list(eps=matEPS,rmin=matRMIN)
    return(list)
}