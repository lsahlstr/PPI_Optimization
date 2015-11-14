#######################################################################
# R code containing a subroutine to construct the 210 inteaction energy 
# and 210 interaction radii as part of a Genetic Algorithm-based force-field 
# optimization for coarse-grained protein-protein binding simulations.
# The parameters are generated by considering 2 x 210 parameters (i.e., 
# each amino acid pair is considered separately) or 2 x 20 parameters 
# (i.e., all possible pairwise energies are constructed based upon the 
# geometric mean).
#
# 	The initial version of the routine was written to optimize energy parameters 
#	with different fitness functions
#		Aaron Frank, U. Michigan, c. 2015
#	The subroutine was further developed to consider energies and distances 
#	simultaneously
#		Logan S. Ahlstrom, U. Michigan c. 2015
#######################################################################


#######################################################################
# Get the 210 interaciton energies (eps_ij) from 20 eps_i values
# 	eps_ij = (eps_i * eps_j)**0.5
#	rmin_ij = (rmin_i + rmin_j)/2
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
					matEPS <- c(matEPS,sqrt(eps_mask[i]*eps_mask[j]))
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