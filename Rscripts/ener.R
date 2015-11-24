#######################################################################
# Inter-protein potential functions
#######################################################################
## 12-6 potential ##
ener_lj <- function(pars) {
	
	# Make large eps and rmin data structures	
  	list <- big_pars (pars)
  	big_eps_mat <- list$big_eps_mat
  	big_rmin_mat <- list$big_rmin_mat
		
	# r^12 term
	rmin_12 <- big_rmin_mat**12
	rij_12 <- rij**12
	tmp1 <- rmin_12/rij_12
	tmp1[is.infinite(tmp1)] <- 0
	
	# r^6 term
	rmin_6 <- big_rmin_mat**6
	rij_6 <- rij**6
	tmp2 <- rmin_6/rij_6
	tmp2[is.infinite(tmp2)] <- 0
	
	# Energy per i,j interaction
	tmp3 <- big_eps_mat*(tmp1 - 2*tmp2)
	
	# Get back system info
	tmp <- sysinfo
	
	# Append energy per conformer to tmp data structure
	tmp$ener <- rowSums(tmp3)
  	
  	return(tmp)
  	
}

## 12-10-6 potential ("ETEN") ##
ener_eten <- function(pars) {
	
	# Make large eps and rmin data structures
  	list <- big_pars (pars)
  	big_eps_mat <- list$big_eps_mat
  	big_rmin_mat <- list$big_rmin_mat
		
	# r^12 term
	rmin_12 <- big_rmin_mat**12
	rij_12 <- rij**12
	tmp1 <- rmin_12/rij_12
	tmp1[is.infinite(tmp1)] <- 0

	# r^10 term
    rmin_10 <- big_rmin_mat**10
    rij_10 <- rij**10
    tmp2 <- rmin_10/rij_10
    tmp2[is.infinite(tmp2)] <- 0
	
	# r^6 term
	rmin_6 <- big_rmin_mat**6
	rij_6 <- rij**6
	tmp3 <- rmin_6/rij_6
	tmp3[is.infinite(tmp3)] <- 0
	
	# Energy per i,j interaction
	tmp4 <- big_eps_mat*(13*tmp1 - 18*tmp2 + 4*tmp3)
	
	# Get back system info
	tmp <- sysinfo
	
	# Append energy per conformer to tmp data structure
	tmp$ener <- rowSums(tmp4)
  	
  	return(tmp)
  	
}

## Short-range 12-10-6 potential ("ETSR") ##
ener_etsr <- function(pars) {
	
	# epsilon
  	big_eps <- big_eps(pars)
  		
	# r^12 term
	rmin_12 <- big_rmin(pars,12)
	tmp1 <- rmin_12*rij_12

	# r^10 term
    rmin_10 <- big_rmin(pars,10)
	tmp2 <- rmin_10*rij_10
	
	# r^6 term
	rmin_6 <- big_rmin(pars,6)
	tmp3 <- rmin_6*rij_6
	
	# Energy per i,j interaction
	tmp4 <- big_eps*(13*tmp1 - 18*tmp2 + 4*tmp3)
	
	# f(r)
	tmp4 <- 1/(1 + (1/(rij_12*(129.7463*rmin_12))))

	# Energy per i,j interaction
	tmp5 <- big_eps*((13*(tmp1*tmp4) - 18*(tmp2*tmp4) + 4*(tmp3*tmp4)))
	
	# Get back system info
	tmp <- sysinfo
	
	# Append energy per conformer to tmp data structure
	tmp$ener <- rowSums(tmp5)
  	
  	return(tmp)
  	
}