#######################################################################
# Inter-protein potential functions
#######################################################################
## 12-6 potential ##
ener_lj <- function(pars) {
	
	# epsilon
  	big_eps <- big_eps(pars)
		
	# r^12 term
	rmin_12 <- big_rmin(pars,12)
	tmp1 <- rmin_12*rij_12
	
	# r^6 term
	rmin_6 <- big_rmin(pars,6)
	tmp2 <- rmin_6*rij_6
	
	# Energy per i,j interaction
	tmp3 <- big_eps*(tmp1 - 2*tmp2)
	
	# Get back system info
	tmp <- sysinfo
	
	# Append energy per conformer to tmp data structure
	tmp$ener <- rowSums(tmp3)
  	
  	return(tmp)
  	
}

## 12-10-6 potential ("ETEN") ##
ener_eten <- function(pars) {
	
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