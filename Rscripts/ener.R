#######################################################################
# Inter-protein potential functions
#######################################################################
## 12-6 potential ##
ener_lj <- function(pars){
	
	# Make large eps and rmin data structures; rij defined globally in main.R	
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