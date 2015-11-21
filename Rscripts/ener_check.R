#######################################################################
# Check energies computed in R with CHARMM energies
#######################################################################
ener_check <- function(tol=0.01){
		
	ener_ref_file <- read.table("charmm.ener")
	ener_ref <- ener_ref_file$V1
	
	if (potFlag == "lj") {
		check <- ener_lj(ipars)
	} else if (potFlag == "eten"){
		check <- ener_eten(ipars)
	} else if (potFlag == "etsr"){
		check <- ener_etsr(ipars)
	}
	ener_comp <- check$ener
	
	diff <- abs(ener_comp - ener_ref)
	
	check_tol <- sum(diff > tol)
		
	if (check_tol == 0) {
		cat("Passed energy test with tol =", tol, ".\n\n")
	} else {
		cat("Failed energy check with tol =", tol, "\nExiting.\n\n")
		cat(sprintf("%s\n\n",date()))
		#stop()
	}
	
	plot(ener_ref,ener_comp)
	
}