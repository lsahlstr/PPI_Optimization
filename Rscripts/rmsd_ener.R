#######################################################################
# Compute energies with initial ("old") and optimized ("new") parameters
# and combine with RMSD and system information. Organize energy and RMSD 
# information and old and new parameters, and then write to output files.
#######################################################################
rmsd_ener <- function() {

	if (potFlag == "lj") {
		rmsd_ener <- ener_lj(bestPars)
		rmsd_ener$oldenergy <- ener_lj(ipars)$ener
	} else if (potFlag == "eten"){
		rmsd_ener <- ener_eten(bestPars)
		rmsd_ener$oldenergy <- ener_eten(ipars)$ener
	} else if (potFlag == "etsr"){
		rmsd_ener <- ener_etsr(bestPars)
		rmsd_ener$oldenergy <- ener_etsr(ipars)$ener
	}

	rmsd_ener <- rmsd_ener[order(rmsd_ener$system,rmsd_ener$ener),]

	# Old and new parameters
	old_new <- data.frame(respair=respairs, 
		eij_initial=ipars[1:(len/2)], 
		eij_optimized=bestPars[1:(len/2)],
		rij_initial=ipars[((len/2)+1):len],
		rij_optimized=bestPars[((len/2)+1):len])
	
	write.table(old_new,file="OldNewComp.txt",row.names=F,quote=F,col.names=T)
	write.table(rmsd_ener,file="rmsdEnerComp.txt",row.names=F,quote=F,col.names=F)
	
	return(rmsd_ener)
	
}