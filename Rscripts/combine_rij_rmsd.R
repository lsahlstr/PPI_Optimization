#######################################################################
# Combine distance information
#######################################################################
combineRijRMSD <- function(pdbs,lsize=1000) {		

	lengthm <- 4 + (210 * lsize)
	superM2 <- matrix(ncol=lengthm)
	
	for (p in seq_along(pdbs)){
		pdb <- pdbs[p]
		superM2 <- rbind(superM2,formatRijRMSD("n",pdb,p,lengthm,lsize))
		superM2 <- rbind(superM2,formatRijRMSD("nn",pdb,p,lengthm,lsize))
	}
	
	return(superM2[-1,])

}