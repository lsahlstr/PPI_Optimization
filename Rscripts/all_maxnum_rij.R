#######################################################################
# Get maximum number of i,j interactions over all systems
#######################################################################
all_maxnum_rij <- function(pdbs) {		

	#cat(sprintf("all_maxnum_rij\n"))
	
	lsize <- 0
	for (p in seq_along(pdbs)){
		pdb <- pdbs[p]
		#cat(sprintf("%s\n",pdb))
		
		max_rij <- maxnum_rij("n",pdb)
		if (max_rij > lsize) {
			lsize <- max_rij
		}
		
		max_rij <- maxnum_rij("nn",pdb)
		if (max_rij > lsize) {
			lsize <- max_rij
		}
		
	}
	
	return(lsize)

}