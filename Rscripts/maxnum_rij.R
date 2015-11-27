#######################################################################
# Find maximum number of i,j interactions for a single system
#######################################################################
maxnum_rij <- function(pool="n",pdb) {
		
	poolRMSD <- read.table(paste(pdb,"/",pool,"_rmsd.dat",sep=""))$V2

	maxnum <- 0
	for (conf in 1:length(poolRMSD)) {

		# read in data
		R <- read.table(paste(pdb,"/",pool,"_rij.",conf,".txt",sep=""),col.names=c("pair","rij"))

		# order based on pair (resnames of the ij pair)
		R <- R[order(R$pair),]
		R$pair <- as.character(R$pair)
		
		# Reference with all 210 aimino acid pairs
		names <- as.character(respairs)

		# Loop over ij interactions
		for (i in seq_along(names)) {
			name <- names[i]
			tmp <- c(subset(R,pair==name)$rij)
			m <- length(tmp)
			cat(sprintf("%s %s%d: %s => %d\n",pdb,pool,conf,name,m))
			if (m > maxnum) {
				maxnum <- m
			}		
		}
	}
	return(maxnum)	
}