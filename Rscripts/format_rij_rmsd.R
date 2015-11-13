#######################################################################
# Format distance information
#######################################################################
formatRijRMSD <- function(pool="n",pdb,p,lengthm,lsize) {
	
	superM <- matrix(ncol=lengthm)
		
	poolRMSD <- read.table(paste(pdb,"/",pool,"_rmsd.dat",sep=""))$V2

	for (conf in 1:length(poolRMSD)) {

		# read in data
		R <- read.table(paste(pdb,"/",pool,"_rij.",conf,".txt",sep=""),col.names=c("pair","rij"))

		# Get system info		
		info <- NULL
		info[1] <- 0
		if (pool == "n") {
			info[1] <- 1
		}
		info[2] <- p
		info[3] <- poolRMSD[conf]
		info[4] <- conf

		# order based on pair (resnames of the ij pair)
		R <- R[order(R$pair),]
		R$pair <- as.character(R$pair)
		
		# Reference with all 210 aimino acid pairs
		names <- as.character(respairs)

		# Loop over ij keys and create list (and matrix)
		# initial list l and matrix m
		l <- list(NULL)
		m <- matrix(nrow=length(names),ncol=lsize)

		for (i in seq_along(names)) {
			name <- names[i]
			tmp <- c(subset(R,pair==name)$rij)
			nzeros <- lsize-length(tmp)

			if(!is.null(tmp)){
				tmp <- as.vector(padzeros(tmp, nzeros, side="right"))
			} else {
				tmp <- rep(0,lsize)
			}						
			l[[name]] <- tmp
			m[i,] <- tmp					
		}
		m <- c(info,as.vector(t(m)))
		superM <- rbind(superM,m)
	}
		
	return(superM[-1,])
	
}