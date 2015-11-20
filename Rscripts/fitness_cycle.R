#######################################################################
# Get the mean Z-score and SLR at each cycle
#######################################################################
fitness_cycle <- function() {

	zscore <- NULL
	slr <- NULL

	if (potFlag == "lj") {
		for (i in 1:iters) {
			zscore <- c(zscore,-1*fitnessZ_lj(as.vector(GAReal@bestSol[[i]][1,])))
			slr <- c(slr,fitnessSLR_lj(as.vector(GAReal@bestSol[[i]][1,])))
		}
	} else if (potFlag == "eten") {
		for (i in 1:iters) {
			zscore <- c(zscore,-1*fitnessZ_eten(as.vector(GAReal@bestSol[[i]][1,])))
			slr <- c(slr,fitnessSLR_eten(as.vector(GAReal@bestSol[[i]][1,])))
		}
	} else if (potFlag == "etsr") {
		for (i in 1:iters) {
			zscore <- c(zscore,-1*fitnessZ_etsr(as.vector(GAReal@bestSol[[i]][1,])))
			slr <- c(slr,fitnessSLR_etsr(as.vector(GAReal@bestSol[[i]][1,])))
		}
	}

	zscore <- data.frame(iters=1:iters,score=zscore)
	slr <- data.frame(iters=1:iters,score=slr)
	
	write.table(zscore,file="zscore.dat",row.names=F,quote=F,col.names=F)
	write.table(slr,file="slr.dat",row.names=F,quote=F,col.names=F)

}