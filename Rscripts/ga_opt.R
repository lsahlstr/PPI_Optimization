ga_opt <- function() {
	# Define min and max
	min_eps <- NULL
	max_eps <- NULL
	min_rmin <- NULL
	max_rmin <- NULL

	for (i in 1:len) {
		# eps
		if (i <= (len/2)) {	
			if (opttypeFlag == "gentle") {
				min_eps[i] <- (ipars[i] - (0.5*ipars[i]))
				max_eps[i] <- (ipars[i] + (0.5*ipars[i]))
			} else {
				min_eps[i] <- eps_min_val
				max_eps[i] <- eps_max_val
			}
		# rmin		
		} else {
			min_rmin[(i-(len/2))] <- (ipars[i] - (1.0*rminSD[(i-(len/2))]))
			max_rmin[(i-(len/2))] <- (ipars[i] + (1.0*rminSD[(i-(len/2))]))
		}
	}
	min <- c(min_eps,min_rmin) 
	max <- c(max_eps,max_rmin)

	# GA to assign weights; fitness function = Z-score for single system
	ffunc <- NULL
	if (potFlag == "lj") {
		if (fitFlag == "zscore") {
			ffunc <- fitnessZ_lj
		} else {
			ffunc <- fitnessSLR_lj
		}
	} else if (potFlag == "eten"){
		if (fitFlag == "zscore") {
			ffunc <- fitnessZ_eten
		} else {
			ffunc <- fitnessSLR_eten
		}
	} else if (potFlag == "etsr"){
		if (fitFlag == "zscore") {
			ffunc <- fitnessZ_etsr
		} else {
			ffunc <- fitnessSLR_etsr
		}
	}

	# Run GA
	cat(sprintf("%s\n\n",date()))
	GAReal <- ga(type = "real-valued", fitness=ffunc, min=min, max=max, popSize=popSize, maxiter=iters, suggestions=jitter(initialSolution), keepBest=T, parallel=TRUE) 
	cat(sprintf("%s\n\n",date()))

	# Save data from all GA cycles to R data structure
	save(GAReal,file="GA.RData")

	#list <- list(GAReal=GAReal,bestPars=bestPars)
	return(GAReal)

}