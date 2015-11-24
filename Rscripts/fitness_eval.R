#######################################################################
# Routines for evaluating fitness functions
#######################################################################
## Z-score with 12-6 potential ##
fitnessZ_lj <- function(pars){
	tmp <- ener_lj(pars)
  	-1.0*weighted.mean(ddply(.data=tmp,.var=c("system"),.fun=Zscore)$V1,weights)
}

## SLR with 12-6 potential ##
fitnessSLR_lj <- function(pars){
	tmp <- ener_lj(pars)
	tmp <- tmp[order(tmp$ener),]  	
  	weighted.mean(ddply(.data=tmp,.var=c("system"),.fun=SLR)$V1,weights)  	
}

## Z-score with 12-10-6 potential ("ETEN") ##
fitnessZ_eten <- function(pars){	
	tmp <- ener_eten(pars)
  	-1.0*weighted.mean(ddply(.data=tmp,.var=c("system"),.fun=Zscore)$V1,weights) 	
}

## SLR with 12-10-6 potential ("ETEN") ##
fitnessSLR_eten <- function(pars){	
	tmp <- ener_eten(pars)
	tmp <- tmp[order(tmp$ener),]  	
  	weighted.mean(ddply(.data=tmp,.var=c("system"),.fun=SLR)$V1,weights) 	  	
}

## Z-score with short-range 12-10-6 potential ("ETSR")##
fitnessZ_etsr <- function(pars){
	tmp <- ener_etsr(pars)
  	-1.0*weighted.mean(ddply(.data=tmp,.var=c("system"),.fun=Zscore)$V1,weights)  	
}

## SLR with short-range 12-10-6 potential ("ETSR") ##
fitnessSLR_etsr <- function(pars){
	tmp <- ener_etsr(pars)  	
	tmp <- tmp[order(tmp$ener),]
  	weighted.mean(ddply(.data=tmp,.var=c("system"),.fun=SLR)$V1,weights)
}