#######################################################################
# Fitness functions
#######################################################################
# Z-score
Zscore <- function(tmp) {
    ((mean(tmp$ener[tmp$flag==1])-mean(tmp$ener[tmp$flag==0]))/sd(tmp$ener[tmp$flag==0]))
}

# Sum of Logarithms of Ranks (SLR) early recognition metric
SLR <- function(tmp){
    ri <- which(tmp$flag==1)
    N <- nrow(tmp)
    i <-  1:length(ri)
    SLRmax <- -sum(log(i/N))
    return(-sum(log(ri/N))/SLRmax)
}

# Z-score with 12-6 potential
fitnessZ_lj <- function(pars){
	tmp <- ener_lj(pars)
  	-1.0*weighted.mean(ddply(.data=tmp,.var=c("system"),.fun=Zscore)$V1,weights)
}

# SLR with 12-6 potential
fitnessSLR_lj <- function(pars){
	tmp <- ener_lj(pars)
	tmp <- tmp[order(tmp$ener),]  	
  	weighted.mean(ddply(.data=tmp,.var=c("system"),.fun=SLR)$V1,weights)  	
}

# Z-score with 12-10-6 potential ("ETEN")
fitnessZ_eten <- function(pars){	
	tmp <- ener_eten(pars)
  	-1.0*weighted.mean(ddply(.data=tmp,.var=c("system"),.fun=Zscore)$V1,weights) 	
}

# SLR with 12-10-6 potential ("ETEN")
fitnessSLR_eten <- function(pars){	
	tmp <- ener_eten(pars)
	tmp <- tmp[order(tmp$ener),]  	
  	weighted.mean(ddply(.data=tmp,.var=c("system"),.fun=SLR)$V1,weights) 	  	
}

# Z-score with short-range 12-10-6 potential ("ETSR")
fitnessZ_etsr <- function(pars){
	tmp <- ener_etsr(pars)
  	-1.0*weighted.mean(ddply(.data=tmp,.var=c("system"),.fun=Zscore)$V1,weights)  	
}

# SLR with short-range 12-10-6 potential ("ETSR")
fitnessSLR_etsr <- function(pars){
	tmp <- ener_etsr(pars)  	
	tmp <- tmp[order(tmp$ener),]
  	weighted.mean(ddply(.data=tmp,.var=c("system"),.fun=SLR)$V1,weights)
}