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