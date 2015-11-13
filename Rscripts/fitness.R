#######################################################################
# Fitness functions
#######################################################################
# Z-score
Zscore <- function(tmp) {
    ((mean(tmp$ener[tmp$flag==1])-mean(tmp$ener[tmp$flag==0]))/sd(tmp$ener[tmp$flag==0]))
    #((mean(tmp$ener[tmp$flag==1])-mean(tmp$ener[tmp$flag==0]))/sd(tmp$ener[tmp$flag==0])) + 1/sd(tmp$ener[tmp$flag==1])
}

# Sum of Logarithms of Ranks (SLR) early recognition metric
SLR <- function(result){
    ri <- which(result$flag==1)
    N <- nrow(result)
    i <-  1:length(ri)
    SLRmax <- -sum(log(i/N))
    return(-sum(log(ri/N))/SLRmax)
    #return(-sum(log(ri/N))/SLRmax + 1/sd(results$ener[tmp$flag==1]))
}