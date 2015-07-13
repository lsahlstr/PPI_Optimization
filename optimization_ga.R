# load required libraries
library(GA)
library(plyr)

# define fitness functions
SLR <- function(result){
    #The Sum of Logarithms of Ranks (SLR) metric [Early recognition Metric]
    #see: Comprehensive Comparison of Ligand-Based Virtual Screening Tools Against the DUD Data set Reveals Limitations of Current 3D Methods
    ri <- which(result$flag==1)
    N <- nrow(result)
    i <-  1:length(ri)
    SLRmax <- -sum(log(i/N))
    return(-sum(log(ri/N))/SLRmax)
}

Zscore <- function(tmp){
    #((mean(tmp$ener[tmp$flag==1])-mean(tmp$ener[tmp$flag==0]))/sd(tmp$ener[tmp$flag==0])) + 1/sd(tmp$ener[tmp$flag==1])
    ((mean(tmp$ener[tmp$flag==1])-mean(tmp$ener[tmp$flag==0]))/sd(tmp$ener[tmp$flag==0]))
}

combineEPS <- function(eps){
    mat <- NULL
    for (i in 1:210){
        #for (j in 1:210){
            #if (i <= j){
		mat <- c(mat,eps[i])
                #mat <- c(mat,sqrt(eps[i]*eps[j]))
            #}
        #}
    }
    return(mat)
}


fitness <- function(eps){
  tmp <- comb[,c(-1,-2,-3)]
  #eps <- combineEPS(eps)
  eps <- matrix(eps,nrow=nrow(tmp),ncol=ncol(tmp),byrow=T)
  tmp <- data.frame(flag=comb[,1],system=comb[,2],rmsd=comb[,3],ener=rowSums(eps*tmp))
  tmp <- tmp[order(tmp$ener),]
  
  mean(ddply(.data=tmp,.var=c("system"),.fun=SLR)$V1)
  #+(cor(tmp$ener,tmp$rmsd))
}


calZ <- function(eps){
  tmp <- comb[,c(-1,-2,-3)]
  #eps <- combineEPS(eps)
  eps <- matrix(eps,nrow=nrow(tmp),ncol=ncol(tmp),byrow=T)
  tmp <- data.frame(flag=comb[,1],system=comb[,2],rmsd=comb[,3],ener=rowSums(eps*tmp))
  mean(ddply(.data=tmp,.var=c("system"),.fun=Zscore)$V1)
}

fitnessZ <- function(eps){
	eps <- matrix(eps,nrow=m,ncol=n,byrow=T)
	tmp <- tmp_comb_A
	tmp$ener <- rowSums(eps*tmp_comb)
  # the next line is the bottleneck in this function
  # will try to resolve this but for now if do optimization over a single system
  # use fitnessZ2 below
  -1.0*mean(ddply(.data=tmp,.var=c("system"),.fun=Zscore)$V1)
}

fitnessZ2 <- function(eps){
	# only use when optimizing parameters for SINGLE system
	eps <- matrix(eps,nrow=m,ncol=n,byrow=T)
	tmp <- tmp_comb_A
	tmp$ener <- rowSums(eps*tmp_comb)
	-1.0*Zscore(tmp)
}

calE <- function(eps){
  tmp <- comb[,c(-1,-2,-3)]
  eps <- combineEPS(eps)
  eps <- matrix(eps,nrow=nrow(tmp),ncol=ncol(tmp),byrow=T)
  tmp <- data.frame(flag=comb[,1],system=comb[,2],rmsd=comb[,3],ener=rowSums(eps*tmp))
}

# read in data 
pdbs <- c("~/go_interact/3gbn/gopair/optimize/GA_n50_nn147_nearn56_eij_Zscore/optfiles")
for (p in seq_along(pdbs)){
    pdb <- pdbs[p]
    if (p==1){
        ipars <- read.table(paste(pdb,"/eij.txt",sep=""))
        iparsRes <- ipars$V1
        ipars <- ipars$V2
        
        native <- t(read.table(paste(pdb,"/n_fij.txt",sep="")))
        nonnative <- t(read.table(paste(pdb,"/nn_fij.txt",sep="")))
        nativeRMSD <- read.table(paste(pdb,"/n_rmsd.dat",sep=""))$V2
        nonnativeRMSD <- read.table(paste(pdb,"/nn_rmsd.dat",sep=""))$V2
        
        native <- data.matrix(native[2:nrow(native),])
        nonnative <- data.matrix(nonnative[2:nrow(nonnative),])

        native <- cbind(rep(1,nrow(native)),cbind(rep(p,nrow(native)),cbind(nativeRMSD,native)))
        nonnative <- cbind(rep(0,nrow(nonnative)),cbind(rep(p,nrow(nonnative)),cbind(nonnativeRMSD,nonnative)))
        comb <- as.matrix(rbind(native,nonnative))
        comb <- matrix(as.numeric(comb),ncol=ncol(comb),nrow=nrow(comb),byrow=F)
    } else {
        native <- t(read.table(paste(pdb,"/n_fij.txt",sep="")))
        nonnative <- t(read.table(paste(pdb,"/nn_fij.txt",sep="")))
        nativeRMSD <- read.table(paste(pdb,"/n_rmsd.dat",sep=""))$V2
        nonnativeRMSD <- read.table(paste(pdb,"/nn_rmsd.dat",sep=""))$V2
        native <- data.matrix(native[2:nrow(native),])
        nonnative <- data.matrix(nonnative[2:nrow(nonnative),])

        native <- cbind(rep(1,nrow(native)),cbind(rep(p,nrow(native)),cbind(nativeRMSD,native)))
        nonnative <- cbind(rep(0,nrow(nonnative)),cbind(rep(p,nrow(nonnative)),cbind(nonnativeRMSD,nonnative)))
        tmp <- as.matrix(rbind(native,nonnative))
        tmp <- matrix(as.numeric(tmp),ncol=ncol(tmp),nrow=nrow(tmp),byrow=F)
        comb <- rbind(comb,tmp)
    }
}

test <- comb[comb[,2]==3,]
test <- comb
comb <- comb[comb[,2]!=3,]

# setup parameters for GA
popSize <- 30
iters <- 100000
min <- rep(0,length(ipars))
max <- rep(1,length(ipars))
initialSolution <- matrix(ipars,ncol=length(ipars),nrow=popSize,byrow=F)
tmp_comb <- comb[,c(-1,-2,-3)]
tmp_comb_A <- data.frame(flag=comb[,1],system=comb[,2],rmsd=comb[,3])
m <- nrow(tmp_comb)
n <- ncol(tmp_comb)


# use GA to assign weights
GAReal <- ga(type = "real-valued", fitness=fitnessZ2, min=min, max=max, popSize=popSize, maxiter=iters,suggestions=initialSolution,keepBest=T)

#stop()

Z <- NULL
for (i in 1:iters){
    #get solution
    # modified here: GAReal@bestSol[[i]] to  GAReal@bestSol[[i]][1,]
    Z <- c(Z,calZ(as.vector(GAReal@bestSol[[i]][1,])))
}

# summarize results
Z <- data.frame(iters=1:iters,Z=Z)

train <- comb
comb <- test

# modified here: GAReal@bestSol[[i]] to  GAReal@bestSol[[i]][1,]
bestPars <- as.vector(GAReal@bestSol[[iters]][1,])
check <- calE(bestPars)
check$oldenergy <- calE(ipars)$ener

check <- check[order(check$system,check$ener),]
results <- data.frame(resname= iparsRes, eij_initial=ipars,eij_optimized=as.matrix(as.vector(GAReal@bestSol[[iters]][1,])))
epars <- read.table(paste(pdb,"/eij.txt",sep=""),col.names=c("ijs","old_eij"))
# modified here: GAReal@bestSol[[i]] to  GAReal@bestSol[[i]][1,]
epars$new_eij <- as.vector(GAReal@bestSol[[iters]][1,])
write.table(epars,file="parametersGA.txt",row.names=F,quote=F,col.names=F)
write.table(results,file="oldNewComparison.txt",row.names=F,quote=F,col.names=T)
write.table(check,file="rmsdEnergyComparison.txt",row.names=F,quote=F,col.names=T)
write.table(Z,file="Zscore.txt",row.names=F,quote=F,col.names=T)

