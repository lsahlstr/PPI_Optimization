enrichment <- function(tt,thresholds=seq(0.05,0.20,0.05)){
    enrich <- NULL
    for (threshold in thresholds) {
		nhits <- floor(threshold*nrow(tt))
		tt$index_energy <- rownames(tt)[order(tt$ener)]
		tt$index_rmsd <- rownames(tt)[order(tt$rmsd)]
		enrich <- c(enrich,sum(tt$index_energy[1:nhits] %in% tt$index_rmsd[1:nhits])/((threshold^2)*nrow(tt)))
	}
	enrich <- as.data.frame(t(enrich))
	colnames(enrich) <- paste("ES_",thresholds,sep="")
	enrich$pdb <- as.character(pdbs[unique(tt$system)])
	return(enrich)
}

ddply(.dat=rmsd_ener,.var=c("system"),.fun=enrichment)