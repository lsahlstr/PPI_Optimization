#######################################################################
# Routine to compute enrichment score. data = rmsd_ener from main.R
#######################################################################
enrichment <- function(data,thresholds=seq(0.05,0.20,0.05)){
    enrich <- NULL
    for (threshold in thresholds) {
		nhits <- floor(threshold*nrow(data))
		data$index_energy <- rownames(data)[order(data$ener)]
		data$index_rmsd <- rownames(data)[order(data$rmsd)]
		enrich <- c(enrich,sum(data$index_energy[1:nhits] %in% data$index_rmsd[1:nhits])/((threshold**2)*nrow(data)))
	}
	enrich <- as.data.frame(t(enrich))
	colnames(enrich) <- paste("ES_",thresholds,sep="")
	enrich$pdb <- as.character(pdbs[unique(data$system)])
	return(enrich)
}