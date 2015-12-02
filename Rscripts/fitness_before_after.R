#######################################################################
# Z-score and SLR for each system before and after optimization
#######################################################################
fitness_before_after <- function() {

	zscore_indiv_old <- NULL
	zscore_indiv_new <- NULL
	slr_indiv_old <- NULL
	slr_indiv_new <- NULL

	if (potFlag == "lj") {

		zscore_indiv_old <- fitnessZ_lj_indiv(ipars)
		zscore_indiv_new <- fitnessZ_lj_indiv(bestPars)
		zscore_indiv_write <- data.frame(old=zscore_indiv_old,new=zscore_indiv_new)
		zscore_indiv_write$pdb <- pdbs
	
		slr_indiv_old <- fitnessSLR_lj_indiv(ipars)
		slr_indiv_new <- fitnessSLR_lj_indiv(bestPars)
		slr_indiv_write <- data.frame(old=slr_indiv_old,new=slr_indiv_new)
		slr_indiv_write$pdb <- pdbs

	} else if (potFlag == "eten") {

		zscore_indiv_old <- fitnessZ_eten_indiv(ipars)
		zscore_indiv_new <- fitnessZ_eten_indiv(bestPars)
		zscore_indiv_write <- data.frame(old=zscore_indiv_old,new=zscore_indiv_new)
		zscore_indiv_write$pdb <- pdbs
	
		slr_indiv_old <- fitnessSLR_eten_indiv(ipars)
		slr_indiv_new <- fitnessSLR_eten_indiv(bestPars)
		slr_indiv_write <- data.frame(old=slr_indiv_old,new=slr_indiv_new)
		slr_indiv_write$pdb <- pdbs

	} else if (potFlag == "etsr") {

		zscore_indiv_old <- fitnessZ_etsr_indiv(ipars)
		zscore_indiv_new <- fitnessZ_etsr_indiv(bestPars)
		zscore_indiv_write <- data.frame(old=zscore_indiv_old,new=zscore_indiv_new)
		zscore_indiv_write$pdb <- pdbs
	
		slr_indiv_old <- fitnessSLR_etsr_indiv(ipars)
		slr_indiv_new <- fitnessSLR_etsr_indiv(bestPars)
		slr_indiv_write <- data.frame(old=slr_indiv_old,new=slr_indiv_new)
		slr_indiv_write$pdb <- pdbs
	
	}
	
	write.table(zscore_indiv_write,file="zscore_indiv.dat",row.names=F,quote=F,col.names=F)
	write.table(slr_indiv_write,file="slr_indiv.dat",row.names=F,quote=F,col.names=F)

}