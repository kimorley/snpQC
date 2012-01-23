# [1] read KING output for UNrelated (because no pedigree info in PLINK file)
# [2] reformat IDs to sanger_id format and remove FAMIDs as not in PLINK file
# [3] subset to IDs, number of SNPs, kinship and IBS0 estimates
# Author: km5
###############################################################################

king.input <- function(FILE){
	kin = read.table(FILE, header=T)
	kin$ID1 = do.call(rbind,lapply(as.character(kin$ID1), function(x){paste(unlist(strsplit(x,"_"))[length(strsplit(x,"_")[[1]])-1],unlist(strsplit(x,"_"))[length(strsplit(x,"_")[[1]])],sep="_")}))
	kin$ID2 = do.call(rbind,lapply(as.character(kin$ID2), function(x){paste(unlist(strsplit(x,"_"))[length(strsplit(x,"_")[[1]])-1],unlist(strsplit(x,"_"))[length(strsplit(x,"_")[[1]])],sep="_")}))
	return(subset(kin, select=c(ID1,ID2,N_SNP,Kinship,IBS0)))
}

