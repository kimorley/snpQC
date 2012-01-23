# [1] Read in PLINK files from GAPI
# [2] Split their "ID" into constituent plate, well, and sample ID information
# [3] Return standard snpStats list from read.plink with extra "info" object containing above info
# Author: km5
###################################################################################################

plink.input <- function(PREFIX){
	require(snpStats)
	a = read.plink(paste(PREFIX,"bed",sep="."))
	info = do.call(rbind,lapply(a$fam$member, function(x){cbind( member = x, id=paste(unlist(strsplit(x,"_"))[length(strsplit(x,"_")[[1]])-1],unlist(strsplit(x,"_"))[length(strsplit(x,"_")[[1]])],sep="_"), plate=unlist(strsplit(x,"_"))[1], well=unlist(strsplit(x,"_"))[2] )} ) )
	a$info = data.frame(info)
	return(a)
}


