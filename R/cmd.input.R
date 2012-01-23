# Collects command-line arguments for family checking
# Author: km5
###############################################################################

cmd.input <- function(){
	ARGS <-Sys.getenv(c("K","KPATH","D","DPATH","M","MPATH"))
	FAM = data.frame(cbind(rbind(as.character(ARGS["K"]),as.character(ARGS["D"]),as.character(ARGS["M"])),rbind(as.character(ARGS["KPATH"]),as.character(ARGS["DPATH"]),as.character(ARGS["MPATH"]))))
	return(FAM)	
}

