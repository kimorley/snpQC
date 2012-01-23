# Family checks
# Produces:
#	[1] plot of kinship coefficient vs IBS0
#	[2] table summarising pair-wise comparisons
#	[3] flag:
#		0 = OK
#		1 = Father unrelated
#		2 = Mother unrelated
#		3 = All unrelated
#		4 = Potential inbreeding
#		5 = Sample swap
# These are returned as objects in a list 
# Author: km5
###############################################################################

famCheck <- function(FAM){	# Assuming ID list order is offspring, father, mother
	data(prunedSnps)
	count = 0
	FAM = cbind(FAM,IN=0)
	# Extract each family member from PLINK files
	for (i in 1:3){	
		if (FAM[i,1]!="0"){
			count = count+1
			FAM$IN[i] = count
			CMD <- paste("grep ",FAM[i,1]," ",FAM[i,2],".fam | awk '{ print $1,$2}'",sep="")
			ID <- read.table(pipe(CMD), colClasses="character")
			write.table(ID, file="keeplist", col.names=F, row.names=F, quote=F)
			write.table(pruned.snps, file="extractlist", col.names=F, row.names=F, quote=F)
			CMD = paste("plink --noweb --allow-no-sex --bfile ",FAM[i,2]," --extract extractlist --keep keeplist --make-bed --out ",FAM[i,2],"-",FAM[i,1],sep="")
			system(CMD)
		}
	}
	# Merge PLINK files of family members
	if (count == 2){	# Duo
		mergelist = data.frame(rbind(paste(FAM[which(FAM$IN==2),2],"-",FAM[which(FAM$IN==2),1],".bed",sep=""),paste(FAM[which(FAM$IN==2),2],"-",FAM[which(FAM$IN==2),1],".bim",sep=""),paste(FAM[which(FAM$IN==2),2],"-",FAM[which(FAM$IN==2),1],".fam",sep="")))
		write.table(t(mergelist), file="mergelist", col.names=F, row.names=F, quote=F)
		CMD = paste("plink --noweb --allow-no-sex --bfile ",FAM[1,2],"-",FAM[1,1]," --merge-list mergelist --make-bed --out temp",sep="")
		system(CMD)
		file.remove(paste(FAM[1,2],"-",FAM[1,1],".bed",sep=""),paste(FAM[1,2],"-",FAM[1,1],".bim",sep=""),paste(FAM[1,2],"-",FAM[1,1],".fam",sep=""),paste(FAM[1,2],"-",FAM[1,1],".log",sep=""),paste(FAM[1,2],"-",FAM[1,1],".nof",sep=""),paste(FAM[1,2],"-",FAM[1,1],".nosex",sep=""))
		file.remove(paste(FAM[which(FAM$IN==2),2],"-",FAM[which(FAM$IN==2),1],".bed",sep=""),paste(FAM[which(FAM$IN==2),2],"-",FAM[which(FAM$IN==2),1],".bim",sep=""),paste(FAM[which(FAM$IN==2),2],"-",FAM[which(FAM$IN==2),1],".fam",sep=""),paste(FAM[which(FAM$IN==2),2],"-",FAM[which(FAM$IN==2),1],".log",sep=""),paste(FAM[which(FAM$IN==2),2],"-",FAM[which(FAM$IN==2),1],".nof",sep=""),paste(FAM[which(FAM$IN==2),2],"-",FAM[which(FAM$IN==2),1],".nosex",sep=""))
	}else if (count == 3){ # Trio
		mergelist = data.frame(cbind(rbind(paste(FAM[2,2],"-",FAM[2,1],".bed",sep=""),paste(FAM[2,2],"-",FAM[2,1],".bim",sep=""),paste(FAM[2,2],"-",FAM[2,1],".fam",sep="")),rbind(paste(FAM[3,2],"-",FAM[3,1],".bed",sep=""),paste(FAM[3,2],"-",FAM[3,1],".bim",sep=""),paste(FAM[3,2],"-",FAM[3,1],".fam",sep=""))))
		write.table(t(mergelist), file="mergelist", col.names=F, row.names=F, quote=F)
		CMD = paste("plink --noweb --allow-no-sex --bfile ",FAM[1,2],"-",FAM[1,1]," --merge-list mergelist --make-bed --out temp",sep="")
		system(CMD)
		for (i in 1:3){
			file.remove(paste(FAM[i,2],"-",FAM[i,1],".bed",sep=""),paste(FAM[i,2],"-",FAM[i,1],".bim",sep=""),paste(FAM[i,2],"-",FAM[i,1],".fam",sep=""),paste(FAM[i,2],"-",FAM[i,1],".log",sep=""),paste(FAM[i,2],"-",FAM[i,1],".nof",sep=""),paste(FAM[i,2],"-",FAM[i,1],".nosex",sep=""))
		}
	}
	# Calculate kinship
	CMD = paste("/nfs/ddd0/software/src/king -b temp.bed --kinship",sep="")
	system(CMD)
	kin = king.input("king.kin0")
	kin = cbind(kin,Rel=0,Err=0)
	data(kinBounds)
	IDLIST = as.vector(FAM[,1])
	for (j in 1:length(kin[,1])){ 
		if (sum(kin[j,1:2] %in% IDLIST[2:3]) == 2){ # This is the parent-parent estimate
			kin$Rel[j] = 1	
			if (kin$Kinship[j] > kinBounds$Kinship[4]){	# If the kinship coeff is greater than the lower bound for 2nd degree relatives, flag it
				kin$Err[j] = 1
			}
		}else{	# These are parent-offspring estimates
			if (kin$Kinship[j] < kinBounds$Kinship[2]){	# If the kinship coeff is less than the lower bound for P-O relationships, flag it
				kin$Err[j] = 1
			}
		}
	}
	# Check for errors
	FLAG = 0
	if (sum(kin$Err) > 0){									# There are errors with pedigree relationships
		if (1 %in% kin$Rel){ # Trio
			if (kin$Err[which(kin$Rel==1)] == 1){				# Parents are too related
				if (sum(kin$Err[which(kin$Rel==0)]) == 0){		# Appropriate P-O relationships
					FLAG = 4									# Inbreeding
				}else{
					FLAG = 5									# Potential within-family sample swaps
				}
			}else{
				if (sum(kin$Err[which(kin$Rel==0)]) == 2){		# Both parents unrelated to child
					FLAG = 3
				}else if(sum(kin$Err[which(kin$Rel==0)]) == 1){	# One parent unrelated
					if (IDLIST[2] %in% kin[which(kin$Rel==0 & kin$Err==1),1:2]){	# Father unrelated
						FLAG = 1
					}else{										# Mother unrelated
						FLAG = 2
					}
				}
			}
		}else{	# Duo
			if (IDLIST[2] %in% kin[which(kin$Rel==0 & kin$Err==1),1:2]){	# Father unrelated
				FLAG = 1
			}else{										# Mother unrelated
				FLAG = 2
			}
		}
	}
	write.table(FLAG, "flag.txt",col.names=F, row.names=F, quote=F)
	# Plot
	require(ggplot2)
	th = theme_bw()
	th$panel.background = theme_rect(fill = "white", colour = NA)
	theme_set(th)
	PP = data.frame(cbind(Kinship=c(-0.25,kinBounds$Kinship[4]),IBS0=c(0,1)))
	PO = data.frame(cbind(Kinship=c(kinBounds$Kinship[1],kinBounds$Kinship[2]),IBS0=c(0,kinBounds$IBS0[2])))
	f1 = qplot(Kinship,IBS0, data=kin, xlim=c(-0.25,0.5), ylim=c(0,1)) + geom_ribbon(data=PO,aes(ymin=0,ymax=kinBounds$IBS0[2]),alpha=0.1,fill="#C12869") + geom_ribbon(data=PP,aes(ymin=0,ymax=1),alpha=0.1,fill="#6698FF") + geom_hline(yintercept=kinBounds$IBS0,colour="grey70",linetype=2) + geom_vline(xintercept=kinBounds$Kinship,colour="grey70",linetype=2) + geom_point(aes(colour=as.factor(Rel),size=3,alpha=0.4)) + xlab("\n Kinship coefficient") + opts(legend.position="none")  
	file.remove("keeplist","extractlist","temp.nof","temp.fam","temp.bim","temp.bed","temp.log","kingTMP.ped","kingTMP.dat","king.kin","king.kin0")
	return(list(tab1 = kin, fig1 = f1, flag = FLAG))
}


