# Individual, SNP, and plate checks
# Produces:
#	[1] plot of call rate vs heterozygosity (plate colour coding)
#	[2] boxplot of heterozygosity by plate
#	[3] boxplot of call rate by plate
#	[4] table of those individuals outside thresholds for plot [1]
#	[5] comparison of average call rate and heterozygosity by plate
# These are returned as objects in a list 
# Author: km5
###############################################################################

dataCheck <- function(a){
	b = row.summary(a$genotypes)
	b = cbind(b, member = row.names(b))
	c = merge(b, a$info, by = "member")
	c = cbind(c, log10call = log10(c$Call.rate))
	require(ggplot2)
	# Changing ggplot2 plotting defaults
	th = theme_bw()
	th$panel.background = theme_rect(fill = "white", colour = NA)
	theme_set(th)
	square = data.frame(cbind(log10call=c(log10(0.9),log10(1)),Heterozygosity=c(mean(c$Heterozygosity)-3*sd(c$Heterozygosity),mean(c$Heterozygosity)+3*sd(c$Heterozygosity))))
	f1 = qplot(log10call, Heterozygosity, data=c) + geom_point(aes(colour=plate)) + ylab("Heterozygosity \n") + xlab("\n Call rate (log10)") + opts(axis.text.x=theme_text(size=12)) + opts(axis.text.y=theme_text(size=12)) + opts(axis.title.x=theme_text(size=12)) + opts(axis.title.y=theme_text(size=12,angle=90)) + geom_vline(xintercept=log10(0.9),colour="#800517",linetype=2) + geom_hline(yintercept=(mean(c$Heterozygosity)-3*sd(c$Heterozygosity)),colour="#800517",linetype=2) + geom_hline(yintercept=(mean(c$Heterozygosity)+3*sd(c$Heterozygosity)),colour="#800517",linetype=2) + geom_ribbon(data=square,aes(ymin=(mean(c$Heterozygosity)-3*sd(c$Heterozygosity)),ymax=(mean(c$Heterozygosity)+3*sd(c$Heterozygosity))),alpha=0.1,fill="#41A317") + opts(axis.text.x=theme_text(size=8)) + opts(axis.text.y=theme_text(size=8)) + opts(axis.title.x=theme_text(size=8)) + opts(axis.title.y=theme_text(size=8,angle=90))
	f2 = qplot(plate, Heterozygosity, data = c, geom = "boxplot", fill = plate) + coord_flip() + opts(legend.position="none") + ylab("\n Heterozygosity") + xlab("Plate")
	f3 = qplot(plate, Call.rate, data = c, geom = "boxplot", fill = plate) + coord_flip() + opts(legend.position="none") + ylab("\n Call rate") + xlab("Plate")
	t1 = subset(c, Heterozygosity < (mean(c$Heterozygosity)-3*sd(c$Heterozygosity)) | Heterozygosity > (mean(c$Heterozygosity)+3*sd(c$Heterozygosity)) | Call.rate < 0.9)
	t1 = subset(t1, select=c(id, plate, well, Call.rate, Heterozygosity))	
	t2 = data.frame(t(rbind(table(c$plate),tapply(c$Call.rate, c$plate, mean), tapply(c$Heterozygosity, c$plate, mean))))
	names(t2) = c("samples","ind.call.rate","ind.heterozygosity")
	temp = data.frame()
	for (x in unique(a$info$plate)){
		d = subset(a$info, plate == x, select = member)
		use = a$fam$member %in% d$member
		e = a$genotypes[use,]
		f = col.summary(e)
		f = cbind(f,plate = x)
		temp = rbind(temp,f)
	}
	t0 = data.frame(tapply(temp$Call.rate, temp$plate, mean))
	names(t0) = c("snp.call.rate")
	t2 = cbind(t2,t0)
	return(list(fig1 = f1, fig2 = f2, fig3 = f3, tab1 = t1, tab2 = t2))
}


