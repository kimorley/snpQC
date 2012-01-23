# Helper script for running SNP array data check
# This will produce a Sweave document in the working directory
# Author: km5
###############################################################################

ARGS <-Sys.getenv(c("PREFIX"))
PREFIX <- as.character(ARGS["PREFIX"])

require(tools)
require(snpQC)

FILE = system.file("snpQCdata.Rnw",package = "snpQC")
Sweave(FILE)
texi2dvi("snpQCdata.tex",pdf=T)

file.remove("snpQCdata-fig1.pdf","snpQCdata-fig2.pdf","snpQCdata-fig3.pdf","snpQCdata.tex","snpQCdata.aux","snpQCdata.log")
q(save = "no")

