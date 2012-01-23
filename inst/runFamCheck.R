# Helper script for running SNP array family check
# This will produce a Sweave document in the working directory
# Author: km5
###############################################################################

require(tools)
require(snpQC)

FAM = cmd.input()

FILE = system.file("snpQCfam.Rnw",package = "snpQC")
Sweave(FILE)
texi2dvi("snpQCfam.tex",pdf=T)

file.remove("snpQCfam-fig1.pdf","snpQCfam.tex","snpQCfam.aux","snpQCfam.log")
q(save = "no")



