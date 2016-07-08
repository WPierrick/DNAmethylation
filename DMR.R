### RUV Analysis for Parkinson disease methylation samples ###
### Modified by Pierrick Wainschtein from Beben Benyamin from Maksimovic et al 2015 ####

#source("http://bioconductor.org/biocLite.R")

#biocLite("sva")
#biocLite("ROCR")
#biocLite("isva")
#biocLite("FlowSorted.Blood.450k")


library(missMethyl)
library(limma)
library(minfi)
#library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#library(sva)
#library(ruv)
#library(ROCR)
#library(RColorBrewer)
#library(matrixStats)
#library(isva)
#library(FlowSorted.Blood.450k)
#source("/clusterdata/uqpwains/ibscratch/workdir/scripts/QCscript/utils.R")
#source("/clusterdata/uqpwains/ibscratch/workdir/scripts/QCscript/functions.R")
#ann.450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#Xreact = read.csv(file="/clusterdata/uqpwains/ibscratch/workdir/scripts/QCscript/48639-non-specific-probes-Illumina450k.csv", 
#                  stringsAsFactors=FALSE)

## Set directory
directory="/clusterdata/uqpwains/ibscratch/QC-IDATS"

#Set to data directory
setwd(directory) 
load("targets_Flt.RObject")
load("gmSetQFlt.RObject")
#load("invN.RObject")

M <- getM(gmSetQFlt, betaThreshold = 0.001)
dmp <- dmpFinder(M, pheno = targets$Phenotype, type = "categorical")
head(dmp)
save(dmp, file="dmp.RObject")
cpgs <- rownames(dmp)

save(cpgs, file = 'DMR_cpgs.RObject')
jpeg("09_DMR.jpg")
par (mfrow = c(2,2))
plotCpg(gmSetQFlt, cpg = cpgs, pheno = targets$Phenotype)
dev.off()