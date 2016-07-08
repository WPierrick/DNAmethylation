#Getbeta

library(missMethyl)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(sva)
library(ruv)
library(ROCR)
library(RColorBrewer)
library(matrixStats)
library(isva)
library(FlowSorted.Blood.450k)
source("/clusterdata/uqpwains/ibscratch/workdir/scripts/QCscript/utils.R")
source("/clusterdata/uqpwains/ibscratch/workdir/scripts/QCscript/functions.R")
ann.450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
Xreact = read.csv(file="/clusterdata/uqpwains/ibscratch/workdir/scripts/QCscript/48639-non-specific-probes-Illumina450k.csv", 
                  stringsAsFactors=FALSE)

## Set directory

setwd ("/clusterdata/uqpwains/ibscratch/QC-IDATS")
load("/ibscratch/wrayvisscher/uqpwains/QC-IDATS/betaval.RObject")
roh <- rownames(betaval)
roh <- t(roh)
roh <- t(roh)
betavalFlt <- as.data.frame(cbind(roh, betaval))
save(betavalFlt, file = "betavalFlt.RObject")