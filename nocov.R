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

## Set directory
data="/clusterdata/uqpwains/ibscratch/RAW-IDATS"
directory="/clusterdata/uqpwains/ibscratch/QC-IDATS"

#Set to data directory
setwd(data)


# Change directory to the analysis directory
setwd(directory) 
#save(RGset, file="RGset.RObject")
#save(detP_raw, file="detP_raw.RObject")
#save(targets, file="targets_raw.RObject")
#load("RGset.RObject")
#load("detP_raw.RObject")
#load("detP.RObject")
#load("targets_raw.RObject")
#load("targets.RObject")
#load("mSet_raw.RObject")
#load("mSetSwFlt.RObject")
#load("gmSetSwFlt.RObject")
#load("gmSetQFlt.RObject")
#load("mValsSw.RObject")
load("mValsSq.RObject")
load("targets_Flt.RObject")



#Inverse normal transformation
invN=apply(mValsSq,1, function(x) InvNorm(x))
invN=t(invN)

#Change invN into mValsSq
mValsSq <- invN
colnames(mValsSq) <- targets$Sample.ID
all(targets$Sample.ID == colnames(mValsSq))

dmps = list(swan=list(),sqn=list())
Disease = factor(targets$Phenotype)
design = model.matrix(~Disease) ## Design matrix


#Compare different analysis model
##Linear model
assoc_results = apply(mValsSq,1,function(x) case_ctrl_assoc(x, targets$Phenotype));
assoc_results=as.data.frame(matrix(unlist(assoc_results), ncol=5, byrow=T));
colnames(assoc_results)=c("adjR2","effect","se","t" ,"pval");
assoc_results$Probe=rownames(mValsSq)
#save(assoc_results, file="assoc_results.RObject");

p <- assoc_results$pval
png("LM-no adj.png")
qqplot(p,"LM on M values")
dev.off()

##limma
dmps$sqn$limma = limmaFit(data=mValsSq,design=design,coef=2)
limma <- dmps$sqn$limma
save(limma, file="limma.RObject")

p <- limma$P.Value
png("limma-no adj.png")
qqplot(p,"LIMMA on M values")
dev.off()

