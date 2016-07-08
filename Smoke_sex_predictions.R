### RUV Analysis for Parkinson disease methylation samples ###
### Modified by Pierrick Wainschtein from Beben Benyamin from Maksimovic et al 2015 ####

#source("http://bioconductor.org/biocLite.R")

#biocLite("sva")
#biocLite("ROCR")
#biocLite("isva")
#biocLite("FlowSorted.Blood.450k")

.libPaths("/ibscratch/wrayvisscher/uqpwains/R/x86_64-pc-linux-gnu-library/3.2")

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
data="/clusterdata/uqpwains/ibscratch/RAW-IDATS2"
directory="/clusterdata/uqpwains/ibscratch/QC-IDATS2"

#Set to data directory
setwd(data)

## read in targets and idats
targets <- read.csv("PD_Targets2_20160701_1008indiv_pred.csv", fill = TRUE, header=T, stringsAsFactor=F)
RGset <- read.450k.exp(targets=targets, recursive=FALSE);
 ## calculate detection p-values

# Change directory to the analysis directory
setwd(directory) 
save(RGset, file="RGset.RObject")

# Script for PCA methylation and cell count / smoking predictors

#source("http://bioconductor.org/biocLite.R")
#biocLite("FlowSorted.Blood.450k")

.libPaths("/ibscratch/wrayvisscher/uqpwains/R/x86_64-pc-linux-gnu-library/3.2")


library(FlowSorted.Blood.450k)
library(minfi)

setwd("/clusterdata/uqpwains/ibscratch/QC-IDATS2")


###################################
#1. Sex prediction

load("RGset.RObject") # RGset (raw data) is obtained using minfi
#RGset$Slide <- as.numeric(RGset$Slide)

### Load Combined target data
#sample_file="PD_Targets2_20160628_1008indiv.csv"
#targets=read.csv(sample_file, header=T, stringsAsFactor=T)
#targets$"Sample.ID" <- as.factor(targets$"Sample.ID")


###
#setwd("/ibscratch/wrayvisscher/uqpwains/workdir/output/02-Predictions/Cellcount")

#MSet <- preprocessRaw(RGset)
#GMsetEx <- mapToGenome(MSet)
#estSex <- getSex(GMsetEx)
#write.table(estSex, file = "predictedSex.txt", col.names=T, row.names=T, qu=F)
#GMsetEx <- addSex(GMsetEx, sex = estSex)
#save(GMsetEx, file="GenomicMethylSetRaw.RObject")
#png("Sex_plot.png")
#plotSex(estSex, id = targets$SAMPLE.ID)
#dev.off()

####################################
#2.Cell counts proportions (Houseman et al 2012)

counts <- estimateCellCounts(RGset, compositeCellType = "Blood",
                             cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Gran"),
                             returnAll = FALSE, meanPlot = FALSE, verbose = TRUE)
round(counts, 2)
save(counts, file = "Cell_Counts.txt", col.names=T, row.names=T, qu=F) # save the data
write.table(counts, file = "Cellcount.txt", col.names=T, row.names=T, qu=F)
#####################################
#3. Smoking prediction

## R object of methylation beta values after QCs (e.g. background corrected and normalised and removing bad probes/samples (i.e call rate < 0.95)) to be used for prediction in your sample
#load("gmSetSwFlt.RObject")
#betavalFlt <- getBeta(gmSetSwFlt)
#save(betavalFlt, file = "betavalFlt.RObject")
#betavalFlt="/clusterdata/uqpwains/ibscratch/QC-IDATS/betavalFlt.RObject"

##Output file 
#outfile="smoking_prediction_in_ChMND_based_on_LBC36.txt"

## load mehtylation data in prediction samples
#load(betavalFlt)

## rename the meth data object to 'pred_data'
## So if the methylation data is called 'beta_filt' then change as follows:
#pred_data=betavalFlt
#rm(meth_data)

## Load the regression object on Lothian Birth Cohort 1936 (LBC36) data 
## This data was used to train smoking predictor
## cg05575921 is the single most associated probe with smoking status
#load("/ibscratch/wrayvisscher/sonia/methylation/Final_2460_Data/Smoking_Prediction/LBC36_cg05575921_train_reg.RObject")

#pred_data=t(pred_data)
#pred_data=as.data.frame(pred_data[,"cg05575921"])
#pred_data$smoke_never_current=NA
#colnames(pred_data)[1]="cg05575921"

#pred_data$smoke_never_current=as.integer(pred_data$smoke_never_current)
#pred_data$smoke_never_current=predict(reg, pred_data, type="response")
#pred_data$ID=rownames(pred_data)

#colnames(pred_data)[2]="predicted_prob_of_being_smoker"
#write.table(pred_data,outfile, quote=F, row.names=F, sep="\t")

#Remove sex and Chr 6
autosomes = !(featureNames(gmSetSwFlt) %in% ann.450k$Name[ann.450k$chr %in% c("chrX","chrY", "chr6")])
mSet_trimmed = gmSetSwFlt[autosomes,]


#MDS plot on Trimmed MSet
png("01_MDS_plot_trimmedMset_Sample.png")
mdsPlot(mSet_trimmed , sampNames = targets$SAMPLE.SET, sampGroups = targets$SAMPLE.SET, 
        pch = 1, pal = brewer.pal(8, "Dark2"), legendPos = "bottomleft",
        main = NULL)

png("01_bisMDS_plot_trimmedMset_Phenotype.png")
par(mfrow=c(1,1))
plotMDS(mSet_trimmed , labels=targets$PHENOTYPE, col=as.integer(factor(targets$PHENOTYPE)))
legend("topleft", legend=c("PD", "Control"), pch=16,cex=1.2,col=1:2)
dev.off()







