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
#Xreact = read.csv(file="/clusterdata/uqpwains/ibscratch/workdir/scripts/QCscript/48639-non-specific-probes-Illumina450k.csv", 
#                  stringsAsFactors=FALSE)

## Set directory
data="/clusterdata/uqpwains/ibscratch/RAW-IDATS"
directory="/clusterdata/uqpwains/ibscratch/QC-IDATS"

#Set to data directory
setwd(data)

### Load Combined target data
#sample_file="PD_targets_plus_pheno_17May16_955.csv"

## read in targets and idats
#targets=read.csv(sample_file, header=T, stringsAsFactor=F)
#targets$"Sample.ID" <- as.factor(targets$"Sample.ID")


# Change directory to the analysis directory
setwd(directory) 
load("RGset.RObject");
#load("detP_raw.RObject")
load("targets.RObject")
#load("targets_Flt.RObject")
#load("detP.RObject")
load("mSet_raw.RObject")
load("mSetSwFlt.RObject")
#load("gmSetSwFlt.RObject")
#load("gmSetQFlt.RObject")
load("mValsSw.RObject")
load("mValsSq.RObject")
#load("mSetQTFlt.RObject")
#load("betaval.RObject")


########################################################################################

# Density plot 
#jpeg("01_densityplotgender.jpg")
#densityPlot(RGset, sampGroups = targets$Sex, main = "Beta sorted by gender", xlab = "Beta")
#dev.off()
#jpeg("01_densityplotsampleset.jpg")
#densityPlot(RGset, sampGroups = targets$SAMPLE.SET, main = "Beta by sample set", xlab = "Beta")
#dev.off()
#jpeg("01_densityplotphenotype.jpg")
#densityPlot(RGset, sampGroups = targets$Phenotype, main = "Beta sorted by phenotype", xlab = "Beta")
#dev.off()

#Remove Sex and Chromosome 6 from Mset
mSetSwFltpreclean = !(featureNames(mSetSwFlt) %in% ann.450k$Name[ann.450k$chr %in% "chr6"])
mSetSwFlt = mSetSwFlt[mSetSwFltpreclean,]
mSetSwFltpreclean = !(featureNames(mSetSwFlt) %in% ann.450k$Name[ann.450k$chr %in% c("chrX","chrY")])
mSetSwFlt = mSetSwFlt[mSetSwFltpreclean,]

mSet_rawpreclean = !(featureNames(mSet_raw) %in% ann.450k$Name[ann.450k$chr %in% "chr6"])
mSet_raw = mSet_raw[mSet_rawpreclean,]
mSet_rawpreclean = !(featureNames(mSet_raw) %in% ann.450k$Name[ann.450k$chr %in% c("chrX","chrY")])
mSet_raw = mSet_raw[mSet_rawpreclean,]


#MDS plot on raw beta
png("01_MDS_msetrawminus6andsex_sample.png")
mdsPlot(mSet_raw, numPositions = 10000, sampNames = targets$SAMPLE.SET, sampGroups = targets$SAMPLE.SET, 
        pch = 1, pal = brewer.pal(8, "Dark2"), legendPos = "bottomleft",
        main = NULL)

png("01_MDS_msetrawminus6andsex_pheno.png")
mdsPlot(mSet_raw, numPositions = 10000, sampNames = targets$Phenotype, sampGroups = targets$Phenotype, 
        pch = 1, pal = brewer.pal(8, "Dark2"), legendPos = "bottomleft",
        main = NULL)
#
png("01_MDS_msetrawminus6andsex_sex.png")
mdsPlot(mSet_raw, numPositions = 10000, sampNames = targets$Sex, sampGroups = targets$Sex, 
        pch = 1, pal = brewer.pal(8, "Dark2"), legendPos = "bottomleft",
        main = NULL)
#
png("01_MDS_mSetSwFltminus6andsex_sample.png")
mdsPlot(mSetSwFlt, numPositions = 10000, sampNames = targets$SAMPLE.SET, sampGroups = targets$SAMPLE.SET, 
        pch = 1, pal = brewer.pal(8, "Dark2"), legendPos = "bottomleft",
        main = NULL)
#
png("01_MDS_mSetSwFltminus6andsex_pheno.png")
mdsPlot(mSetSwFlt, numPositions = 10000, sampNames = targets$Phenotype, sampGroups = targets$Phenotype, 
        pch = 1, pal = brewer.pal(8, "Dark2"), legendPos = "bottomleft",
        main = NULL)
#
png("01_MDS_mSetSwFltminus6andsex_sex.png")
mdsPlot(mSetSwFlt, numPositions = 10000, sampNames = targets$Sex, sampGroups = targets$Sex, 
        pch = 1, pal = brewer.pal(8, "Dark2"), legendPos = "bottomleft",
        main = NULL)

#
#
#png("01_MDS_betaval_sample.png")
#mdsPlot(betaval, sampNames = targets$SAMPLE.SET, sampGroups = targets$SAMPLE.SET, 
#       main = NULL)
#
#png("01_MDS_betaval_pheno.png")
#mdsPlot(betaval, sampNames = targets$Phenotype, sampGroups = targets$Phenotype, 
#        pch = 1, pal = brewer.pal(8, "Dark2"), legendPos = "bottomleft",
#        main = NULL)
#
#png("01_MDS_betaval_sex.png")
#mdsPlot(betaval, sampNames = targets$Sex, sampGroups = targets$Sex, 
##        pch = 1, pal = brewer.pal(8, "Dark2"), legendPos = "bottomleft",
 #       main = NULL)
#png("01_MDS_mValsSw_sample.png")
#mdsPlot(mValsSw, sampNames = targets$SAMPLE.SET, sampGroups = targets$SAMPLE.SET, 
#        pch = 1, pal = brewer.pal(8, "Dark2"), legendPos = "bottomleft",
#        main = NULL)
#
#png("01_MDS_mValsSw_pheno.png")
#mdsPlot(mValsSw, sampNames = targets$Phenotype, sampGroups = targets$Phenotype, 
#        pch = 1, pal = brewer.pal(8, "Dark2"), legendPos = "bottomleft",
#        main = NULL)
#
#png("01_MDS_mValsSw_sex.png")
#mdsPlot(mValsSw, sampNames = targets$Sex, sampGroups = targets$Sex, 
#        pch = 1, pal = brewer.pal(8, "Dark2"), legendPos = "bottomleft",
#        main = NULL)
#
#
##MDS plot on clean individuals
#png("04_MDS_plot_samplesetFlt.png")
#par(mfrow=c(1,1))
#plotMDS(mValsSw, labels=targets$SAMPLE.SET, col=as.integer(factor(targets$SAMPLE.SET)))
#legend("topleft", legend=c("George_Mellick", "Martin_Kennedy", "CONTROL:_QIMR", "CONTROL:_George Mellick", "CONTROL:_Martin Kennedy"), pch=16,cex=1.2,col=1:2)
#dev.off()

#png("04_MDS_plotFlt.png")
#par(mfrow=c(1,1))
#plotMDS(mValsSw, labels=targets$Phenotype, col=as.integer(factor(targets$Phenotype)))
#legend("topleft", legend=c("PD", "Control"), pch=16,cex=1.2,col=1:2)
#dev.off()


#pdf("02-detection-p-values-after-removal.pdf", width=14)
#barplot(apply(detP,2,mean), main="Mean detection p-values after renoval", col=as.factor(targets$X2D.MATRIX.ID), xaxt="none")
#abline(h=0.001,col="red")
#dev.off()



#jpeg("06_MSetrawdensitybeangroup.jpg")
#densityBeanPlot(mSet_raw, sampGroups = targets[, which(targets$SAMPLE.SET=='George_Mellick' & targets$SAMPLE.SET=='Martin_Kennedy') ], sampNames = targets[, which(targets$SAMPLE.SET=='George_Mellick' & targets$SAMPLE.SET=='Martin_Kennedy') ], main = "Beta per group", pal = brewer.pal(8, "Dark2"), numPositions = 100000)
#dev.off()

#jpeg("06_MSetrawdensitybeanpheno.jpg")
#densityBeanPlot(mSet_raw, sampGroups = targets$Phenotype, sampNames = targets&Phenotype, main = "Beta per phenotype", pal = brewer.pal(8, "Dark2"), numPositions = 100000)
#dev.off()

#jpeg("06_MSetrawdensitybeansex.jpg")
#densityBeanPlot(mSet_raw, sampGroups = targets$Sex, sampNames = targets$Sex, main = "Beta per sex", pal = brewer.pal(8, "Dark2"), numPositions = 100000)
#dev.off()

#jpeg("07_mSetSwFltdensitybeangroup.jpg")
#densityBeanPlot(RGset, sampGroups = targets$SAMPLE.SET, sampNames = targets$SAMPLE.SET, main = "RGset per sample set", pal = brewer.pal(8, "Dark2"), numPositions = 100000)
#dev.off()

#jpeg("07_mSetSwFltdensitybeangroup.jpg")
#densityBeanPlot(RGset, sampGroups = targets$Sex, sampNames = targets$Sex, main = "RGSet per sex", pal = brewer.pal(8, "Dark2"), numPositions = 100000)
#dev.off()

#jpeg("07_mSetSwFltdensitybeangroup.jpg")
#densityBeanPlot(RGset, sampGroups = targets$Phenotype, sampNames = targets$Phenotype, main = "RGSet per phenotype", pal = brewer.pal(8, "Dark2"), numPositions = 100000)
#dev.off()

#jpeg("07_mSetSwFltdensitybeangroup.jpg")
#densityBeanPlot(mSetSwFlt, sampGroups = targets$SAMPLE.SET, sampNames = targets$SAMPLE.SET, main = "Beta per group", pal = brewer.pal(8, "Dark2"), numPositions = 100000)
#dev.off()

#jpeg("07_mSetSwFltdensitybeanpheno.jpg")
#densityBeanPlot(mSetSwFlt, sampGroups = targets$Phenotype, sampNames = targets$Phenotype, main = "Beta per phenotype", pal = brewer.pal(8, "Dark2"), numPositions = 100000)
#dev.off()

#jpeg("07_mSetSwFltdensitybeansex.jpg")
#densityBeanPlot(mSetSwFlt, sampGroups = targets$Sex, sampNames = targets$Sex, main = "Beta per sex", pal = brewer.pal(8, "Dark2"), numPositions = 100000)
#dev.off()

#jpeg("08_mValsSwdensitybeangroup.jpg")
#densityBeanPlot(mValsSw, sampGroups = targets$SAMPLE.SET, sampNames = targets$SAMPLE.SET, main = "Beta per group", pal = brewer.pal(8, "Dark2"), numPositions = 100000)
#dev.off()

#jpeg("08_mValsSwdensitybeanpheno.jpg")
#densityBeanPlot(mValsSw, sampGroups = targets$Phenotype, sampNames = targets$Phenotype, main = "Beta per phenotype", pal = brewer.pal(8, "Dark2"), numPositions = 100000)
#dev.off()

#jpeg("08_mValsSwdensitybeansex.jpg")
#densityBeanPlot(mValsSw, sampGroups = targets$Sex, sampNames = targets$Sex, main = "Beta per sex", pal = brewer.pal(8, "Dark2"), numPositions = 100000)
#dev.off()


#jpeg("05_PCA-mvalQT.jpg")
#princomp(mValsSq)
#dev.off()

#jpeg("05_PCA-mvalswan.jpg")
#princomp(mValsSw)
#dev.off()

#png("03_Normalization_plotQT.png")
#par(mfrow=c(1,2))
#plotBetasByType(MSet.raw [,1], main = "Raw")
#plotBetasByType(mSetQTFlt[,1], main = "QT")
#dev.off()

#png("03_Normalization_plotNORM.png")
#par(mfrow=c(1,2))
#plotBetasByType(mSetQTFlt [,1], main = "QT")
#plotBetasByType(mSetSwFlt[,1], main = "SWAN")
#dev.off()
