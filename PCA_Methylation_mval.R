### Perform PC Analysis ###

setwd("/clusterdata/uqpwains/ibscratch/QC-IDATS")
#load("mValsSw.RObject")

## Set directory
data="/clusterdata/uqpwains/ibscratch/RAW-IDATS"
directory="/clusterdata/uqpwains/ibscratch/QC-IDATS"

#Set to data directory
setwd(data)

### Load Combined target data
#sample_file="PD_targets_final_30May16_954_pred.csv"


## Specify covariates for which methylation data should be adjusted
covariates=c("Phenotype", "SAMPLE.SET", "Sex", "Sentrix.Barcode", "Sample.Section", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "Smoking_prediction", "AGE_crono_plus_bio")
# Add as factor on target file

## read in targets and idats
targets <- read.csv("/ibscratch/wrayvisscher/uqpwains/RAW-IDATS/PD_targets_final_02June16_954_pred_complete.txt", fill = TRUE, header=T, stringsAsFactor=F)

targets$"Sample.ID" <- as.factor(targets$"Sample.ID")



PCobj = prcomp(na.omit(targets[,covariates]), scale. = T, retx=T, center=T)
save(PCobj,file="PCobjCov.RObject")
PCs = PCobj$rotation
save(PCs,file="PCsCov.RObject")

#load("PCsSw.RObject")
#PCsCov <- PCs
#jpeg("04PCA_covariates.jpg")
#plot(PCsCov[,1:2], col="blue")
#dev.off()
#write.table(round(PCs, digits = 4), "PCs.csv",  col.names=T, row.names=T, qu=F, sep=",")











