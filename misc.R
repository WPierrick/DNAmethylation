data="/clusterdata/uqpwains/ibscratch/RAW-IDATS2"
directory="/clusterdata/uqpwains/ibscratch/QC-IDATS2"

#Set to data directory
setwd(data)

### Load Combined target data
sample_file="PD_Targets2_20160701_1008indiv_pred.csv"

## read in targets and idats
targets=read.csv(sample_file, header=T, stringsAsFactor=F)
#targets$"Sample.ID" <- as.factor(targets$"Sample.ID")


# Change directory to the analysis directory
setwd(directory) 
load("RGset.RObject");
#load("detP_raw.RObject")
load("targets.RObject")
#load("targets_Flt.RObject")
#load("detP.RObject")
load("mSet_raw.RObject")

load("gmSetSwFlt.RObject") # MethylSet with SWAN normalization
directory="/clusterdata/uqpwains/ibscratch/QC-IDATS2"


library(RColorBrewer)

autosomes = !(featureNames(gmSetSwFlt) %in% ann.450k$Name[ann.450k$chr %in% c("chrX","chrY", "chr6")])
mSet_trimmed = gmSetSwFlt[autosomes,]


png("01_MDS_plot_Msettrimm_Sample.png")
mdsPlot(mSet_trimmed , sampNames = targets$SAMPLE.SET, sampGroups = targets$SAMPLE.SET, 
        pch = 1, pal = brewer.pal(8, "Dark2"), legendPos = "bottomleft",
        main = NULL)
dev.off()

targetsG$SAMPLE.ID <- paste(targetsG$Array.Info.Sentrix.ID, targetsG$Array.Info.Sentrix.Position, sep = '_')
outg <- which(as.array(IDremoved)%in%targetsG$SAMPLE_ID)
combi <- merge(test, targetsG, by.x  = "IDremoved", by.y = "SAMPLE.ID", all.y = FALSE)
#combined <- merge(hovprob,annot,by.x="ProbeID",by.y="Name",all.y=TRUE)

targets$Basename <- paste('/ibscratch/wrayvisscher/uqpwains/RAW-IDATS2', sep = '/', targets$SAMPLE.ID)
targets2 <- which(!as.array(IDremoved)%in%targets$SAMPLE.ID)
test2 <- merge(test, targets, by.x = 'IDremoved', by.y = 'SAMPLE.ID', all.y = TRUE)
targets2 <- targets[!(targets$SAMPLE.ID %in% test$IDremoved),]

## read in targets and idats
predsex <- read.table("predictedSex.txt", header=T, stringsAsFactor=F)
targets$predsex <- paste(predsex$predictedSex)

#Change RGset basename as character
RGset@phenoData@data$Basename <- as.character(RGset@phenoData@data$Basename)

#Replace M and F by Male and Female and compare the predicted sex and the actual one
targets$predsex[targets$predsex == "M"] <- "Male"
targets$predsex[targets$predsex == "F"] <- "Female"
targets$GENDER %in% targets$predsex

#Load smoking predictions and add them to target file
smokpred <- read.table('smoking_prediction_in_ChMND_based_on_LBC36.txt', header = TRUE)
targets$smokepred <- paste(smokpred$predicted_prob_of_being_smoker)

#Remove X col
targets$X <- NULL

#Remove X at the begining of the sampleID
agepred$SampleID <- gsub("X","",agepred$SampleID,fixed=TRUE)

celpred <- read.table("Cellcount.txt", header = TRUE)
agepred <- read.csv("PDmeth2_probes_horvath.output.csv")

targets$X <- NULL
targets$dnamage <- agepred$DNAmAge
targets[c('CD8T', 'CD4T','NK', 'Bcell', 'Mono', 'Gran')] <- celpred[c('CD8T', 'CD4T','NK', 'Bcell', 'Mono', 'Gran')]


write.csv(targets, "PD_Targets2_20160701_1008indiv_predfull.csv")

#Histogram
plotsmoking <- ggplot(targets, aes( x = smokepred, fill=GENDER))+ geom_histogram(colour = "black", binwidth = .02) + 
  geom_vline(aes(xintercept = mean(smokepred)), color="red", linetype="dashed", size=1) +
    labs(title="Smoking prediction per phenotype \n using methylation level for cg05575921") + labs(x="Predicted smoking status", y="Count")  + facet_grid( . ~ PHENOTYPE)
plotsmoking

# Boxplot on cell predictions
library(reshape2)
targets5 <- c( "GENDER", "PHENOTYPE", "CD8T","CD4T","NK", "Bcell", "Mono", "Gran")
targets5 <- targets[targets5]
datacell <- melt(targets5, id.var = c("PHENOTYPE", "GENDER"))
plotcell <- ggplot(datacell, aes(x = variable, y = value, color = GENDER)) + geom_boxplot() +    facet_grid(. ~ PHENOTYPE) +  xlab("Cell type") +
  ylab("Proportion")  + ggtitle("Cell type proportion per phenotype") 
plotcell