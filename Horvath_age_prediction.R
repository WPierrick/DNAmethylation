# Script for file curation to use with Horvath online calculator
# Pierrick Wainschtein

library(data.table)

# 1) Set working directory and loading object
# Starting files are annotation file and a beta values RObject (matrix) directly outputed from the getBeta() in Minfi package

setwd ("/ibscratch/wrayvisscher/uqpwains/QC-IDATS2")
load("/ibscratch/wrayvisscher/uqpwains/QC-IDATS/betavalraw.RObject")
datMiniAnnotation=read.csv(file = "/clusterdata/uqpwains/ibscratch/workdir/scripts/02-Predictions/Hovarth script/datMiniAnnotation.csv")

# 2) Change the matrix to make it a data frame and make sure the first column are cpg sites
# Method 1 :
betavalAge <- as.data.frame(setDT(as.data.frame(BetaSwFlt), keep.rownames = TRUE)[])
colnames(betavalAge)[1] <- "ProbeID"

# Method 2, use only if data.table is not working
#betamodif <- rownames(betaval)
#betamodif <- t(betamodif)
#betamodif <- t(betamodif)
#betavalAge <- as.data.frame(cbind(betamodif, betaval))
#names(betavalAge) <- "Sample.ID"

# 3) Format the annotation file according to what is required by Horvath Calculator

# Method 1 (recommended):
# Starting from betavalAge
head(betavalAge[,1:5])
# ProbeID 100973270083_R01C01 100973270083_R01C02 100973270083_R02C01 100973270083_R02C02
# 1 cg00050873          0.80381106           0.4749340           0.5139241           0.5907928
# 2 cg00212031          0.04156652           0.8060606           0.7307692           0.4681373
datReduced <- betavalAge[which(betavalAge[,"ProbeID"]%in%as.character(datMiniAnnotation[,1])),]
# Output file number 1:
write.table(datReduced, "PD2meth_probes_predictor_horvath.csv", row.names=F, sep="," )


## Drop hte "X" in colnames
cnames <- colnames(datReduced)
cnames <- gsub("X","",cnames,fixed=TRUE)
colnames(datReduced) <- cnames

combined <- merge(datReduced,datMiniAnnotation,by.x="ProbeID",by.y="Name",all.y=TRUE)
cnames <- colnames(combined)
cnames <- gsub("X","",cnames,fixed=TRUE)
colnames(combined) <- cnames


#table(sort(combined[,"ProbeID"]==sort(datMiniAnnotation[,"Name"])))
#class(combined[,1])
#table(sort(as.character(combined[,"ProbeID"])==sort(as.character(datMiniAnnotation[,"Name"]))))
#summary(sort(as.character(datMiniAnnotation[,1]))==sort(as.character(combined[,1])))
#combinedFiltered <- combined[,colnames(datReduced)]
#head(combinedFiltered[,1:5])
#cnames <- colnames(combinedFiltered)
#> cnames <- gsub("X","",cnames,fixed=TRUE)
#> colnames(combinedFiltered) <- cnames
#> head(combinedFiltered[,1:5])
#pheno <- read.csv("LBC1921_annot.csv")
#> head(pheno)
#c("ProbeID",as.character(pheno[,"SampleID"]))
#length(intersect(as.character(pheno[,"SampleID"]),colnames(combinedFiltered)))
#cnames <- c("ProbeID",as.character(pheno[,"SampleID"]))
#head(colnames(combined))
#combined <- merge(datReduced,datMiniAnnotation,by.x="ProbeID",by.y="Name",all.y=TRUE)
#> cnames <- colnames(combined)
#> cnames <- gsub("X","",cnames,fixed=TRUE)
#> colnames(combined) <- cnames
#> head(combined[,1:5])
#table(cnames%in%colnames(combined))
#> which(!cnames%in%colnames(combined))
#> cnames[which(!cnames%in%colnames(combined))]
#> "101070190063_R03C01"%in%pheno[,"SampleID"]
#> "101070190063_R03C01"%in%colnames(combined)

## Verif
## table(sort(as.character(datMiniAnnotation[,1]))==sort(as.character(combined[,1])))
## If everything worked well, it should be TRUE for everything
#combinedFiltered <- combined[,colnames(combinedFiltered)]
#cnames <- colnames(combinedFiltered)
#cnames <- gsub("X","",cnames,fixed=TRUE)
#colnames(combinedFiltered) <- cnames

cnames <- c("ProbeID",as.character(annotfile[,"SampleID"]))
combinedFiltered <- combined[,cnames]
write.csv(combinedFiltered, file = "PDmeth2_probes_predictor_horvath_filtered.csv",col.names=TRUE,row.names=TRUE,quote=FALSE)


#Method 2 from Horvath tutorial, not sure it's working:
#dat0 <- betavalAge
#match1=match(datMiniAnnotation[,1], rownames(dat0) )
#dat0Reduced=dat0[match1,]
#dat0Reduced[,1]=as.character(dat0Reduced[,1])
#dat0Reduced[is.na(match1),1]= as.character(datMiniAnnotation[is.na(match1),1])
#datout=data.frame(dat0Reduced)
# make sure you output numeric variables...
#for (i in 2:dim(datout)[[2]]  ){datout[,i]= as.numeric(as.character(gsub(x=datout[,i],pattern="\"",replacement=""))) }
#replace "MethylationData" with a filename of your choice
#write.table(datout, "PDmeth_probes_predictor_horvath.csv", row.names=F, sep="," )


# 4) Create a file containing phenotype using annotation file from Horvath website:

annotfile <- data.frame(SampleID = targets$SAMPLE.ID, Age = 'NA', tissue = "BloodWB", Female = targets$GENDER)
annotfile$Female <- as.character(annotfile$Female)

annotfile$Female[annotfile$Female == "Male"] <- "0"
annotfile$Female[annotfile$Female == "Female"] <- "1"

# Out put file number 2:
write.table(annotfile, file="PD2_annotfile.csv", sep=",", row.names=F, quote=F)

# 5) Sent output files 1 and 2 on the online calculator
# 6) Enjoy Australia
# 7) When the results are back, load targets and prediction from Horvath calculator files

targetsfull=read.csv("PD_targets_final_02June16_954_pred_complete.csv", header=T, stringsAsFactor=F)
targets$"Sample.ID" <- as.factor(targets$"Sample.ID")
predic=read.csv("PDmeth_probes_predictor_horvath_filtered.output.csv", header=T, stringsAsFactor=F)

# 8) Make sure the name of first column of predict is Sample.ID
names(predic)[1] <- "Sample.ID"

# 9) Remove the X at the beginning of sample IDs (R can add a X)
cnames <- predic[,1]
cnames <- gsub("X","",cnames,fixed=TRUE)
predic[,1] <- cnames

# 10) Merge the predictions on the target file
targets <- merge(targets, predic[,c(1,2)], by="Sample.ID", sort = FALSE, all.x = TRUE)
#write.csv(targets, "newTargetname.csv")
targets <- targets[-235, ]

# 11) Plot correlations between chronological age and predicted age
corrcomplete <-  cor(targets$Age_both, targets$DNAmAge, use="complete.obs")
plot(targets$Age_both,targets$DNAmAge, main = "Predicted age vs chronological age", xlab = "Chronological age", ylab = "Predicted age (Horvath)")
abline(0,1,lwd=3)
abline(lm(targets$Age_both ~ targets$DNAmAge), col = 4 )

targets.sub <- subset(targets, targets$Phenotype == "PD")
corrpd <- cor(targets.sub$Age_both,targets.sub$DNAmAge, use="complete.obs")
plot(targets.sub$Age_both,targets.sub$DNAmAge, main = "Predicted age vs chronological age in PD", xlab = "Chronological age", ylab = "Predicted age (Horvath)", text(corrpd))
abline(0,1,lwd=3)
abline(lm(targets.sub$Age_both ~ targets.sub$DNAmAge), col = 4 )

targets.control <- subset(targets, targets$Phenotype == "Control")
corrcontrol <- cor(targets.control$Age_both, targets.control$DNAmAge, use="complete.obs")
plot(targets.control$Age_both, targets.control$DNAmAge, main = "Predicted age vs chronological age in Control", xlab = "Chronological age", ylab = "Predicted age (Horvath)")
abline(0,1,lwd=3)
abline(lm(targets.control$Age_both ~ targets.control$DNAmAge), col = 4 )

targets.male <- subset(targets, targets$Sex == "Male")
corrmale <- cor(targets.male$Age_both, targets.male$DNAmAge, use="complete.obs")
plot(targets.male$Age_both, targets.male$DNAmAge, main = "Predicted age vs chronological age in Male", xlab = "Chronological age", ylab = "Predicted age (Horvath)")
abline(0,1,lwd=3)
abline(lm(targets.male$Age_both ~ targets.male$DNAmAge), col = 4 )

targets.female <- subset(targets, targets$Sex == "Female")
corrfemale <- cor(targets.female$Age_both, targets.female$DNAmAge, use="complete.obs")
plot(targets.female$Age_both, targets.female$DNAmAge, main = "Predicted age vs chronological age in Male", xlab = "Chronological age", ylab = "Predicted age (Horvath)")
abline(0,1,lwd=3)
abline(lm(targets.female$Age_both ~ targets.female$DNAmAge), col = 4 )

targets.mellick <- subset(targets, targets$SAMPLE.SET == "George_Mellick")
corrmell <- cor(targets.mellick$Age_both, targets.mellick$DNAmAge, use="complete.obs")
plot(targets.mellick$Age_both, targets.mellick$DNAmAge, main = "Predicted age vs chronological age in Mellick", xlab = "Chronological age", ylab = "Predicted age (Horvath)")
abline(0,1,lwd=3)
abline(lm(targets.mellick$Age_both ~ targets.mellick$DNAmAge), col = 4 )

targets.kennedy <- subset(targets, targets$SAMPLE.SET == "Martin_Kennedy")
corrken <- cor(targets.kennedy$Age_both, targets.kennedy$DNAmAge, use="complete.obs")
plot(targets.kennedy$Age_both, targets.kennedy$DNAmAge, main = "Predicted age vs chronological age in Kennedy", xlab = "Chronological age", ylab = "Predicted age (Horvath)")
abline(0,1,lwd=3)
abline(lm(targets.kennedy$Age_both ~ targets.kennedy$DNAmAge), col = 4 )
targets <- as.factor(targets$SAMPLE.SET)


# Barplot of regression
install.packages("ggthemes")
require(ggthemes)
regre <- data.frame(names = c("Control","PD","Male","Female","Mellick", "Kennedy", "Global") , corrval = c(corrcontrol, corrpd, corrmale, corrfemale, corrmell, corrken, corrcomplete ), nombr = c(1:7))
regre$names <- factor(regre$names, levels = regre$names)
plotcorr <- ggplot( regre, aes(x = names, y = corrval,  color= names, fill = names)) +
  geom_bar(stat='identity', position='dodge', width = .5) + 
  ggtitle("Correlations values between predicted age and chronological age") +  ylim(0,1) +
  xlab("Category") + ylab("Correlation") + scale_fill_hc() +scale_color_hc() +theme_hc()
plotcorr


# Plot of chronological age versus predicted age
plotage <- ggplot(targets, aes(x = Age_both, y = DNAmAge, color = Sex, shape = Phenotype) )+ geom_point() +
  ggtitle("Predicted age versus chronological age")+  geom_smooth(method = "lm",se=FALSE) +
  facet_grid(. ~ Phenotype) + geom_abline(intercept = 0, slope = 1, linetype="dotted", colour="black") +
  xlab("Chronological age") + ylab("Predicted age") + scale_fill_hc() +scale_color_hc() +theme_hc()
plotage 

# 12) Calculate the age acceleration
ageaccel <- targets$DNAmAge - targets$Age_both
mean(ageaccel, na.rm= TRUE)
plot(ageaccel)
plot(targets$Age_both, ageaccel, xlab = "Chronological age", ylab = "Predicted age", main = "Predicted age versus chronological age")

# Acceleration per cohort
tapply(ageaccel, targets$SAMPLE.SET, na.rm = TRUE, mean)
tapply(ageaccel, targets$Sex, na.rm = TRUE, mean)
tapply(ageaccel, targets$Phenotype, na.rm = TRUE, mean)

targetsage <- as.data.frame(cbind(targets, ageaccel)) 
plotaccel<- ggplot(targetsage, aes(x = Age_both, y = ageaccel, color = Phenotype)) + geom_point() +
  ggtitle("Age acceleration") +  geom_smooth(method = "lm",se=FALSE) +
  facet_grid(. ~ Phenotype) + geom_abline(intercept = 0, slope = 0, linetype="dotted", colour="black") +
  xlab("Chronological age") + ylab("Predicted age - Chronological age")
plotaccel


plotaccel2<- ggplot(targetsage, aes(x = Age_both, y = ageaccel, color = Sex)) + geom_point() +
  ggtitle("Age acceleration") +  geom_smooth(method = "lm",se=FALSE) +
  facet_grid(. ~ Sex) + geom_abline(intercept = 0, slope = 0, linetype="dotted", colour="black") +
  xlab("Chronological age") + ylab("Predicted age - Chronological age")
plotaccel2


# Boxplot on cell predictions
library(reshape2)
targets5 <- c( "Sex", "Phenotype", "CD8T","CD4T","NK", "Bcell", "Mono", "Gran")
targets5 <- targets[targets5]
datacell <- melt(targets5, id.var = c("Phenotype", "Sex"))
plotcell <- ggplot(datacell, aes(x = variable, y = value, color = Sex)) + geom_boxplot() +    facet_grid(. ~ Phenotype) +  xlab("Cell type") +
  ylab("Proportion")  + ggtitle("Cell type proportion per phenotype") 
plotcell

#Histogram
targetsfull <-read.csv("PD_targets_final_02June16_954_pred_complete.csv", fill = TRUE, header = T)
plotsmoking <- ggplot(targetsfull, aes( x = Smoking_prediction, fill=Sex))+ geom_histogram(colour = "black", binwidth = .02) + 
  geom_vline(aes(xintercept = mean(Smoking_prediction)), color="red", linetype="dashed", size=1) +
  scale_fill_hc() +scale_color_hc() +theme_hc() +
  labs(title="Smoking prediction per phenotype \n using methylation level for cg05575921") + labs(x="Predicted smoking status", y="Count")  + facet_grid( . ~ Phenotype)
  
plotsmoking

