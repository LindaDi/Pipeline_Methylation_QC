################################################################################
######################## Pipeline for QC methylation ###########################
################################################################################

################################
# 10/2020 
# Linda Dieckmann, Darina Czamara
# linda.dieckmann@psych.mpg.de
################################

#run everything on cluster due to memory
# please see Maksimovic et al. for reference of the main steps.
https://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html

###############################################
# load R packages
library(minfi) 
library(minfiData)
library(minfiDataEPIC)
library(RColorBrewer)	
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) 
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(sva)
library(RPMM)
library(ggplot2)
library(missMethyl)
library(matrixStats)
library(wateRmelon)
library(reshape)
library(Hmisc)

sessionInfo <- sessionInfo()

# set your working directory for the analysis
setwd("/folder/")

## should be folder with 
  # samplesheet
  # folder RData (empty)
  # folder Reports (empty)
  # folder addFiles 
    #-> for files with: sex of the samples, cross-reactive probes file from Chen et al., cross-hibridizing and polymorphic probes on EPIC from McCartney et al. 
  # BMIQ normalization function from Teschendorff
  # folder files_for_mixup (empty)
  # folder finalData (empty)

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

# please write down the samples and probe size for the different steps
# you can use an excel file sample_steps for documentation

#################################################################################
#### Read in iDATs and prepare data formats ####
#################################################################################

#######################################################
##### samplesheet:
# Sample_Name	Sample_Well	Sample_Plate	Slide	Array	Basename
  # pathfolder/slide/slide_array
# the column for the path must be named "Basename"
# most easiest to make in excel a new column "path" with the path to the files
# then make a new column named "Basename" and use a formula: path & slide & "/" & slide & "_" & array
# e.g. (=A2&B2&"/" ...)
# pull the entries down, so that this is filled in for all
#########################################################
##### Read in sheet with batch and ID info:
targets = read.csv("samplesheet_itu_lcs_for_minfi.csv")

RGSet	= read.metharray.exp(targets = targets)
dim(RGSet) # probe level with red/green intensities

save(RGSet, file = "RData/RGSet_original.Rdata")
######################################################### 
# load("RData/RGSet_original.Rdata")
##### Generate batch file
pd_orig <- pData(RGSet) # phenotype data
save(pd_orig, file="RData/pd_original.Rdata")
#########################################################
##### Get raw beta values:
Mset = preprocessRaw(RGSet) # cpG locus level, with 2 channels methylated/ unmethylated
save(Mset,file="RData/Mset_original.Rdata")

RatioSet = ratioConvert(Mset, what = "both", keepCN = TRUE)# CpG locus level, but not mapped to a genome, Beta and M values
save(RatioSet, file="RData/RatioSet_original.Rdata") 

RawBetas = getBeta(RatioSet) # Beta value matrix
save(RawBetas,file="RData/RawBetas_original.Rdata") 

gRatioSet<-mapToGenome(RatioSet,mergeManifest=TRUE) 
save(gRatioSet, file = "RData/gRatioSet_original.Rdata") 
########################################################

#################################################################################
#### Quality checks and exclusion of samples ####
#################################################################################
##### detection p values
# calculate the detection p-values 
detP <- detectionP(RGSet) 
head(detP) 
save(detP, file = "RData/detP.Rdata") 

# examine mean detection p-values across all samples 
pdf("Reports/detectionP.pdf", width = 8, height = 3)
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2)) 
barplot(colMeans(detP), col=pal[factor(targets$Sample_Name)], las=2, cex.names=0.4,ylab="Mean detection p-values") 
abline(h=0.01,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Name)), fill=pal, bg="white") 
barplot(colMeans(detP), col=pal[factor(targets$Sample_Name)], las=2, cex.names=0.4, ylim = c(0,0.002), ylab="Mean detection p-values") 
legend("topleft", legend=levels(factor(targets$Sample_Name)), fill=pal, bg="white") 
dev.off()

# remove poor quality samples with a mean detection p-value >.05
keep <- colMeans(detP) < 0.05 
RGSet_qual <- RGSet[,keep]
RGSet_qual
save(RGSet_qual, file = "RData/RGSet_qual.Rdata") 

# remove poor quality samples from raw betas
load("RData/RawBetas_original.Rdata")
RawBetas_qual <- RawBetas[,keep]
dim(RawBetas_qual)
save(RawBetas_qual, file = "RData/RawBetas_qual.Rdata") 

# remove poor quality samples from targets data = targets_qual
targets_qual <- targets[keep,] 
targets_qual[,1:5] 

# remove poor quality samples from detection p-value table 
detP_qual <- detP[,keep] 
dim(detP_qual) 
save(detP_qual, file = "RData/detP_qual.Rdata") 

# -> note how many samples were excluded

# run minfi QC Report
qcReport(RGSet_qual, sampGroups=targets_qual$Slides, pdf="Reports/qcReport.pdf") 
#1) density plots of the methylation Beta values for all samples, typically colored by sample group. 
  #While the density plots are useful for identifying deviant samples, it is not easy to identify the specific problem sample
#2) bean plots show each sample in its own section. 
  #While the shape of the distribution for “good” samples will differ from experiment to experiment, many conditions have methylation profiles characterized by two modes - one with close to 0% methylation, and a second at close to 100% methylation.
#3) controlStripPlot function allows plotting of individual control probe types

###############################################################################
##### distribution artefacts
# this is added to pipeline based on practical experience

#check individual densities
pdf("Reports/beta_densities.pdf")
for (i in 1:ncol(RGSet_qual))
{titel<-paste(rownames(pData(RGSet_qual))[i])
densityPlot(as.matrix(RawBetas_qual[,i]),main=titel)
print(i)}
dev.off() 
# -> note slides with different distributions
# -> note IDs to exclude
################################################################################

################################################################################
#check sex matches - added 
################################################################################
# note that this is done on whole sample, first ignoring that we already exclude based on detect p & distribution
load("RData/gRatioSet_original.Rdata")

# predict sex with methylation data
predictedSex<-getSex(gRatioSet,cutoff=-2)
sex<-cbind(sampleNames(RGSet),as.character(pd_orig$Sample_Name),predictedSex$predictedSex,predictedSex$yMed,predictedSex$xMed,predictedSex$yMed-predictedSex$xMed)
sex<-as.data.frame(sex)
names(sex)<-c("ArrayID","ID","predictedSex","yMed","xMed","yMed-xMed")
sex$yMed<-as.numeric(as.character(sex$yMed))
sex$xMed<-as.numeric(as.character(sex$xMed))
sex$ID <- as.character(sex$ID)
save(sex, file="RData/sex_predicted.Rdata")

###
# if necessary, preprocess txt file with gender info in RStudio (depends on the kind of file you have)
gender<-read.table("addFiles/pheno_gender_corrected.txt",sep="\t",header=T, stringsAsFactors=F)
gender <- gender[ , c("SampleName","Gender")]
gender$SampleName <- gsub("-", "_", gender$SampleName)
gender$Gender <- gsub("W", "F", gender$Gender)
gender$SampleName <- gsub("i", "I", gender$SampleName)
save(gender, file="RData/gender.RData")
write.table(gender,"addFiles/gender.txt",sep="\t",row.names=FALSE)
###

load("RData/sex_predicted.Rdata")
load("RData/gender.RData")
test<-merge(sex,gender,by.x="ID",by.y="SampleName")
NoIDMatchSexGender <- subset(sex$ID, !(sex$ID %in% gender$SampleName))
NoIDMatchSexGender
# -> check if empty, there should be no ID mismatches

write.table(test,"Reports/predicted_vs_childsex.txt",sep="\t",quote=F,row.names=F)
table(test$Gender,test$predictedSex)

test[test$Gender=="M" & test$predictedSex=="F",]
# -> note IDs and ArrayIDs

test[test$Gender=="F" & test$predictedSex=="M",]
# -> note IDs and ArrayIDs

# recheck if gender is correct
# which are true mismatches
# -> exclude these before ComBat
##########################################################
# sample exclusions:
  # distribution artefacts
  # poor quality (detection p-value)
  # sex mismatches
####### Individuals to exclude (BEFORE NORMALIZATION) :
# we already excluded samples based on detection p-value
# we need to further exclude the distribution artefacts and sex mismatches from the files

# -> write them in out_index (ArrayIDs):
out_index<- c("")
duplicated(out_index) # see if an ID was found in both steps
i1 <- intersect(out_index, colnames(RGSet_qual)) # see how many are still in RGSet_qual
i2 <- intersect(out_index, colnames(RGSet)) # re-check that all identified IDs are in the full RGset (no misspelling)
setdiff(i2, i1) # ID(s) identified in detP step and again in later quality control

RGSet_clean <- RGSet_qual[ , !(colnames(RGSet_qual) %in% out_index)]
RGSet_clean
save(RGSet_clean, file = "RData/RGSet_clean.Rdata")

RawBetas_clean = RawBetas_qual[ ,!(colnames(RawBetas_qual) %in% out_index)]
save(RawBetas_clean, file = "RData/RawBetas_clean.Rdata") 

detP_clean <- detP_qual[,!(colnames(detP_qual) %in% out_index)] 
save(detP_clean, file = "RData/detP_clean.Rdata") 

pd_clean <-pData(RGSet_clean)
save(pd_clean, file = "RData/pd_clean.Rdata") 

## Get annotations probes:
annot = getAnnotation(RGSet_clean)
save(annot,file="RData/annot_clean.Rdata")

#################################################################################
### Normalization       (This step takes a while and requires lots of memory):
##################################################################################
# normalization is done before probe filtering, see answer to reviewer at the end of Maskimovic pipeline
# there are a lot of different normalization methods, the best depends also on study design and question
# noob is for background correction, and often a good first step, see e.g. https://www.biostars.org/p/149628/)
# Quantile normalization, is used in  the pipeline 
  # and also recommended for studies without big differences in samples (e.g. case-control) from Fortin
  # IF you have quite different samples (e.g. different tissues, experimental design etc) it is better to use FunNorm from minfi
# BMIQ is also very common and was suggested in combination with other methods (e.g. noob: Liu et al. 2016, or QN (from lumi package):Marabita et al. 2013, Wang et al. 2015 -> but from R pckg lumi)
  # note regarding BMIQ: BMIQ from wateRmelon might be not stable, especially after minfi quantile -> used an original version from Andrew Teschendorff
source("BMIQ_1.6_Teschendorff.R")

load("RData/RGSet_clean.Rdata")
load("RData/detP_clean.Rdata")
load("RData/annot_clean.Rdata")
load("RData/RawBetas_clean.Rdata")
load("RData/pd_clean.Rdata")

# 1. if you want noob normalization
  # noob = preprocessNoob(RGSet_clean)
  # save(noob, file="RData/noob.Rdata")
# -> performs a background correction based on the negative control probes (Triche, 2013)
# -> output methylSet

# 1.2 preprocessQuantile after noob (e.g. recommended from Fortin et al.):
  # GSetquantile.norm = preprocessQuantile(noob)
  # save(GSetquantile.norm, file="RData/GSetquantile.norm")
  # -> output GenomicRatioSet
  # GSetquantileBetas = getBeta(GSetquantile.norm) # get betas

# 1.3 BMIQ after noob:
  # BMIQ.norm = BMIQ(noob)
  # -> output is a matrix with beta values
  # save(BMIQ.norm, file="RData/BMIQ.norm.Rdata")

# 2. if you need functional normalization
  # Funnorm = preprocessFunnorm(RGSet_clean, sex = NULL, bgCorr = T, dyeCorr = T, nPCs = 2, verbose = TRUE)
# -> output GenomicRatioSet
  # save(Funnorm, file="RData/Funnorm.Rdata")
  # Funnorm_Betas = getBeta(Funnorm)
  # save(Funnorm_Betas, file="RData/Funnorm_Betas.Rdata")

# -> based on pipeline, literature and data inspection we decided to implement
#  stratified quantile normalization (followed by BMIQ) as a standard:

# Quantile normalization (our standard):
quantileN = preprocessQuantile(RGSet_clean)
save(quantileN, file="RData/quantileN.Rdata")
# -> output: GenomicRatioSet
quantileNBetas = getBeta(quantileN) # get betas
save(quantileNBetas,file="RData/quantileNBetas.Rdata")
quantileNMs = getM(quantileN)  
save(quantileNMs,file="RData/quantileNMs.Rdata")

  # -> stratified quantile normalization 
  # The distributions of probe intensities for different samples are made identical. Often used in microarray analysis.
  # Stratified quantile normalisation: Probes are stratified by genomic region then quantile normalised (The distributions of probe intensities for different samples are made identical.)
    # with sex chromosomes normalised separately when male and female samples are present. 
    # No background correction, zeros removed by outlier function. Not recommended for cancer-normal comparisons or other groups with global differences.

# BMIQ after quantile normalization:
probeType = as.data.frame(annot[rownames(quantileNBetas),c("Name","Type")])
probeType$probeType = ifelse(probeType$Type %in% "I",1,2)

BMIQ.quantileN = apply(quantileNBetas[,1:length(colnames(quantileNBetas))],2,function(a) BMIQ(a,probeType$probeType,plots=FALSE)$nbeta)

length(which(is.nan(BMIQ.quantileN))) # should be 0
save(BMIQ.quantileN, file="RData/BMIQ.quantileN.Rdata")

#
# #Check distributions before and after normalization
# png(file="Reports/Beta_Distributions_Norm.png")
# par(mfrow=c(3,2))
# densityPlot(RawBetas_clean, sampGroups = pd_clean$Slide, legend=FALSE, main = "Raw Betas", xlab = "Beta")
# densityPlot(quantileNBetas, sampGroups = pd_clean$Slide, legend=FALSE, main = "QuanNorm adjusted Betas", xlab = "Beta")
# densityPlot(GSetquantileBetas, sampGroups = pd_clean$Slide, legend=FALSE, main = "noobQuanNorm adjusted Betas", xlab = "Beta")
# densityPlot(BMIQ.norm, sampGroups = pd_clean$Slide, legend=FALSE, main = "noobBMIQ adjusted Betas", xlab = "Beta")
# densityPlot(BMIQ.quantileN, sampGroups = pd_clean$Slide, legend=FALSE, main = "quantileBMIQ adjusted Betas", xlab = "Beta")
# densityPlot(Funnorm_Betas, sampGroups = pd_clean$Slide, legend=FALSE, main = "Funnorm adjusted Betas", xlab = "Beta")
# dev.off()
# 

## Check distributions before and after chosen normalization
png(file="Reports/BetaValue_Distributions_Norm_Quantile.png",width=1400,height=700,pointsize=12)
par(mfrow=c(1,2))
densityPlot(RawBetas_clean, sampGroups = pd_clean$Slide, legend=FALSE, main = "Raw Betas", xlab = "Beta")
densityPlot(quantileNBetas, sampGroups = pd_clean$Slide, legend=FALSE, main = "Quantile Adjusted Betas", xlab = "Beta")
densityPlot(BMIQ.quantileN, sampGroups = pd_clean$Slide, legend=FALSE, main = "Quantile-BMIQ Adjusted Betas", xlab = "Beta", ylim=c(0,5))
dev.off()

# 
# ### to see potential outliers in the distribution
# # I) raw batas
# names_RawBetas <- colnames(RawBetas_clean)
# pdf("Reports/Beta_Densities_RawBetas.pdf")
# for (i in 1:ncol(RawBetas_clean)) {
#   i_mat <- as.matrix(RawBetas_clean[ ,i])
#   densityPlot(i_mat, main=names_RawBetas[i])
#   }
# dev.off()
# # II) Quantile normalized betas
# names_quantileBetas <- colnames(quantileNBetas)
# pdf("Reports/Beta_Densities_QuantileBetas.pdf")
# for (i in 1:ncol(quantileNBetas)) {
#   i_mat <- as.matrix(quantileNBetas[ ,i])
#   name <- colnames(quantileNBetas[,i])
#   densityPlot(i_mat, main=names_quantileBetas[i])
# }
# dev.off()

############################################################################################
### Filtering probes
############################################################################################
# NOTE: we will filter both in the BMIQ beta object and quantileN Genomic ratio object
# as GenomicRatio Sets contain more info and can be also used for M values
# further, this makes it easier to switch to only stratified quantile normalization, if wanted or needed.

##### filter out probes that failed in one or more samples based on detection p < .01
# ensure probes are in the same order in the BMIQ beta and detP objects
detP_clean_f <- detP_clean[match(rownames(BMIQ.quantileN),rownames(detP_clean)),]
# ensure probes are in the same order in the quantileN (GenomicRatio) and detP objects
detP_clean_ff <- detP_clean[match(featureNames(quantileN),rownames(detP_clean)),]

# remove any probes that have failed in one or more samples
keep <- rowSums(detP_clean_f < 0.01) == ncol(BMIQ.quantileN) 
table(keep) # -> note output (excel)
keep_ff <- rowSums(detP_clean_ff < 0.01) == ncol(quantileN) 
table(keep_ff) 
#  FALSE   TRUE 

# exclude these probes from our normalized data
BMIQ.quantileN_filtered <- BMIQ.quantileN[keep,]
# matrix
dim(BMIQ.quantileN_filtered)
# 
quantileN_filtered <- quantileN[keep_ff,]
# genomicRatioSet
quantileN_filtered

##### filter out probes on sex chromosomes
keep <- !(rownames(BMIQ.quantileN_filtered) %in% annot$Name[annot$chr %in% c("chrX","chrY")])
table(keep) # -> note output (excel)
BMIQ.quantileN_filtered <- BMIQ.quantileN_filtered[keep,]
dim(BMIQ.quantileN_filtered)
# 
keep_ff <- !(featureNames(quantileN_filtered) %in% annot$Name[annot$chr %in% c("chrX","chrY")])
table(keep_ff) 
# 16916 795677
quantileN_filtered <- quantileN_filtered[keep_ff,]
quantileN_filtered

##### removal of probes where common SNPs may affect the CpG. 
# You can either remove all probes affected by SNPs (default), or only those with minor allele frequencies greater than a specified value.
quantileN_filtered <- dropLociWithSnps(quantileN_filtered) # this NEEEDs a GenomicRatioSet, beta matrix does not work
quantileN_filtered # -> note again dims
# 
BMIQ.quantileN_filtered <-  BMIQ.quantileN_filtered[rownames(BMIQ.quantileN_filtered) %in% featureNames(quantileN_filtered),]
dim(BMIQ.quantileN_filtered)
# 

# load Chen Probe annotation file and exclude SNP and Cross Hybridizing probes:
load("addFiles/ChenProbeIDs.rdata") #### Data from Chen et al (2013) PMID: 23314698
# loaded as annot2
annot2$SNPs = annot2[,"EUR"] #### Can change Global to "AFR", "AMR", "ASN", or "EUR" to match content specific allelic frequency; Use "Global" if population is very diverse
index<-which(annot2$sex=="Exclude" | annot2$CRSH=="Exclude" | annot2$EUR=="Exclude")
length(index) #62520
exclude_Chen <-annot2[index,]

#not in pipeline, but additionally use data from McCartney et al which are based on EPIC# PMID: 27330998
exclude_crosshyb <-read.table("addFiles/CpGs_crosshybridizing_EPIC.txt",sep="",header=F)
keep <- !(rownames(BMIQ.quantileN_filtered) %in% exclude_crosshyb$V1) 
table(keep) # -> note output
BMIQ.quantileN_filtered <- BMIQ.quantileN_filtered[keep,] 
dim(BMIQ.quantileN_filtered)
#
keep_ff <- !(featureNames(quantileN_filtered) %in% exclude_crosshyb$V1) 
table(keep_ff) # 
quantileN_filtered <- quantileN_filtered[keep_ff,] 
quantileN_filtered # 

exclude_poly <-read.table("addFiles/CpGs_polymorphic_EPIC.txt",sep="",header=T)
index<-which(exclude_poly$EUR_AF>0.05) #n=10971
exclude_polym <- exclude_poly[index,]

keep <- !(rownames(BMIQ.quantileN_filtered) %in% exclude_polym$IlmnID) 
table(keep) # -> note output
BMIQ.quantileN_filtered <- BMIQ.quantileN_filtered[keep,] 
dim(BMIQ.quantileN_filtered)
keep_ff <- !(featureNames(quantileN_filtered) %in% exclude_polym$IlmnID) 
table(keep_ff) 
quantileN_filtered <- quantileN_filtered[keep_ff,] 

### ### ###
save(BMIQ.quantileN_filtered, file="RData/BMIQ.quantileN_filtered.Rdata") # QN+BMIQ normalized, probe filtered & sample-qced Betas.
save(quantileN_filtered, file="RData/quantileN_filtered.Rdata") # GenomicRatio of filtered QN normalized data
Betas_quantileN_filtered <- getBeta(quantileN_filtered) # beta values
save(Betas_quantileN_filtered, file="RData/Betas_quantileN_filtered.Rdata")
Ms_quantileN_filtered <- getM(quantileN_filtered)
save(Ms_quantileN_filtered, file="RData/Ms_quantileN_filtered.Rdata")

## Check density plots after excluding the poorly-detected probes:
png(file="Reports/BetaValue_Distributions_Norm_quantile_Filter.png")
densityPlot(BMIQ.quantileN_filtered,sampGroups = pd_clean$Slide, legend=FALSE, main = "Post Filtering - Normalized Beta", xlab = "Beta")
dev.off() 


################################################################################
# control for batch effects with combat          
################################################################################
##### This is added to the pipeline
# We use combat to remove batch effects
# see e.g. Wen Bin Goh et al. 2017 for importance of batch effects removement 

#combat to remove batch effects
## function that will calculate the variance of each row
rowVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE, SumSquares=FALSE, twopass=FALSE) {
  if (SumSquares) return(rowSums(x^2, na.rm, dims))
  N <- rowSums(!is.na(x), FALSE, dims)
  Nm1 <- if (unbiased) N-1 else N
  if (twopass) {x <- if (dims==0) x - mean(x, na.rm=na.rm) else
    sweep(x, 1:dims, rowMeans(x,na.rm,dims))}
  (rowSums(x^2, na.rm, dims) - rowSums(x, na.rm, dims)^2/N) / Nm1
}

mval <- apply(BMIQ.quantileN_filtered, 2, function(x) log2((x)/(1-x))) # M values

## Calculate the variance of each probe and remove any with a variance of 0 prior to Combat.
vars = as.matrix(rowVars(mval))
which(vars==0) # 0

## Replace all probes with no variance with NA and remove them from the normalized data set
vars[vars == 0] = NA
vars = na.omit(vars)
intersect = intersect(rownames(vars), rownames(mval))
print(length(intersect)) # probes without variance == 0

BMIQ.quantileN_filtered_batch = BMIQ.quantileN_filtered[intersect,]
mval = mval[intersect,]

## Ensure Objects are aligned
table(ifelse(rownames(pd_clean) == colnames(mval),"Match","Off")) # All should match

## Check variation in array data associated with batch (ie. Slide/plate/box)
## Run a principle component analysis to determine if there are any remaining batch effects following data normalization.

PCobj = prcomp(t(mval), retx = T, center = T, scale. = T)
save(PCobj, file="RData/PCobj.Rdata")

pdf("Reports/boxplot_PCA.pdf")
boxplot(PCobj$x,col="grey",frame=F)
dev.off()

# Can use Scree plot to determine number of PCs to keep
pdf("Reports/screeplot_PCA.pdf")
plot(PCobj,type="line",cex.lab=1.5, cex.main=1.5) 
dev.off()

# Extract the PCs from the PCobj object
PCs = PCobj$x

# Extract the proportion of variability and cumulative proportion of 
# varibility explained by the top R PCs.
R = 5
propvar = summary(PCobj)$importance["Proportion of Variance", 1:R] # -> write in xlsx
cummvar = summary(PCobj)$importance["Cumulative Proportion", 1:R] # -> write in xlsx

R = 4 # choose eg 4 PCs
## Generate plots of the resulting PCs
# Plot of the proportion of variability explained by the top R PCs
# Plot of the cummulative proportion of variability explained by the top R PCs
pdf("Reports/PCA_variance_explained.pdf")
par(mfrow=c(1,2))	
par(mar = c(5,5,4,2))
barplot(propvar*100, xlab = paste("Top", R, "PCs", sep = " "), ylab = "Variation Explained (%)", cex.axis = 1.5, cex.lab = 1.8, cex.names = 1.5)
par(mar = c(5,5,4,2))
barplot(cummvar*100, xlab = paste("Top", R, "PCs", sep = " "), ylab = "Cumulative Variation Explained (%)",cex.axis = 1.5, cex.lab = 1.8, cex.names = 1.5)
dev.off()

# Plot of PCX and PCY; by Batch     
PCs = PCobj$x
PCs =PCs[,1:R]
Prin.comp<-merge(PCs,pd_clean, by = "row.names",all=T) 

pdf(file="Reports/PC_Variation_by_batch.pdf")
par(mfrow=c(3,1))
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Sample_Plate), xlab = "PC1", ylab = "PC2", main="Plate")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Slide), xlab = "PC1", ylab = "PC2", main="Slide")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Array), xlab = "PC1", ylab = "PC2", main="Array")
dev.off()

# extreme outliers?
o1 <- 3*sd(Prin.comp$PC1)
o2 <- 3*sd(Prin.comp$PC2)
which(abs(Prin.comp$PC1) > o1 && abs(Prin.comp$PC2) > o2)


## Can  test via ANOVA, to look for strong variations in PCs by Batch:
anova(lm(Prin.comp$PC1~Prin.comp$Sample_Plate)) 
anova(lm(Prin.comp$PC2~Prin.comp$Sample_Plate)) 
anova(lm(Prin.comp$PC3~Prin.comp$Sample_Plate)) 
anova(lm(Prin.comp$PC4~Prin.comp$Sample_Plate)) 

#for slide
anova(lm(Prin.comp$PC1~Prin.comp$Slide)) 
anova(lm(Prin.comp$PC2~Prin.comp$Slide)) 
anova(lm(Prin.comp$PC3~Prin.comp$Slide))
anova(lm(Prin.comp$PC4~Prin.comp$Slide)) 


#for Array
anova(lm(Prin.comp$PC1~Prin.comp$Array)) 
anova(lm(Prin.comp$PC2~Prin.comp$Array))
anova(lm(Prin.comp$PC3~Prin.comp$Array)) 
anova(lm(Prin.comp$PC4~Prin.comp$Array)) 

#first correct for plate
## Model matrix for batch-corrections (May need to adjust model matrix to 'protect' coefficients (study specific)):
mod <- model.matrix(~1, data=pd_clean)

## Run ComBat to remove most significant batch effects itertively
M_combat_1plate = ComBat(mval,batch = pd_clean$Sample_Plate, mod = mod)
save(M_combat_1plate,file="RData/M_combat_1plate.Rdata")

## Check to see if batch effect was succesfully removed
PCobj = prcomp(t(M_combat_1plate), retx = T, center = T, scale. = T)
PCs = PCobj$x
PCs =PCs[,1:R]
Prin.comp<-merge(PCs,pd_clean, by = "row.names",all=T) 

#### Check whether batches are still distinguished by first and second PC:
pdf(file="Reports/PC_Variation_by_batch_afterCombat1.pdf")
par(mfrow=c(3,1))
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Sample_Plate), xlab = "PC1", ylab = "PC2", main="Plate")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Slide), xlab = "PC1", ylab = "PC2", main="Slide")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Array), xlab = "PC1", ylab = "PC2", main="Array")
dev.off()

## Can further test via ANOVA, to look for strong variations in PCs by Batch:
anova(lm(Prin.comp$PC1~Prin.comp$Sample_Plate)) 
anova(lm(Prin.comp$PC2~Prin.comp$Sample_Plate)) 
anova(lm(Prin.comp$PC3~Prin.comp$Sample_Plate)) 
anova(lm(Prin.comp$PC4~Prin.comp$Sample_Plate)) 

#for slide
anova(lm(Prin.comp$PC1~Prin.comp$Slide)) 
anova(lm(Prin.comp$PC2~Prin.comp$Slide)) 
anova(lm(Prin.comp$PC3~Prin.comp$Slide)) 
anova(lm(Prin.comp$PC4~Prin.comp$Slide)) 

#for array
anova(lm(Prin.comp$PC1~Prin.comp$Array))
anova(lm(Prin.comp$PC2~Prin.comp$Array)) 
anova(lm(Prin.comp$PC3~Prin.comp$Array)) 
anova(lm(Prin.comp$PC4~Prin.comp$Array)) 

#second correction for slide
mod <- model.matrix(~1, data=pd_clean)
M_combat_2slide = ComBat(M_combat_1plate,batch = pd_clean$Slide, mod = mod)
save(M_combat_2slide,file="RData/M_combat_2slide.Rdata")       

## Check to see if batch effect was succesfully removed
PCobj = prcomp(t(M_combat_2slide), retx = T, center = T, scale. = T)
PCs = PCobj$x
PCs =PCs[,1:R]
Prin.comp<-merge(PCs,pd_clean, by = "row.names",all=T) 

#### Check whether batches are still distinguished by first and second PC:
pdf(file="Reports/PC_Variation_by_batch_afterCombat2.pdf")
par(mfrow=c(3,1))
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Sample_Plate), xlab = "PC1", ylab = "PC2", main="Plate")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Slide), xlab = "PC1", ylab = "PC2", main="Slide")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Array), xlab = "PC1", ylab = "PC2", main="Array")
dev.off()

## Can further test via ANOVA, to look for strong variations in PCs by Batch:
# for plate
anova(lm(Prin.comp$PC1~Prin.comp$Sample_Plate)) 
anova(lm(Prin.comp$PC2~Prin.comp$Sample_Plate)) 
anova(lm(Prin.comp$PC3~Prin.comp$Sample_Plate)) 
anova(lm(Prin.comp$PC4~Prin.comp$Sample_Plate)) 

#for slide
anova(lm(Prin.comp$PC1~Prin.comp$Slide))
anova(lm(Prin.comp$PC2~Prin.comp$Slide)) 
anova(lm(Prin.comp$PC3~Prin.comp$Slide)) 
anova(lm(Prin.comp$PC4~Prin.comp$Slide)) 

#for array
anova(lm(Prin.comp$PC1~Prin.comp$Array)) 
anova(lm(Prin.comp$PC2~Prin.comp$Array)) 
anova(lm(Prin.comp$PC3~Prin.comp$Array)) 
anova(lm(Prin.comp$PC4~Prin.comp$Array)) 

#third correction for array
mod <- model.matrix(~1, data=pd_clean)
M_combat_3array = ComBat(M_combat_2slide,batch = pd_clean$Array, mod = mod)
save(M_combat_3array,file="RData/M_combat_3array.Rdata")       

## Check to see if batch effect was succesfully removed
PCobj = prcomp(t(M_combat_3array), retx = T, center = T, scale. = T)
PCs = PCobj$x
PCs =PCs[,1:R]
Prin.comp<-merge(PCs,pd_clean, by = "row.names",all=T) 

#### Check whether batches are still distinguished by first and second PC:
pdf(file="Reports/PC_Variation_by_batch_afterCombat3.pdf")
par(mfrow=c(3,1))
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Sample_Plate), xlab = "PC1", ylab = "PC2", main="Plate")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Slide), xlab = "PC1", ylab = "PC2", main="Slide")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Array), xlab = "PC1", ylab = "PC2", main="Array")
dev.off()

## Can further test via ANOVA, to look for strong variations in PCs by Batch:

# for plate
anova(lm(Prin.comp$PC1~Prin.comp$Sample_Plate)) 
anova(lm(Prin.comp$PC2~Prin.comp$Sample_Plate)) 
anova(lm(Prin.comp$PC3~Prin.comp$Sample_Plate)) 
anova(lm(Prin.comp$PC4~Prin.comp$Sample_Plate)) 

#for slide
anova(lm(Prin.comp$PC1~Prin.comp$Slide))
anova(lm(Prin.comp$PC2~Prin.comp$Slide)) 
anova(lm(Prin.comp$PC3~Prin.comp$Slide)) 
anova(lm(Prin.comp$PC4~Prin.comp$Slide)) 

#for array
anova(lm(Prin.comp$PC1~Prin.comp$Array)) 
anova(lm(Prin.comp$PC2~Prin.comp$Array)) 
anova(lm(Prin.comp$PC3~Prin.comp$Array)) 
anova(lm(Prin.comp$PC4~Prin.comp$Array)) 


#-> if nothing significant ok, leave it like that.


#################################################################################
## Convert the batch-adjusted M-values back into betas:

expit2 = function(x) 2^x/(1+2^x)
Betas_combated = expit2(M_combat_3array)
dim(Betas_combated) # for now final betas!
# 716331    264

## Save normalized and batch-adjusted beta values
save(Betas_combated,file = "RData/Betas_combated.Rdata")

#plot final densities
png(file="Reports/BetaValue_Distributions_afterNormCombat.png",width=700,height=700,pointsize=12)
densityPlot(Betas_combated, sampGroups = pd_clean$Slide, legend=FALSE, main = "PostQC - Normalized and Batch Corrected Beta", xlab = "Beta")
dev.off() 

all.equal(colnames(Betas_combated),rownames(pd_clean)) #TRUE
annotated_pd_clean =new("AnnotatedDataFrame", data= as.data.frame(pd_clean)) #extend to AnnotatedDataFrame (required for ESet)

Betas_combated_ExprSet = new("ExpressionSet", exprs= as.matrix(Betas_combated), phenoData=annotated_pd_clean)
save(Betas_combated_ExprSet,file="RData/Betas_combated_ExprSet.Rdata")

#################################################################################
#################################################################################
### We also combat our data set with no probes removed, 
### in case all CpGs are needed (e.g. for epigenetic clocks), i.e. when you need to avoid missing values (CpGs) even at costs.

load("RData/BMIQ.quantileN.Rdata")

#combat to remove batch effects
## function that will calculate the variance of each row
rowVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE, SumSquares=FALSE, twopass=FALSE) {
  if (SumSquares) return(rowSums(x^2, na.rm, dims))
  N <- rowSums(!is.na(x), FALSE, dims)
  Nm1 <- if (unbiased) N-1 else N
  if (twopass) {x <- if (dims==0) x - mean(x, na.rm=na.rm) else
    sweep(x, 1:dims, rowMeans(x,na.rm,dims))}
  (rowSums(x^2, na.rm, dims) - rowSums(x, na.rm, dims)^2/N) / Nm1
}

# ComBat on BMIQ.quantileN (no probes were removed there)
mval2 <- apply(BMIQ.quantileN, 2, function(x) log2((x)/(1-x))) # M values

## Calculate the variance of each probe and remove any with a variance of 0 prior to Combat.
vars2 = as.matrix(rowVars(mval2))
which(vars2==0) # 0

## Replace all probes with no variance with NA and remove them from the normalized data set
vars2[vars2 == 0] = NA
vars2 = na.omit(vars2)
intersect2 = intersect(rownames(vars2), rownames(mval2))
print(length(intersect2)) # probes without variance == 0

BMIQ.quantileN_batch = BMIQ.quantileN[intersect2,]
mval2 = mval2[intersect2,]

## Ensure Objects are aligned
table(ifelse(rownames(pd_cleaned) == colnames(mval2),"Match","Off")) # All should match

## Check variation in array data associated with batch (ie. Slide/plate/box)
## Run a principle component analysis to determine if there are any remaining batch effects following data normalization.

PCobj2 = prcomp(t(mval2), retx = T, center = T, scale. = T)
save(PCobj2, file="RData/PCobj2.Rdata")

pdf("Reports/boxplot_PCA_2.pdf")
boxplot(PCobj2$x,col="grey",frame=F)
dev.off()

# Can use Scree plot to determine number of PCs to keep
pdf("Reports/screeplot_PCA_2.pdf")
plot(PCobj2,type="line",cex.lab=1.5, cex.main=1.5) 
dev.off()

# Extract the PCs from the PCobj object
PCs2 = PCobj2$x

# Extract the proportion of variability and cumulative proportion of 
# varibility explained by the top R PCs.
R = 5
propvar = summary(PCobj2)$importance["Proportion of Variance", 1:R] 
cummvar = summary(PCobj2)$importance["Cumulative Proportion", 1:R] # 

R = 5 # choose PCs
## Generate plots of the resulting PCs
# Plot of the proportion of variability explained by the top R PCs
# Plot of the cummulative proportion of variability explained by the top R PCs
pdf("Reports/PCA_variance_explained_2.pdf")
par(mfrow=c(1,2))	
par(mar = c(5,5,4,2))
barplot(propvar*100, xlab = paste("Top", R, "PCs", sep = " "), ylab = "Variation Explained (%)", cex.axis = 1.5, cex.lab = 1.8, cex.names = 1.5)
par(mar = c(5,5,4,2))
barplot(cummvar*100, xlab = paste("Top", R, "PCs", sep = " "), ylab = "Cumulative Variation Explained (%)",cex.axis = 1.5, cex.lab = 1.8, cex.names = 1.5)
dev.off()

# Plot of PCX and PCY; by Batch     
PCs2 = PCobj2$x
PCs2 =PCs2[,1:R]
Prin.comp2 <-merge(PCs2,pd_clean, by = "row.names",all=T) 

pdf(file="Reports/PC_Variation_by_batch_2.pdf")
par(mfrow=c(3,1))
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Sample_Plate), xlab = "PC1", ylab = "PC2", main="Plate")
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Slide), xlab = "PC1", ylab = "PC2", main="Slide")
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Array), xlab = "PC1", ylab = "PC2", main="Array")
dev.off()

# extreme outliers?
o1 <- 3*sd(Prin.comp2$PC1)
o2 <- 3*sd(Prin.comp2$PC2)
which(abs(Prin.comp2$PC1) > o1 && abs(Prin.comp2$PC2) > o2)

## Can  test via ANOVA, to look for strong variations in PCs by Batch:
anova(lm(Prin.comp2$PC1~Prin.comp2$Sample_Plate)) 
anova(lm(Prin.comp2$PC2~Prin.comp2$Sample_Plate)) 
anova(lm(Prin.comp2$PC3~Prin.comp2$Sample_Plate)) 
anova(lm(Prin.comp2$PC4~Prin.comp2$Sample_Plate)) 
anova(lm(Prin.comp2$PC5~Prin.comp2$Sample_Plate)) 

#for slide
anova(lm(Prin.comp2$PC1~Prin.comp2$Slide)) 
anova(lm(Prin.comp2$PC2~Prin.comp2$Slide)) 
anova(lm(Prin.comp2$PC3~Prin.comp2$Slide))
anova(lm(Prin.comp2$PC4~Prin.comp2$Slide)) 
anova(lm(Prin.comp2$PC5~Prin.comp2$Slide)) 

#for Array
anova(lm(Prin.comp2$PC1~Prin.comp2$Array)) 
anova(lm(Prin.comp2$PC2~Prin.comp2$Array))
anova(lm(Prin.comp2$PC3~Prin.comp2$Array)) 
anova(lm(Prin.comp2$PC4~Prin.comp2$Array)) 
anova(lm(Prin.comp2$PC5~Prin.comp2$Array))

#first correct for plate
## Model matrix for batch-corrections (May need to adjust model matrix to 'protect' coefficients (study specific)):
mod2 <- model.matrix(~1, data=pd_cleaned)

## Run ComBat to remove plate batch effects
M_combat2_1plate = ComBat(mval2,batch = pd_cleaned$Sample_Plate, mod = mod2)
save(M_combat2_1plate,file="RData/M_combat2_1plate.Rdata")

## Check to see if batch effect was succesfully removed
PCobj2 = prcomp(t(M_combat2_1plate), retx = T, center = T, scale. = T)
PCs2 = PCobj2$x
PCs2 =PCs2[,1:R]
Prin.comp2 <-merge(PCs2,pd_cleaned, by = "row.names",all=T) 

#### Check whether batches are still distinguished by first and second PC:
pdf(file="Reports/PC_Variation_by_batch_2afterCombat1.pdf")
par(mfrow=c(3,1))
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Sample_Plate), xlab = "PC1", ylab = "PC2", main="Plate")
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Slide), xlab = "PC1", ylab = "PC2", main="Slide")
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Array), xlab = "PC1", ylab = "PC2", main="Array")
dev.off()

## Can further test via ANOVA, to look for strong variations in PCs by Batch:
anova(lm(Prin.comp2$PC1~Prin.comp2$Sample_Plate)) 
anova(lm(Prin.comp2$PC2~Prin.comp2$Sample_Plate)) 
anova(lm(Prin.comp2$PC3~Prin.comp2$Sample_Plate)) 
anova(lm(Prin.comp2$PC4~Prin.comp2$Sample_Plate)) 
anova(lm(Prin.comp2$PC5~Prin.comp2$Sample_Plate)) 

#for slide
anova(lm(Prin.comp2$PC1~Prin.comp2$Slide)) 
anova(lm(Prin.comp2$PC2~Prin.comp2$Slide)) 
anova(lm(Prin.comp2$PC3~Prin.comp2$Slide)) 
anova(lm(Prin.comp2$PC4~Prin.comp2$Slide)) 
anova(lm(Prin.comp2$PC5~Prin.comp2$Slide)) 


#for array
anova(lm(Prin.comp2$PC1~Prin.comp2$Array))
anova(lm(Prin.comp2$PC2~Prin.comp2$Array)) 
anova(lm(Prin.comp2$PC3~Prin.comp2$Array)) 
anova(lm(Prin.comp2$PC4~Prin.comp2$Array)) 
anova(lm(Prin.comp2$PC5~Prin.comp2$Array)) 

#second correction for array
mod2 <- model.matrix(~1, data=pd_cleaned)
M_combat2_2array = ComBat(M_combat2_1plate,batch = pd_cleaned$Array, mod = mod2)
save(M_combat2_2array,file="RData/M_combat2_2array.Rdata")       

## Check to see if batch effect was succesfully removed
PCobj2 = prcomp(t(M_combat2_2array), retx = T, center = T, scale. = T)
PCs2 = PCobj2$x
PCs2 =PCs2[,1:R]
Prin.comp2<-merge(PCs2,pd_cleaned, by = "row.names",all=T) 

#### Check whether batches are still distinguished by first and second PC:
pdf(file="Reports/PC_Variation_by_batch2_afterCombat2.pdf")
par(mfrow=c(3,1))
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Sample_Plate), xlab = "PC1", ylab = "PC2", main="Plate")
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Slide), xlab = "PC1", ylab = "PC2", main="Slide")
plot(Prin.comp2$PC1,Prin.comp2$PC2,pch=16, col=as.factor(Prin.comp2$Array), xlab = "PC1", ylab = "PC2", main="Array")
dev.off()


## Can further test via ANOVA, to look for strong variations in PCs by Batch:

# for plate
anova(lm(Prin.comp2$PC1~Prin.comp2$Sample_Plate)) 
anova(lm(Prin.comp2$PC2~Prin.comp2$Sample_Plate)) 
anova(lm(Prin.comp2$PC3~Prin.comp2$Sample_Plate)) 
anova(lm(Prin.comp2$PC4~Prin.comp2$Sample_Plate)) 
anova(lm(Prin.comp2$PC5~Prin.comp2$Sample_Plate)) 

#for slide
anova(lm(Prin.comp2$PC1~Prin.comp2$Slide))
anova(lm(Prin.comp2$PC2~Prin.comp2$Slide)) 
anova(lm(Prin.comp2$PC3~Prin.comp2$Slide)) 
anova(lm(Prin.comp2$PC4~Prin.comp2$Slide)) 
anova(lm(Prin.comp2$PC5~Prin.comp2$Slide)) 

#for array
anova(lm(Prin.comp2$PC1~Prin.comp2$Array)) 
anova(lm(Prin.comp2$PC2~Prin.comp2$Array)) 
anova(lm(Prin.comp2$PC3~Prin.comp2$Array)) 
anova(lm(Prin.comp2$PC4~Prin.comp2$Array)) 
anova(lm(Prin.comp2$PC5~Prin.comp2$Array)) 

#-> ok, leave it like that.

#######################################################
## Convert the batch-adjusted M-values back into betas:

expit2 = function(x) 2^x/(1+2^x)
BMIQ.quantileN_combated = expit2(M_combat2_2array)
dim(BMIQ.quantileN_combated) 
#

## Save normalized and batch-adjusted beta values
save(BMIQ.quantileN_combated,file = "RData/BMIQ.quantileN_combated.Rdata")

#plot final densities
png(file="Reports/BetaValue_Distributions_afterNormCombat2.png",width=700,height=700,pointsize=12)
densityPlot(BMIQ.quantileN_combated, sampGroups = pd_cleaned$Slide, legend=FALSE, main = "PostQC - Normalized and Batch Corrected Beta", xlab = "Beta")
dev.off() 

all.equal(colnames(BMIQ.quantileN_combated),rownames(pd_cleaned)) #TRUE
annotated_pd_clean =new("AnnotatedDataFrame", data= as.data.frame(pd_cleaned)) #extend to AnnotatedDataFrame (required for ESet)

BMIQ.quantileN_combated_ExprSet = new("ExpressionSet", exprs= as.matrix(BMIQ.quantileN_combated), phenoData=annotated_pd_clean)
save(BMIQ.quantileN_combated_ExprSet,file="RData/BMIQ.quantileN_combated_ExprSet.Rdata")


#################################################################################
## checks with MixedUpMapper
#https://github.com/molgenis/systemsgenetics/wiki/Methylation-mixupmapping
#################################################################################
# If you have genotype data available, checks for further mismatches can be done.
# We prepare our methylation data to have it in the correct format.

# I) prepare your final beta file
#attention: format for methylation data: IDs in columns, probes in rows, first column-name is an empty tab!
write.table(Betas_combated, file="files_for_mixup/Betas_combat_mixedupmapper.txt", sep="\t", quote=F, row.names=T, col.names=T)
#use sed to get it in right format (in shell)
# sed 's/NamedererstenIDalsoersteColumn/\tNameIDersteColumn/' alterName.txt > neuerName.txt
# sed 's/203038300047_R01C01/\t203038300047_R01C01/' Betas_combated_mixedupmapper.txt  > Betas_combated_mixedupmapper_final.txt


# II) we need a file genotypemethylationcoupling 
# should be only 2 columns: Sample_Name & ArrayID
genomeID <- data.frame(pd_clean$Sample_Name, pd_clean$Slide, pd_clean$Array)
colnames(genomeID) <- c("Sample_Name", "Slide", "Array")
# columns to paste together
cols <- c('Slide' , 'Array')
# create a new column `ArrayID` with the three columns collapsed together
genomeID$ArrayID <- apply(genomeID[ , cols ] , 1 , paste , collapse = "_" )
genomeID <- genomeID[ , !( names(genomeID) %in% cols ) ]

genotypemethylationcoupling <- genomeID
write.table(genotypemethylationcoupling, file="files_for_mixup/genotypemethylationcoupling.txt", row.names=F, col.names=T, quote =F, sep="\t")

## check with MixedupMapper...

# exclude mixed up probes that were mismatches both re methylation & genotype

#################################################################################

# IF YOU HAVE CORD BLOOD: CHECK CONTAMINATION
# maybe do other data-specific quality checks


#################################################################################
# add if you want:

### GAPHUNTER ###

# Step 4: R-code for detecting outliers using gaphunter
# prepared by: Sinjini Sikdar and adapted by: N Fernandez-Jimenez - nora.fernandez@ehu.eus (16/03/2018)
# Line 74 - modified by MKL(mikyeong.lee@nih.gov) to keep the digits after the decimal ## Dec 12, 2017 
## A final outlier-removed beta file ("Beta4ewas.Rda")

# Required input & parameters
N <- ncol(Betas_combated_mixc) ## No. of samples in your data. 
gapsize <- 0.3    ## Gap size between the outliers and all other values. We have used 0.3 in our analysis.
cutoff <- ifelse(5>(0.25/100)*N,5/N,0.0025)   ## This cutoff is chosen for detecting probes with outliers. We have chosen this
## cutoff such that a probe can have a maximum of 5 or 0.0025 of the total number                                                
## of samples (whichever is larger) as outliers. Can change it if required.

## Note: If gap hunter returns error as it may not detect any outlier, please reduce the gap size or increase the cutoff.
beta_matrix<-Betas_combated_mixc
# Program for gap-hunter
## log files for identifying missing values (NAs) before running gaphunter
probes_log_before <- data.frame(probes=rownames(beta_matrix),no_missing_before=apply(beta_matrix, 1, function(x) sum(is.na(x))))
samples_log_before <- data.frame(samples=colnames(beta_matrix),no_missing_before=apply(beta_matrix, 2, function(x) sum(is.na(x))))

## imputation (as gaphunter cannot handle missing values)
beta1 <- t(apply(beta_matrix,1,function(x) ifelse(is.na(x),median(x,na.rm=T),x)))          

## gaphunter
library(minfi)
gapres <- gaphunter(beta1,keepOutliers=TRUE,threshold=gapsize,outCutoff=cutoff)
# note in excel
gapres1 <- gaphunter(beta1,keepOutliers=FALSE,threshold=gapsize,outCutoff=cutoff)


all_signals_probes <- gapres$proberesults
all_signals_samples <- gapres$sampleresults
all_signals_all <- merge(all_signals_probes,all_signals_samples,by.x="row.names",by.y="row.names")
without_outlier_signals_probes <- gapres1$proberesults
outlier <- setdiff(rownames(all_signals_probes),rownames(without_outlier_signals_probes))  ## probes with outliers
outlier_signals_all <- all_signals_all[which(all_signals_all[,1]%in%outlier),]          
no_incl <- max(outlier_signals_all$Groups)+2
colnames(outlier_signals_all) <- c("probes",colnames(outlier_signals_all)[2:no_incl],colnames(beta1))

new_beta <- NULL
for(p in 1:nrow(outlier_signals_all)){
  group_number <- as.numeric(strsplit(names(which.max(outlier_signals_all[p,3:no_incl])),"")[[1]][6])
  df2 <- outlier_signals_all[p,c(rep(TRUE,no_incl),outlier_signals_all[p,-c(1:no_incl)]!=group_number)]
  sub_int <- colnames(df2)[-c(1:no_incl)]
  beta_out <- beta_matrix[which(rownames(beta_matrix)%in%df2$probes),]
  beta_out[names(beta_out)%in%sub_int] <- NA
  new_beta <- rbind(new_beta,c(probes=df2$probes,beta_out))
}

rownames(new_beta) <- new_beta[,1]
new_beta_1 <- data.frame(new_beta[,-1],check.names=FALSE)
beta2 <- as.data.frame(beta_matrix)
beta2[match(rownames(new_beta_1),rownames(beta2)), ] <- new_beta_1

betas_after_gap <- beta2
# Outputs
## Final beta matrix (saved as betas_after_gap in .Rda format) with all outliers detected as NAs. 
save(betas_after_gap, file="RData/betas_after_gap.Rda")
## We have the same beta matrix as before (input matrix) with additional NAs for extreme outliers.
beta_matrix[is.na(new_beta_1)] <- NA; 
Betas_gapped <- beta_matrix
save(Betas_gapped, file = "RData/Betas_gapped.Rda") ## Use this for further analysis
dim(Betas_gapped) 

load("RData/Betas_gapped.Rda")
## log files in .csv formats. These log files will have columns with:
## probe/sample names: "probes"/"samples", 
## number of missing values (NAs) before running gaphunter: "no_missing_before", 
## number of missing values (NAs including outliers) after running gaphunter: "no_missing after",
## number of outliers identified by gaphunter: "outlier_gaphunter", and 
## the percentage of (non-missing) samples/probes identified as outliers by gaphunter: "percent_outlier_gaphunter".

probes_log_after <- data.frame(probes=rownames(betas_after_gap),no_missing_after=apply(betas_after_gap, 1, function(x) sum(is.na(x))))
log_probes_combined <- merge(probes_log_before,probes_log_after,by.x="probes",by.y="probes",all.x=TRUE,all.y=TRUE)
log_probes_combined$outlier_gaphunter <- log_probes_combined$no_missing_after-log_probes_combined$no_missing_before
log_probes_combined$percent_outlier_gaphunter <- (log_probes_combined$outlier_gaphunter/(ncol(betas_after_gap)-log_probes_combined$no_missing_before))*100
write.table(log_probes_combined,"Gaphunter_CpGs.txt",quote=F,row.names=F,sep="\t")

samples_log_after <- data.frame(samples=colnames(betas_after_gap),no_missing_after=apply(betas_after_gap, 2, function(x) sum(is.na(x))))
log_samples_combined <- merge(samples_log_before,samples_log_after,by.x="samples",by.y="samples",all.x=TRUE,all.y=TRUE)
log_samples_combined$outlier_gaphunter <- log_samples_combined$no_missing_after-log_samples_combined$no_missing_before
log_samples_combined$percent_outlier_gaphunter <- (log_samples_combined$outlier_gaphunter/(nrow(betas_after_gap)-log_samples_combined$no_missing_before))*100
write.csv(log_samples_combined,"Gaphunter_samples.csv")  
write.table(log_samples_combined,"Gaphunter_samples.txt",quote=F,row.names=F,sep="\t")

#######################################################################################################################################

#########################################################################################
# make RGSet with final samples
load("RData/Betas_combated_ExprSet_/CC.Rdata")
laod("RData/RGSet_clean.Rdata")

RGSet_cleaned_final <- RGSet_clean[ , !(colnames(RGSet_qual) %in% sampleNames(Betas_combated_ExprSet))]
RGSet_cleaned_final
save(RGSet_cleaned_final, file = "finalData/RGSet_cleaned_final.Rdata")

#########################################################################################

# NOW YOU CAN ESTIMATE CELL TYPES.

# remember this function if you need M values: mval <- apply(Betas_combated, 2, function(x) log2((x)/(1-x))) # M values
