set.seed(100)
sds=NULL
#uses DEseq2 to do cv and selection

#library(EDASeq)	
#library(corrplot)
library(e1071)
library(DESeq2)
library(MASS)
#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")

###################################################
#  READ  IN DATA          
###################################################

# read in gene expression data (N=120)
geneexpression <- read.table("genetable.txt", header = T, row.names = 1)
metadata <- read.table("Samples104.EDASeqFullBrainFeatures.txt", header = T)
# remove sample outliers (from authors)
geneexpression <- geneexpression[-c(1, 2, 5, 7, 17,21,  26, 54, 58, 83, 87, 92, 94, 95,117, 118)] #################10 samples + 5 samples + 1 (N=104)
#permute to match order of phenotype 
sample_names = c("ba10.s14",	"ba10.s16",	"ba10.s23",	"ba10.s69",	"ba10.s71",	"ba10.s72",	"ba10.s73",	"ba10.s74",	"ba10.s75",	"ba10.s79",	"ba10.s8",	"ba10.s85",	"ba10.s86",	"ba10.s88",	"ba19.s1",	"ba19.s11",	"ba19.s15",	"ba19.s16",	"ba19.s17",	"ba19.s2",	"ba19.s23",	"ba19.s25",	"ba19.s26",	"ba19.s27",	"ba19.s28",	"ba19.s3",	"ba19.s31",	"ba19.s32",	"ba19.s33",	"ba19.s34",	"ba19.s35",	"ba19.s36",	"ba19.s37",	"ba19.s38",	"ba19.s39",	"ba19.s40",	"ba19.s41",	"ba19.s42",	"ba19.s44",	"ba19.s45",	"ba19.s46",	"ba19.s47",	"ba19.s50",	"ba19.s51",	"ba19.s52",	"ba19.s53",	"ba19.s55",	"ba19.s56",	"ba19.s58",	"ba19.s6",	"ba19.s61",	"ba19.s62",	"ba19.s63",	"ba19.s64",	"ba19.s66",	"ba19.s67",	"ba19.s68",	"ba19.s69",	"ba19.s7",	"ba19.s70",	"ba19.s71",	"ba19.s72",	"ba19.s73",	"ba19.s74",	"ba19.s75",	"ba19.s76",	"ba19.s78",	"ba19.s79",	"ba19.s80",	"ba19.s81",	"ba19.s82",	"ba19.s84",	"ba19.s85",	"ba19.s88",	"ba19.s89",	"ba19.s9",	"ba44.s1",	"ba44.s13",	"ba44.s14",	"ba44.s17",	"ba44.s22",	"ba44.s27",	"ba44.s28",	"ba44.s29",	"ba44.s3",	"ba44.s37",	"ba44.s39",	"ba44.s4",	"ba44.s44",	"ba44.s45",	"ba44.s47",	"ba44.s48",	"ba44.s50",	"ba44.s52",	"ba44.s55",	"ba44.s56",	"ba44.s58",	"ba44.s60",	"ba44.s63",	"ba44.s64",	"ba44.s65",	"ba44.s75",	"ba44.s8",	"ba44.s9")
geneexpression<-geneexpression[,sample_names]

###################################################
#  DESEQ2         
###################################################
#make model
dds <- DESeqDataSetFromMatrix(countData=geneexpression, colData = metadata, design =~  brainregion+ Site01 + Age + Sex.01 +Dx.01)
dds$Dx.01 <- factor(dds$Dx.01, levels=c("0","1"))

dds <-estimateSizeFactors(dds)
dds<-estimateDispersions(dds)
dds<-nbinomWaldTest(dds, maxit = 500)
res <- results(dds, alpha = 0.05)
D_DE_gene <- subset(res, res$padj < 0.05)
D_DE_gene = cbind(D_DE_gene, which(res$padj<0.05))
    
write.csv(D_DE_gene, "main_deseq.csv")
