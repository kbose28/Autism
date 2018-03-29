library(DESeq2)
library(MASS)
library(Rcpp)
library(FarmTest) 
library(RcppArmadillo)
#C++ functions that are modifications of the package, for comparison
sourceCpp("Cfunctions.cpp")

###################################################
#  READ  IN DATA          
###################################################

# read in gene expression data (N=120)
geneexpression <- read.table("genetable.txt", header = T, row.names = 1)
metadata <- read.table("Samples104.EDASeqFullBrainFeatures.txt", header = T)
# remove sample outliers (from authors)
geneexpression <- geneexpression[-c(1, 2, 5, 7, 17,21,  26, 54, 58, 83, 87, 92, 94, 95,117, 118)] #################10 samples + 5 samples + 1 (N=104)
#permute to match order of phenotype 
sample_names = c("ba10.s14",    "ba10.s16",     "ba10.s23",     "ba10.s69",     "ba10.s71",     "ba10.s72",     "ba10.s73",     "ba10.s74",     "ba10.s75",     "ba10.s79",     "ba10.s8",      "ba10.s85",     "ba10.s86",     "ba10.s88",     "ba19.s1",      "ba19.s11",     "ba19.s15",     "ba19.s16",     "ba19.s17",     "ba19.s2",      "ba19.s23",     "ba19.s25",     "ba19.s26",     "ba19.s27",     "ba19.s28",     "ba19.s3",      "ba19.s31",     "ba19.s32",     "ba19.s33",     "ba19.s34",     "ba19.s35",     "ba19.s36",     "ba19.s37",     "ba19.s38",     "ba19.s39",     "ba19.s40",     "ba19.s41",     "ba19.s42",     "ba19.s44",     "ba19.s45",     "ba19.s46",     "ba19.s47",     "ba19.s50",     "ba19.s51",     "ba19.s52",     "ba19.s53",     "ba19.s55",     "ba19.s56",     "ba19.s58",     "ba19.s6",      "ba19.s61",     "ba19.s62",     "ba19.s63",     "ba19.s64",     "ba19.s66",     "ba19.s67",     "ba19.s68",     "ba19.s69",     "ba19.s7",      "ba19.s70",     "ba19.s71",     "ba19.s72",     "ba19.s73",     "ba19.s74",     "ba19.s75",     "ba19.s76",     "ba19.s78",     "ba19.s79",     "ba19.s80",     "ba19.s81",     "ba19.s82",     "ba19.s84",     "ba19.s85",     "ba19.s88",     "ba19.s89",     "ba19.s9",      "ba44.s1",      "ba44.s13",     "ba44.s14",     "ba44.s17",     "ba44.s22",     "ba44.s27",     "ba44.s28",     "ba44.s29",     "ba44.s3",      "ba44.s37",     "ba44.s39",     "ba44.s4",      "ba44.s44",     "ba44.s45",     "ba44.s47",     "ba44.s48",     "ba44.s50",     "ba44.s52",     "ba44.s55",     "ba44.s56",     "ba44.s58",     "ba44.s60",     "ba44.s63",     "ba44.s64",     "ba44.s65",     "ba44.s75",     "ba44.s8",      "ba44.s9")
geneexpression<-geneexpression[,sample_names]
Humanfeature <- read.table("P10.70.19.GC.GeneLength.txt", header = TRUE, row.names = 1)
Humanfeature$GC <- Humanfeature$GC/100
geneexpression = geneexpression[,metadata$Sex.01==1]
###################################################
#  DESEQ2         
###################################################
#make model
dds <- DESeqDataSetFromMatrix(countData=geneexpression, colData = metadata, design =~  brainregion+Sex.01+ Site01 + Age  +Dx.01)
dds$Dx.01 <- factor(dds$Dx.01, levels=c("0","1"))
  dds <- estimateSizeFactors(dds)
  Normdata <- counts(dds, normalized = TRUE)
  index <- rowSums(Normdata >= (10))>= (0.9*dim(Normdata)[2])
  Normdata<-  Normdata[index,]
  D_counts = log2(Normdata+1)
  p = NROW(D_counts)
  n = NCOL(D_counts)

   Z <- cbind(matrix(1,n , 1) , metadata$Dx.01,  metadata$Site01, metadata$Age, metadata$brainregion, metadata$Sex.01)

   coef = mu_robust_F(0.5,matrix(D_counts , p, n), matrix(Z, n, NCOL(Z)))
 write.csv(coef, "coefficients.csv")
   y_leftover <- D_counts-t(coef)%*% t(Z)+ (coef[2,])%*%t(metadata$Dx.01)

  X = (y_leftover[,metadata[,"Dx.01"]==1])
  Y = (y_leftover[,metadata[,"Dx.01"]==0])
  
  nx = NCOL(X)
  ny = NCOL(Y)


  muhatx = mu_robust(0.5, matrix(X, p, nx))#the first term is redundant, using CV
#  covx = Cov_Huber(0.6,  matrix(t(X-rowMeans(X)%*%t(rep(1,nx))),nx,p), matrix(rep(0,nx), nx,1))
  covx = cov(X)
print(dim(covx))
  eigs = Eigen_Decomp(covx)
  values = eigs[,nx+1]
  vectors = eigs[,1:nx]  
  values = pmax(values,0)
  ratio=c()
  for(i in 1:(floor(min(nx,p)/2))){
     ratio=append(ratio, values[i+1]/values[i])}
  ratio = ratio[is.finite(ratio)]
  Kx = which.min(ratio)
print(Kx)
  Bx = matrix(NA, p, Kx)


  for (k in 1:Kx){
    Bx[,k] = sqrt(values[k])*(X-rowMeans(X)%*%t(rep(1,nx)))%*%vectors[,k]/sqrt(nx*values[k])
  }
write.csv(Bx, "Xcoefs.csv")
  Bx2 = apply(Bx,1, function(y) sum(y^2))
  thetax = mu_robust(0.5, matrix(X^2, p, nx))#the first term is redundant, using CV
  varhatx_0 = ( thetax - muhatx^2)* ( thetax > muhatx^2) +(thetax)* ( thetax <=muhatx^2)
  varhatx = (varhatx_0 - Bx2)* (varhatx_0 > Bx2) +(varhatx_0)* ( varhatx_0 <=Bx2)
  sehatx = sqrt(varhatx/nx)
  fx = mu_robust_F(0.5, matrix(rowMeans(X),1, p), matrix(Bx, p, Kx))
  
  muhaty = mu_robust(0.5, matrix(Y, p, ny))#the first term is redundant, using CV
#  covy = Cov_Huber(0.6,  matrix(t(Y-rowMeans(Y)%*%t(rep(1,ny))),ny,p), matrix(rep(0,ny), ny,1))
   covy = cov(Y)
  eigs = Eigen_Decomp(covy)
  values = eigs[,ny+1]
  vectors = eigs[,1:ny]
  values = pmax(values,0)
  ratio=c()
  for(i in 1:(floor(min(ny,p)/2))){
     ratio=append(ratio, values[i+1]/values[i])}
  ratio = ratio[is.finite(ratio)]
  Ky = which.min(ratio)
  print(Ky)
  By = matrix(NA, p, Ky)
  for (k in 1:Ky){
    By[,k] = sqrt(values[k])*(Y-rowMeans(Y)%*%t(rep(1,ny)))%*%vectors[,k]/sqrt(ny*values[k])
  }

write.csv(By, "Ycoefs.csv")

  By2 = apply(By,1, function(y) sum(y^2))
  thetay = mu_robust(0.5, matrix(Y^2, p, ny))#the first term is redundant, using CV
  varhaty_0 = ( thetay - muhaty^2)* ( thetay > muhaty^2) +(thetay)* ( thetay <=muhaty^2)
  varhaty = (varhaty_0 - By2)* (varhaty_0 > By2) +(varhaty_0)* ( varhaty_0 <=By2)
  sehaty = sqrt(varhaty/ny)
  fy = mu_robust_F(0.5, matrix(rowMeans(Y),1, p), matrix(By, p, Ky))
means = (muhatx-Bx%*%fx-muhaty +By%*%fy)
  stat=(muhatx-Bx%*%fx-muhaty +By%*%fy)/sqrt(sehatx^2 + sehaty^2)
  pvalue = 2*stats::pnorm(-abs(stat))
  rejected.alldata = farm.FDR(pvalue, alpha = 0.05)
   detected = rejected.alldata$rejected
  genes = rownames(D_counts)
  R_DE_gene = cbind(genes[detected[,1]],detected[,1], detected[,2], detected[,3], means[detected[,1]])

  found = NROW(R_DE_gene)
  print(found)
 write.csv(R_DE_gene, "farm_main_male.csv")







