#setwd("/media/IMID/Projects/Project4-Celgene/Data-Storage/RNA_DATA/Raw_counts/")
#setwd("/mnt/md127/Analysis_Projects/Classes_RNAseq_2019/")
setwd("/home/tonig/Classes_RNAseq_2019/Not_CC_adj/")

######################
###########AIM: Find differentially expressed classes in all IMIDS individually using iteratively the next pipeline
############0.1. Split disease in Hi and Low, establish two sets, half hi + half low + half control
            ###0.2 Normalize each set  independtly
            ###0.3 Find DEG and Var genes independtly
            ###0.4 Find classes using M3C algorithm PAM measure using different variants in selecting gene

##############Toni Gomez IMIDOMICS November 2019


##########################


############Loading libs
library(edgeR)
library(limma)

library(gplots)
library(RColorBrewer)
library(M3C)
library(NMF) # loading for heatmap plotting function
library(ggsci)
library(heatmap.2x)## approach
library(parallel)
library(ggplot2)
library(randomForest)
library(parallel)
library(ggplot2)
library(rmarkdown)
library(glmnet)
library(doParallel)
library(HandTill2001)
library(pROC)
library(HandTill2001)
library(caret)

##loading accessory scripts

source("Predictor_scripts/multi_brier.R")
source("Predictor_scripts/multi_logloss.R")
source("Predictor_scripts/RF_utilities.R")
source("Predictor_scripts/makefolds.R")



##loading counts, samplesheet and cellular component

cts <- readRDS("Raw_counts_no_outliers.rds")
annot<-read.csv("Annotation_RNAseq_Table_no_outliers.csv")
##load cellular component

cell.comp<-readRDS("cell_proportions_no_outliers.rds")

# select relevant cell types
predicted.cellcounts <- cell.comp[,c('Granulocytes',
                                                       'B cells (CD19+)',
                                                       'CD4+ T cells',
                                                       'CD8+ T cells',
                                                       'NK cells (CD3- CD56+)',
                                                       'Monocytes (CD14+)')]
# scale to sum to 100
predicted.cellcounts.scaled <- (predicted.cellcounts/rowSums(predicted.cellcounts))*100



####################################
####################################

##PHASE 1: Discovery Set
#######################################################
#######################################################  

##normalize jointly CTRLs, and HI DIS, subset later to deal with analysis
##normalize jointly CTRLs, and LO DIS, subset later to deal with analysis


##deal with all IMIDs

my.imids<-c("RA","PS","PSA","SLE","CD","UC")

for ( r in 1:length(my.imids)){

  main_dir <- my.imids[r]
  
  if (!dir.exists(main_dir)){
    dir.create(main_dir)
  } else {
    print("Dir already exists!")
  }
  
  imid<-my.imids[r]
  
  message ("Analyzing IMID ", my.imids[r])
  
  ##normalize jointly half of the CTRLs and HI

  message ("Loading set1 and set2 ", my.imids[r])

  
  set1<-readRDS(paste0(main_dir,"/","Set1_",imid,".rds"))
  set2<-readRDS(paste0(main_dir,"/","Set2_",imid,".rds"))
  
  A.set1<-as.character(set1[!(set1$DIS=="CTRL"),"set1"])
  B.hi<-as.character(set1[set1$DIS=="CTRL","set1"])
  
  A.set2<-as.character(set2[!(set2$DIS=="CTRL"),"set2"])
  B.lo<-as.character(set2[set2$DIS=="CTRL","set2"])
  
   outliers <- c("ix0007847","ix0006088", "ix0005117", "ix0001040")

  A.set1<-A.set1[!(A.set1 %in% outliers)]
  A.set2<-A.set2[!(A.set2 %in% outliers)]

  B.hi<-B.hi[!(B.hi %in% outliers)]
  B.lo<-B.lo[!(B.lo %in% outliers)]  

message ("Normalizing set1 ", my.imids[r])
  
  
  ########
  
   
  my.cts<-cts[,colnames(cts) %in% c(A.set1,B.hi,A.set2,B.lo)]

  samples2test<-annot[annot$X %in% c(A.set1,B.hi,A.set2,B.lo),]

  group<-samples2test$imid
  group<-droplevels.factor(group)

  y <- DGEList(my.cts,group = group)

  keep <- rowSums(cpm(y)>2) >= length(c(A.set1,A.set2))#### set number for all samples in group, no matter CTRL because the size is the same,
## and we want to check at least 2 cpms across group
  y <- y[keep, ]

##Normalize using TMM

  y$samples$lib.size <- colSums(y$counts)
  y <- calcNormFactors(y, method="TMM")


 ##get normalized counts in RPKM

  gene.length<-read.delim("Gene_length_hg19.txt",head=T)

 ##same order as counts

  gene.length<-gene.length[match(rownames(y$counts),gene.length$Geneid),]

# reading gene length file as a matrix
  y$genes$Length<- gene.length$Length
# rpkm function
 # my.rpkm <- rpkm(y)
  #my.rpkm<-log2(my.rpkm+1)
  
  my.selection<-c("SCL_logCPM")
  ####selection
  
 # my.selection<-c("SCL_logCPM","SCL_logRPKM")
  my.rpkm <- cpm(y,log = T,prior.count = 2)
  
  my.rpkm <- t(scale(t(my.rpkm))) #The scale function calculates column z-scores, but since it is used as t (transpose), it is actually performed on rows. 
                                  ##Then we scale each row (each gene) to have mean zero and standard deviation one:
  
##remove batch,sex,age and cellular component using linear model
  samples2test$Sex<-ifelse(samples2test$sex=="FEMALE","F","M")

##BATCH

  my.batch<-paste("B", samples2test$new.plate,sep="")
  cell.counts<-predicted.cellcounts.scaled[rownames(predicted.cellcounts.scaled) %in% c(A.set1,B.hi,A.set2,B.lo),]
  cell.counts<-cell.counts[match(samples2test$X,rownames(cell.counts)),]

  ##design <- model.matrix(~group+samples2test$Sex + my.batch  + samples2test$age + as.matrix(cell.counts[,-6]))##monocytes error, out due to next
  design <- model.matrix(~group+samples2test$Sex + my.batch  + samples2test$age ) 
 rownames(design)<-samples2test$X

  fit <- lmFit(my.rpkm[,colnames(my.rpkm) %in% rownames(design)], design)
  mAdj.fit    <- fit$coefficients[,-c(1,2)]##extract coefficients to keep

  all.adj<- as.matrix(my.rpkm[,colnames(my.rpkm) %in% rownames(design)]) - mAdj.fit %*% t(design[,-c(1,2)])##subtract the cofounders
  
  hi.adj<-all.adj[,colnames(all.adj) %in% A.set1]
  lo.adj<-all.adj[,colnames(all.adj) %in% A.set2]
  
  rm(fit)
  rm(mAdj.fit)
  rm(all.adj)
  rm(y)
  rm(samples2test)
  rm(my.cts)
  rm(group)
  
  gc()   
  
  
###RF predictions
  
  ###########RF procedure to test classes
  ####Construct model on Hi and predict on LOW
  message ("Test class in set2 samples ",imid)
  
  ###load 
  my.options<-c("DEG","DEG_VAR","VAR","OT","OT_DEG","OT_DEG_VAR","OT_DEG_IMM","OT_IMM","VAR_IMM","VAR_OT_IMM","OT_VAR","IMM","IMM_DEG")
  ##01 Use classes from M3C
  
   for (i in 1:length(my.options)){ 
     
     sub_dir <-paste0(my.options[i],"_",my.selection)
     
     output_dir <- file.path(imid, sub_dir)
     if (!dir.exists(output_dir)){
       dir.create(output_dir)
     } else {
       print("Dir already exists!")
     }
     
     print(i)
    annon<-read.csv(paste0(main_dir,"/",my.options[i],"_",my.selection,"/","Annon_classes_Set1",imid,"_",my.options[i],".csv"))
  
      classes <- paste("C",annon$consensuscluster,sep="_")
      classes<-as.factor(classes)
  
  ##02 load HI data and get DEG genes between classes
  ###avoid problems with one class only
    annon$ID<-annon$X
    data.tp<-hi.adj[,colnames(hi.adj) %in% annon$ID]
    
    if(min(table(annon$consensuscluster))<2){
      
      class.remov<-names(table(annon$consensuscluster)[table(annon$consensuscluster)==1])
      ###sink in a file the outliers
      outliers<-annon[annon$consensuscluster %in% class.remov,"X"]
      
      sink(paste0(main_dir,"/","Outliers",imid,"_",my.selection,".csv",sep=""),append = TRUE)  
      print(outliers)
      sink()  
      annon<-annon[!(annon$consensuscluster %in% class.remov),]
      data.tp<-data.tp[,colnames(data.tp) %in% annon$X]
      classes <- paste("C",annon$consensuscluster,sep="_")
      classes<-as.factor(classes)
      
      
    }
    
    if(dim(table(classes))==1){
      
      sink(paste0(main_dir,"/",my.options[i],"_",my.selection,"/","Warning",imid,"_",my.options[i],"_",my.selection,".csv",sep=""))  
      print("Only one class predicted!")
      sink()  
      
    }else{
  data.tp<-data.tp[,match(annon$ID,colnames(data.tp))]
  
  ##AOV to identify genes that are variable among and between groups.
  message("Identifying genes that are variable among and between classes groups ",imid,"_",my.options[i]," ...",Sys.time())
  
  pv <- sapply(1:dim(data.tp)[1], function(i) {
    mydataframe <- data.frame(y=data.tp[i,], ig=classes)
    fit <- aov(y ~ ig, data=mydataframe)
    summary(fit)[[1]][["Pr(>F)"]][1]
  })
  names(pv) <- rownames(data.tp)
  pv.sig <- names(pv)[pv < 0.05/1:dim(data.tp)[1]/length(classes)] ## bonferonni
  p<-ifelse(length(pv.sig)>250,250,length(pv.sig))
  
  ##using selected genes diff express
  
  to.train<-data.tp[rownames(data.tp) %in% pv.sig,]
  ##03 Build RF model
  
  y <- classes
  ####select genes in each external fold, select the one with less error rate and construct RF
  
  ntrees <- 999
  
  genes<-list()
  errors<-c()
  rf.list<-list()
  cores<-8
  folds<-5
  folds.genes <- makefolds(y,5)
  
  
  for(K in 1:folds){
    
    message("calculating outer fold for selecting genes",K,"  ...",Sys.time())
    
    fold <- folds.genes[[K]]
    
    rf.scores <- calcultateCVfold.genes(t(to.train),as.factor(classes),fold,p,cores,ntrees)
    genes[[K]]<-rf.scores[[2]]
    errors[K]<-rf.scores[[3]]
    rf.list[[K]]<-rf.scores[[4]]
    
  }
  
  ##04 select those one with less errors
  
  sel<-which(errors==min(errors))[1]##index of selection
  
  ##06 get statistics
  
  to.build<-as.data.frame(to.train[rownames(to.train) %in%  rownames(rf.list[[sel]]$importance),])
  to.build<-as.data.frame(t(to.build))
  to.build$predicted.response<-predict(rf.list[[sel]], to.build, type="response")
  to.build$predicted.response<-ordered(to.build$predicted.response)
  to.build$real<-classes
  to.build$real<-ordered(to.build$real)
  
  
  stats.preds<-confusionMatrix(to.build$predicted.response,  
                               to.build$real)
  write.csv(as.data.frame(stats.preds$overall),paste0(main_dir,"/",my.options[i],"_",my.selection,"/","Stats_RF_set1","_",imid,".csv",sep=""))  
  
  ##05 predict
  
  ##Select set2 samples to predict
  
  
  predictor_data<-t(lo.adj[rownames(lo.adj) %in% rownames(rf.list[[sel]]$importance),])##
  
  RF_predictions_probs<-predict(rf.list[[sel]], predictor_data, type="p")
  
  RF_predictions_responses<-predict(rf.list[[sel]], predictor_data, type="response")
  
  RF_p<-as.data.frame(RF_predictions_probs)
  RF_a<-as.data.frame(RF_predictions_responses)
  
  #####print results
  write.csv(RF_p,paste0(main_dir,"/",my.options[i],"_",my.selection,"/","Probs_prediction_set1",imid,".csv",sep="")) 
  
   ################################################################
      ########################################################
      gc()
      
}
 }
  
}
