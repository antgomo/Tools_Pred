#setwd("/home/tonig/Classes_RNAseq_2020")
#setwd("/Users/agomez/IMIDOMICS/Classes_RNAseq_2020/Not_CC_adj/")

setwd("/mnt/md127/Analysis_Projects/Classes_RNAseq_2019/Not_CC_adj")

######################
###########AIM: Find differentially expressed classes in all IMIDS individually using iteratively the next pipeline
############0.1. Split disease in Hi and Low, establish two sets, half hi + half low + half control
            ###0.2 Normalize each set  independtly
            ###0.3 Find DEG and Var genes independtly
            ###0.4 Find classes using M3C algorithm PAM measure using different variants in selecting gene

##############Toni Gomez IMIDOMICS June 2020


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


#load OT genes

ot.sle<-read.csv("OpenTargets/targets_associated_with_systemic_lupus_erythematosus.csv")
ot.sle<-as.character(ot.sle$target.gene_info.symbol)

ot.cd<-read.csv("OpenTargets/targets_associated_with_Crohn's_disease.csv")
colnames(ot.cd)[1]<-"target.gene_info.symbol"
ot.cd<-as.character(ot.cd$target.gene_info.symbol)

ot.uc<-read.csv("OpenTargets/targets_associated_with_ulcerative_colitis.csv")
colnames(ot.uc)[1]<-"target.gene_info.symbol"
ot.uc<-as.character(ot.uc$target.gene_info.symbol)

ot.ra<-read.csv("OpenTargets/targets_associated_with_rheumatoid_arthritis.csv")
ot.ra<-as.character(ot.ra$target.gene_info.symbol)

ot.psa<-read.csv("OpenTargets/targets_associated_with_psoriatic_arthritis.csv")
colnames(ot.psa)[1]<-"target.gene_info.symbol"
ot.psa<-as.character(ot.psa$target.gene_info.symbol)

ot.ps<-read.csv("OpenTargets/targets_associated_with_psoriasis.csv")
colnames(ot.ps)[1]<-"target.gene_info.symbol"
ot.ps<-as.character(ot.ps$target.gene_info.symbol)

###load imm genes sets

imm<-readRDS("go_immune_genes.rds")


####corticoids

cortic<-readRDS("/media/IMID/Projects/Project4-Celgene/Data-Storage/Phenotype-Data/All_RNA_samples/As_predictors/PHENODATA_PRED.rds")
cortic<-cortic[,c("indivId","IMID","IMID_Corticoids")]


###FUNCTION DEFINITION

source("Function_classes.R")

####Filter one select HI classes OR Hi + Low

####################################
####################################

##PHASE 1: Discovery Set
#######################################################
#######################################################  

##normalize jointly CTRLs, and HI DIS, subset later to deal with analysis
##normalize jointly CTRLs, and LO DIS, subset later to deal with analysis


##deal with all IMIDs

#my.imids<-c("RA","PS","PSA","SLE","CD","UC")
my.imids<-c("PSA")

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

  message ("Normalizing HI and half of the controls", my.imids[r])
  my.selection<-c("RPKM")
  
  
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

  ########
  
   
  my.cts<-cts[,colnames(cts) %in% c(A.set1,B.hi)]

  samples2test<-annot[annot$X %in% c(A.set1,B.hi),]

  group<-samples2test$imid
  group<-droplevels.factor(group)

  y <- DGEList(my.cts,group = group)

  keep <- rowSums(cpm(y)>2) >= length(A.set1)#### set number for all samples in group, no matter CTRL because the size is the same,
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
  my.rpkm <- rpkm(y)
  my.rpkm<-log2(my.rpkm+1)
  
  ####selection
  
 # my.selection<-c("SCL_logCPM","SCL_logRPKM")
  my.selection<-c("RPKM")
  
  #my.rpkm <- cpm(y,log = T,prior.count = 2)
  
 # my.rpkm <- t(scale(t(my.rpkm))) #The scale function calculates column z-scores, but since it is used as t (transpose), it is actually performed on rows. 
                                  ##Then we scale each row (each gene) to have mean zero and standard deviation one:
  
##remove batch,sex,age and cellular component using linear model
  samples2test$Sex<-ifelse(samples2test$sex=="FEMALE","F","M")

##BATCH

  my.batch<-paste("B", samples2test$new.plate,sep="")
  cell.counts<-predicted.cellcounts.scaled[rownames(predicted.cellcounts.scaled) %in% c(A.set1,B.hi),]
  cell.counts<-cell.counts[match(samples2test$X,rownames(cell.counts)),]
  cortic<-cortic[cortic$IMID==imid,] 
  cortic2<-cortic[cortic$indivId %in% c(A.set1,B.hi),]
  cortic2<-cortic2[match(samples2test$X,cortic2$indivId),]

  design <- model.matrix(~group+samples2test$Sex + my.batch  + samples2test$age + cortic2$IMID_Corticoids)##monocytes error, out due to next
  
  rownames(design)<-samples2test$X

  fit <- lmFit(my.rpkm[,colnames(my.rpkm) %in% rownames(design)], design)
  mAdj.fit    <- fit$coefficients[,-c(1,2)]##extract coefficients to keep

  hi.adj<- as.matrix(my.rpkm[,colnames(my.rpkm) %in% rownames(design)]) - mAdj.fit %*% t(design[,-c(1,2)])##subtract the cofounders

  rm(fit)
  rm(mAdj.fit)
  rm(y)
  rm(samples2test)
  rm(my.cts)
  rm(group)
  
  gc()   
  
  ##same procedure to deal with LOW
  
  my.cts<-cts[,colnames(cts) %in% c(A.set2,B.lo)]
  
  samples2test<-annot[annot$X %in% c(A.set2,B.lo),]
  
  group<-samples2test$imid
  group<-droplevels.factor(group)
  
  y <- DGEList(my.cts,group = group)
  
  keep <- rowSums(cpm(y)>2) >= length(A.set2)#### set number for all samples in group, no matter CTRL because the size is the same,
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
  my.rpkm <- rpkm(y)
  my.rpkm<-log2(my.rpkm+1)
  
  ####selection
  
  # my.selection<-c("SCL_logCPM","SCL_logRPKM")
  my.selection<-c("RPKM")
   #my.rpkm <- cpm(y,log = T,prior.count = 2)
  
    #my.rpkm <- t(scale(t(my.rpkm))) #The scale function calculates column z-scores, but since it is used as t (transpose), it is actually performed on rows. 
  ##Then we scale each row (each gene) to have mean zero and standard deviation one:
  
  ##remove batch,sex,age and cellular component using linear model
  samples2test$Sex<-ifelse(samples2test$sex=="FEMALE","F","M")
  ##BATCH
  
  my.batch<-paste("B", samples2test$new.plate,sep="")
  cell.counts<-predicted.cellcounts.scaled[rownames(predicted.cellcounts.scaled) %in% c(A.set2,B.lo),]
  cell.counts<-cell.counts[match(samples2test$X,rownames(cell.counts)),]
  
  cortic2<-cortic[cortic$indivId %in% c(A.set2,B.lo),]
  cortic2<-cortic2[match(samples2test$X,cortic2$indivId),]

  design <- model.matrix(~group+samples2test$Sex + my.batch  + samples2test$age + cortic2$IMID_Corticoids)##monocytes error, out due to next

    
  rownames(design)<-samples2test$X
  
  fit <- lmFit(my.rpkm[,colnames(my.rpkm) %in% rownames(design)], design)
  mAdj.fit    <- fit$coefficients[,-c(1,2)]##extract coefficients to keep
  
  lo.adj<- as.matrix(my.rpkm[,colnames(my.rpkm) %in% rownames(design)]) - mAdj.fit %*% t(design[,-c(1,2)])##subtract the cofounders
  

##remove variables to avoid memory problems

  rm(fit)
  rm(mAdj.fit)
  rm(y)
  rm(samples2test)

  rm(my.cts)
  rm(group)

  gc()   


###DEG and VAR tests
  
  
  message ("Discovery set DEG HI ",my.imids[r])


  res.hi<-read.csv(paste0(main_dir,"/","DEG_Set1_",my.imids[r],".csv",sep=""))
  rownames(res.hi)<-res.hi$X
  res.hi<-res.hi[,-1]

  gc()
  
  ###############################################Now, get DEG LO
  
  message ("Discovery set Set2 ",my.imids[r])
  
  res.lo<-read.csv(paste0(main_dir,"/","DEG_Set2_",my.imids[r],".csv",sep=""))
  rownames(res.lo)<-res.lo$X
  res.lo<-res.lo[,-1]
  
  gc()
  
  
  ##VARIABLE genes
  message ("Finding most variable genes between CTRL and set1 ",my.imids[r])

  
  var.sel.hi<-read.csv(paste0(main_dir,"/","VAR_genes_Set1_",my.imids[r],".csv",sep=""))
  var.sel.hi<-as.character(var.sel.hi$x)


  gc()
  ##collect not DEG
  message ("Finding most variable genes between CTRL and set2 ",my.imids[r])
  
  var.sel.lo<-read.csv(paste0(main_dir,"/","VAR_genes_Set2_",my.imids[r],".csv",sep=""))
  var.sel.lo<-as.character(var.sel.lo$x)
  gc()
  ##Try to find classes using 3 options, 
  #################################################################3

  
  my.options<-c("DEG","DEG_VAR","VAR","OT","OT_DEG","OT_DEG_VAR","OT_DEG_IMM","OT_IMM","VAR_IMM","VAR_OT_IMM","OT_VAR","IMM","IMM_DEG")
  

  for (i in 1:length(my.options)){
  
  
  if(my.options[i]=="DEG"){
    
    ##specify dirs to write output results
    

    sub_dir <-paste0(my.options[i],"_",my.selection)
    output_dir <- file.path(imid, sub_dir)
    
    if (!dir.exists(output_dir)){
      dir.create(output_dir)
    } else {
      print("Dir already exists!")
    }
    
    ###define function to work in
    ##define sets to search classes
    genes.sel.hi<-rownames(res.hi)
    genes.sel.lo<-rownames(res.lo)
    
    
    #run function
    class.pipeline(hi.adj,lo.adj,A.set1,A.set2,imid,genes.sel.hi,genes.sel.lo)
    
gc()
################################################################
    
    gc()
}else if (my.options[i]=="OT"){
    #################################################################
  message ("OT pipeline for class Discovery ",imid)
  
  sub_dir <-paste0(my.options[i],"_",my.selection)
  output_dir <- file.path(imid, sub_dir)
  
  #load OT genes
  
  if(imid=="SLE"){
    
    ot.genes<-ot.sle
  
  }else if(imid=="CD"){
    
    ot.genes<-ot.cd
    
    
  } else if(imid=="UC"){
    
    ot.genes<-ot.uc
    
    
  }  else if(imid=="PS"){
    
    ot.genes<-ot.ps
    
    
  }  else if(imid=="RA"){
    
    ot.genes<-ot.ra
    
    
  }else{
    
    ot.genes<-ot.psa
    
    
  }
  
  
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  } else {
    print("Dir already exists!")
  }
  genes.sel.hi<-ot.genes
  genes.sel.lo<-ot.genes
  
  class.pipeline(hi.adj,lo.adj,A.set1,A.set2,imid,genes.sel.hi,genes.sel.lo)
  
  
  ################################################################
  

}else if (my.options[i]=="OT_DEG"){
  #################################################################
  message ("OT_DEG pipeline for class Discovery ",imid)
  
  sub_dir <-paste0(my.options[i],"_",my.selection)
  output_dir <- file.path(imid, sub_dir)
  
  #load OT genes
  
  if(imid=="SLE"){
    
    ot.genes<-ot.sle
    
  }else if(imid=="CD"){
    
    ot.genes<-ot.cd
    
    
  } else if(imid=="UC"){
    
    ot.genes<-ot.uc
    
    
  }  else if(imid=="PS"){
    
    ot.genes<-ot.ps
    
    
  }  else if(imid=="RA"){
    
    ot.genes<-ot.ra
    
    
  }else{
    
    ot.genes<-ot.psa
    
    
  }
  
  
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  } else {
    print("Dir already exists!")
  }
  
  genes.sel.hi<-c(ot.genes,rownames(res.hi))
  genes.sel.lo<-c(ot.genes,rownames(res.lo))
  
  class.pipeline(hi.adj,lo.adj,A.set1,A.set2,imid,genes.sel.hi,genes.sel.lo)
 ############################################

  
  
}else if (my.options[i]=="OT_DEG_VAR"){
  #################################################################
  message ("OT_DEG_VAR pipeline for class Discovery ",imid)
  
  sub_dir <-paste0(my.options[i],"_",my.selection)
  output_dir <- file.path(imid, sub_dir)
  
  #load OT genes
  
  if(imid=="SLE"){
    
    ot.genes<-ot.sle
    
  }else if(imid=="CD"){
    
    ot.genes<-ot.cd
    
    
  } else if(imid=="UC"){
    
    ot.genes<-ot.uc
    
    
  }  else if(imid=="PS"){
    
    ot.genes<-ot.ps
    
    
  }  else if(imid=="RA"){
    
    ot.genes<-ot.ra
    
    
  }else{
    
    ot.genes<-ot.psa
    
    
  }
  
  
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  } else {
    print("Dir already exists!")
  }
  
  genes.sel.hi<-c(ot.genes,rownames(res.hi),var.sel.hi)
  genes.sel.lo<-c(ot.genes,rownames(res.lo),var.sel.lo)
  
  class.pipeline(hi.adj,lo.adj,A.set1,A.set2,imid,genes.sel.hi,genes.sel.lo)
  
  
}else if (my.options[i]=="OT_DEG_IMM"){
  #################################################################
  message ("OT_DEG_IMM pipeline for class Discovery ",imid)
  
  sub_dir <-paste0(my.options[i],"_",my.selection)
  output_dir <- file.path(imid, sub_dir)
  
  #load OT genes
  
  if(imid=="SLE"){
    
    ot.genes<-ot.sle
    
  }else if(imid=="CD"){
    
    ot.genes<-ot.cd
    
    
  } else if(imid=="UC"){
    
    ot.genes<-ot.uc
    
    
  }  else if(imid=="PS"){
    
    ot.genes<-ot.ps
    
    
  }  else if(imid=="RA"){
    
    ot.genes<-ot.ra
    
    
  }else{
    
    ot.genes<-ot.psa
    
    
  }
  
  
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  } else {
    print("Dir already exists!")
  }
  
  
  genes.sel.hi<-c(ot.genes,rownames(res.hi),imm)
  genes.sel.lo<-c(ot.genes,rownames(res.lo),imm)
  
  class.pipeline(hi.adj,lo.adj,A.set1,A.set2,imid,genes.sel.hi,genes.sel.lo)
  
}else if (my.options[i]=="OT_IMM"){
  #################################################################
  message ("OT_IMM pipeline for class Discovery ",imid)
  
  sub_dir <-paste0(my.options[i],"_",my.selection)
  output_dir <- file.path(imid, sub_dir)
  
  #load OT genes
  
  if(imid=="SLE"){
    
    ot.genes<-ot.sle
    
  }else if(imid=="CD"){
    
    ot.genes<-ot.cd
    
    
  } else if(imid=="UC"){
    
    ot.genes<-ot.uc
    
    
  }  else if(imid=="PS"){
    
    ot.genes<-ot.ps
    
    
  }  else if(imid=="RA"){
    
    ot.genes<-ot.ra
    
    
  }else{
    
    ot.genes<-ot.psa
    
    
  }
  
  
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  } else {
    print("Dir already exists!")
  }
  
  genes.sel.hi<-c(ot.genes,imm)
  genes.sel.lo<-c(ot.genes,imm)
  
  class.pipeline(hi.adj,lo.adj,A.set1,A.set2,imid,genes.sel.hi,genes.sel.lo)
 
}else if (my.options[i]=="VAR_IMM"){
  #################################################################
  message ("VAR_IMM pipeline for class Discovery ",imid)
  
  sub_dir <-paste0(my.options[i],"_",my.selection)
  output_dir <- file.path(imid, sub_dir)
  
  
  
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  } else {
    print("Dir already exists!")
  }
  
  genes.sel.hi<-c(var.sel.hi,imm)
  genes.sel.lo<-c(var.sel.lo,imm)
  
  class.pipeline(hi.adj,lo.adj,A.set1,A.set2,imid,genes.sel.hi,genes.sel.lo)
 
}else if (my.options[i]=="VAR_OT_IMM"){
  #################################################################
  message ("VAR_OT_IMM pipeline for class Discovery ",imid)
  
  sub_dir <-paste0(my.options[i],"_",my.selection)
  output_dir <- file.path(imid, sub_dir)
  
  #load OT genes
  
  if(imid=="SLE"){
    
    ot.genes<-ot.sle
    
  }else if(imid=="CD"){
    
    ot.genes<-ot.cd
    
    
  } else if(imid=="UC"){
    
    ot.genes<-ot.uc
    
    
  }  else if(imid=="PS"){
    
    ot.genes<-ot.ps
    
    
  }  else if(imid=="RA"){
    
    ot.genes<-ot.ra
    
    
  }else{
    
    ot.genes<-ot.psa
    
    
  }
  
  
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  } else {
    print("Dir already exists!")
  }
  
  genes.sel.hi<-c(ot.genes,var.sel.hi,imm)
  genes.sel.lo<-c(ot.genes,var.sel.lo,imm)
  
  class.pipeline(hi.adj,lo.adj,A.set1,A.set2,imid,genes.sel.hi,genes.sel.lo)
  
}else if (my.options[i]=="IMM"){
  #################################################################
  message ("IMM pipeline for class Discovery ",imid)
  
  sub_dir <-paste0(my.options[i],"_",my.selection)
  output_dir <- file.path(imid, sub_dir)
  
  
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  } else {
    print("Dir already exists!")
  }
  
  genes.sel.hi<-c(imm)
  genes.sel.lo<-c(imm)
  
  class.pipeline(hi.adj,lo.adj,A.set1,A.set2,imid,genes.sel.hi,genes.sel.lo)
  
}else if (my.options[i]=="IMM_DEG"){
  #################################################################
  message ("IMM_DEG pipeline for class Discovery ",imid)
  
  sub_dir <-paste0(my.options[i],"_",my.selection)
  output_dir <- file.path(imid, sub_dir)
  
  
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  } else {
    print("Dir already exists!")
  }
  
  genes.sel.hi<-c(rownames(res.hi),imm)
  genes.sel.lo<-c(rownames(res.lo),imm)
  
  class.pipeline(hi.adj,lo.adj,A.set1,A.set2,imid,genes.sel.hi,genes.sel.lo)
  
  
}else if (my.options[i]=="OT_VAR"){
  #################################################################
  message ("OT_VAR pipeline for class Discovery ",imid)
  
  sub_dir <-paste0(my.options[i],"_",my.selection)
  output_dir <- file.path(imid, sub_dir)
  
  #load OT genes
  
  if(imid=="SLE"){
    
    ot.genes<-ot.sle
    
  }else if(imid=="CD"){
    
    ot.genes<-ot.cd
    
    
  } else if(imid=="UC"){
    
    ot.genes<-ot.uc
    
    
  }  else if(imid=="PS"){
    
    ot.genes<-ot.ps
    
    
  }  else if(imid=="RA"){
    
    ot.genes<-ot.ra
    
    
  }else{
    
    ot.genes<-ot.psa
    
    
  }
  
  
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  } else {
    print("Dir already exists!")
  }
  genes.sel.hi<-c(var.sel.hi,ot.genes)
  genes.sel.lo<-c(var.sel.lo,ot.genes)
  
  class.pipeline(hi.adj,lo.adj,A.set1,A.set2,imid,genes.sel.hi,genes.sel.lo)
  
}else if (my.options[i]=="DEG_VAR"){
  #################################################################
  message ("DEG_VAR pipeline for class Discovery ",imid)
  
  sub_dir <-paste0(my.options[i],"_",my.selection)
  output_dir <- file.path(imid, sub_dir)
  
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  } else {
    print("Dir already exists!")
  }
  
  genes.sel.hi<-c(rownames(res.hi),var.sel.hi)
  genes.sel.lo<-c(rownames(res.lo),var.sel.lo)
  
  class.pipeline(hi.adj,lo.adj,A.set1,A.set2,imid,genes.sel.hi,genes.sel.lo)
  
 }else{
      
      #################################################################
      ##option 2,most variable genes vs CTRL
      message ("VAR pipeline for class Discovery ",imid)
      
   sub_dir <-paste0(my.options[i],"_",my.selection)
   output_dir <- file.path(imid, sub_dir)
      
   if (!dir.exists(output_dir)){
     dir.create(output_dir)
   } else {
     print("Dir already exists!")
   }
   
   genes.sel.hi<-c(var.sel.hi)
   genes.sel.lo<-c(var.sel.lo)
   
   class.pipeline(hi.adj,lo.adj,A.set1,A.set2,imid,genes.sel.hi,genes.sel.lo)
  
    }

}

  }



