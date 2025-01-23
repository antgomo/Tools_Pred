library(TCA)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(doParallel)

library(minfi)
library(DMRcate)
library(data.table)

setwd("/mnt/md127/Analysis_Projects/FIS_Meth/P6_validation/")
####load pheno
samplesheet<-read.csv("samplesheets_p6_notinFIS.csv")
colnames(samplesheet)<-gsub("\\.x","",colnames(samplesheet))
###moderates are non responders
samplesheet<-samplesheet[!(samplesheet$response=="MOD"),]##not moderates
samplesheet$RESP<-ifelse(samplesheet$response=="GOOD","R","NR")



##select only baseline
samplesheet<-samplesheet[samplesheet$week=="wk0",]
samplesheet$SEX<-ifelse(samplesheet$Sex=="FEMALE","F","M")
######################load TCA data

#######################
load("/mnt/md127/Analysis_Projects/FIS_Meth/TCA_results/Validation/P6/Betas_cell_subtype_P6_RNR.RData")

# estimate cell-type-specific methylation for the associated site


####[1] "CD8T"  "CD4T"  "NK"    "Bcell" "Mono"  "Neu"

###Interactome hotspots combined with methylation
setwd("/mnt/md127/Analysis_Projects/FIS_Meth/TCA_results/")

library(FEM)
load("../hprdAsigH-13Jun12.Rd")
#####GO enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)
#library(GOstats)
library(GO.db)
library(DOSE)
library(annotate)



analysis<-"P6_antiTNF"

if (!dir.exists(analysis)){
  dir.create(analysis)
} else {
  print("Dir already exists!")
}

names(Z_hat)<-c("CD8T","CD4T","NK","Bcell","Mono","Neu")

for( r in 1:length(Z_hat)){

##with adjusted betas

###arrange the order!! avoid M
  

    sub_dir <-names(Z_hat)[r]
    
    output_dir <- file.path(analysis,sub_dir)
    if (!dir.exists(output_dir)){
      dir.create(output_dir)
    } else {
      print("Dir already exists!")
    }
    


  beta.test<-Z_hat[[r]][,colnames(Z_hat[[r]]) %in% samplesheet$Name]
  beta.test<-beta.test[,match(samplesheet$Name,colnames(beta.test))]


##you must adjust the beta matrix before deal with analyasis!!!
  #design <- model.matrix(~0 + wk0$S0_DAS28 +wk0$Age)##w0
  design <- model.matrix(~0 + samplesheet$RESP +samplesheet$Age +samplesheet$SEX)##w0
  
####generate adjusted beta matrix

  fit <- lmFit(beta.test, design)
  mAdj.fit    <- fit$coefficients[,-1]

  mAdj<- as.matrix(beta.test - mAdj.fit %*% t(design[,-1]))

##for betas
  betaAdj     <- ilogit2(mAdj)

##for Ms

  rm(fit)
  rm(mAdj)
  rm(mAdj.fit)
  rm(beta.test)
  gc()


##check for covariates

##try to find batches in data
##same without batch
#library(ggplot2)

#pca <- prcomp(t(betaAdj))
#pca_data_perc<-round(100*((pca$sdev)^2/sum(pca$sdev^2)),digits=2)



#df_pca_data<-data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], batch=wk0$SEX)

#ggplot(df_pca_data, aes(PC1,PC2, color = batch))+
 # geom_point(size=3)+
  #ggtitle(colnames(cell_comp2)[6])+
  #labs(x=paste0("PC1 (",pca_data_perc[1],")"), y=paste0("PC2 (",pca_data_perc[2],")"))

 # source("../GenStatM_quant.R")
  
  #statM.o <- GenstatM_quant(betaAdj,wk0$S0_DAS28,"EPIC")
  statM.o <- GenStatM(betaAdj,samplesheet$RESP,"EPIC")
  
 # intEpi.o<-DoIntEpi450k(statM.o,hprdAsigH.m,c=2)##quant c=2
  intEpi.o<-DoIntEpi450k(statM.o,hprdAsigH.m,c=1)##
  
  save(intEpi.o,file=paste0(output_dir,"/","int_stats_",analysis,"_",names(Z_hat)[r],".RData"))

  source("../DOEpiMod2.R")##get nominal p-vals

  EpiMod.o<-DoEpiMod2(intEpi.o,nseeds=100,gamma=0.5,nMC=1000,sizeR.v=c(1,100), minsizeOUT=10,writeOUT=F,nameSTUDY=analysis,ew.v=NULL);

  save(EpiMod.o,file=paste0(output_dir,"/","Epimod_","_",names(Z_hat)[r],".RData"))
  rm(betaAdj)
  gc()
  
##plots


library(marray)
library(corrplot)


for(m in 1:length(names(EpiMod.o$topmod))){
  FemModShow(EpiMod.o$topmod[[m]],
             name=paste0(output_dir,"/",names(EpiMod.o$topmod)[m]),EpiMod.o,mode="Epi")}


##GO enrichment for each module



#plot.list<-list()

  for(m in 1:length(names(EpiMod.o$topmod))){
  
  
    genesid<-as.character(EpiMod.o$topmod[[m]][,"Symbol"])
    genesid<-genesid[ genesid != "" ]##remove empty elements from a vector
  
    genesid <- genesid[!is.na(genesid)]
  
    eg<-bitr(genesid, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
    ego2 <- enrichGO(gene         = eg$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.1,
                   qvalueCutoff  = 0.05,
                   readable=T
    )
  
    png(paste0(output_dir,"/","GO_Module","_",analysis,"_",names(Z_hat)[r],"_",names(EpiMod.o$topmod)[m],".png",sep=""),width = 1200, height = 800)
  
      p<-dotplot(ego2, title=names(EpiMod.o$topmod)[m],showCategory=35)
  
      print(p)
  
      dev.off()
  # myplots[[i]] <- p
  }
#print(myplots)
}






