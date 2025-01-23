library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(doParallel)

library(minfi)
library(DMRcate)
setwd("/mnt/md127/Analysis_Projects/FIS_Meth/FIS_Dat/")
####load pheno
load("phenoData.RData")


idat.folder <- "/mnt/md127/Analysis_Projects/FIS_Meth/FIS_Dat/idats/" ###change for wour path dir

targets <- read.metharray.sheet(base=idat.folder)
targets$Basename<-paste(idat.folder,targets$xipName, sep="")

load("FIS_norm_cellcomp_data.RData")
###begin here after normalized beta pipeline is done

setwd("/mnt/md127/Analysis_Projects/FIS_Meth/")


##only betas in blood target

targets<-targets[targets$tissue=="blood",]

betas <-betas[,colnames(betas) %in% targets$xipName]
targets <-targets[targets$xipName %in% colnames(betas),]

betas<-betas[,match(targets$xipName,colnames(betas))]

pD.w0<-targets[grep("S0",targets$donation,perl=T),]
pD.w12<-targets[grep("S12|SRM",targets$donation,perl=T),]

##remove sexual chromosomes and SNPs probes

betas<-rmSNPandCH(betas, dist = 2, mafcut = 0.05, and = TRUE, rmcrosshyb = TRUE, rmXY=TRUE)



##select betas depending week

betas.fis<-betas[,colnames(betas) %in% pD.w12$xipName]

##put in same order as phenodata to get clear Sample names
pD.w12$Code<-gsub("-S12|-SRM","",pD.w12$donation)
#pD.w0$Code<-gsub("-S0","",pD.w0$donation)

betas.fis<-betas.fis[,match(as.character(pD.w12$xipName),colnames(betas.fis))]#w12
#betas.fis<-betas.fis[,match(as.character(pD.w0$xipName),colnames(betas.fis))]#w0

colnames(betas.fis)<-pD.w12$Code
#colnames(betas.fis)<-pD.w0$Code


phenoData<-read.csv("phenodata_FIS_fixed.csv")
phenoData$delta<-phenoData$S0_DAS28-phenoData$SR12_DAS28
phenoData$EULAR2<-ifelse(phenoData$SR12_DAS28<3.2 & phenoData$delta>1.2,"GOOD",
                        ifelse(phenoData$SR12_DAS28<3.2 & phenoData$delta>0.6 & phenoData$delta<1.2 ,"MOD",
                               ifelse(phenoData$SR12_DAS28<3.2 & phenoData$delta<0.6,"NON_RESP",
                                      ifelse(phenoData$SR12_DAS28>3.2 & phenoData$SR12_DAS28<5.1 & phenoData$delta>1.2,"MOD",
                                             ifelse(phenoData$SR12_DAS28>3.2 &  phenoData$SR12_DAS28<5.1 & phenoData$delta>0.6 & phenoData$delta<1.2,"MOD",
                                                    ifelse(phenoData$SR12_DAS28>3.2 & phenoData$SR12_DAS28<5.1 & phenoData$delta<0.6,"NON_RESP",
                                                           ifelse(phenoData$SR12_DAS28>5.1 & phenoData$delta>1.2,"MOD",
                                                                  ifelse(phenoData$SR12_DAS28> 5.1 & phenoData$delta>0.6 & phenoData$delta<1.2 ,"NON_RESP","NON_RESP"))))))))

phenoData$RESP<-ifelse(phenoData$EULAR2=="NON_RESP","NR",ifelse(phenoData$EULAR2=="NA","NA","R"))



library(limma)

##simple limma approach, build design matrix with Cell composition
##avoid Tocilizumab, Tofacinitib and others in treatmen

phenoData3<-phenoData[!(is.na(phenoData$Age)),]
phenoData3<-phenoData3[!(is.na(phenoData3$RESP)),]


#####merge pheno with pD samplesheet
pD.w12<-merge(pD.w12,phenoData3[,c("antiTNF","RESP","SEX","Age","Code","EULAR2")],by="Code")
pD.w12<-pD.w12[!(pD.w12$antiTNF=="Tocilizumab" | pD.w12$antiTNF=="Tofacitinib" | pD.w12$antiTNF=="Abatacept" | pD.w12$antiTNF=="Other"),]


pD.w12<-pD.w12[pD.w12$Code %in% colnames(betas.fis),]


betas.fis<-betas.fis[,colnames(betas.fis) %in% pD.w12$Code]

betas.fis<-betas.fis[,match(pD.w12$Code,colnames(betas.fis))]

cell_comp2<-cell_comp[rownames(cell_comp) %in% pD.w12$xipName,]
cell_comp2<-cell_comp2[match(pD.w12$xipName, rownames(cell_comp2)),]

#design <- model.matrix(~pD.w12$RESP  + pD.w12$Age + pD.w12$SEX + cell_comp2)##only sex
design <- model.matrix(~pD.w12$RESP  + pD.w12$Age + cell_comp2)##only sex

colnames(design)[2]<-"Comp" ##Response UP
rownames(design)<-pD.w12$Code

#make desgin without the variable to test and remove everything


fit <- lmFit(betas.fis, design)
mAdj.fit    <- fit$coefficients[,-c(1,2)]

mAdj<- as.matrix(betas.fis - mAdj.fit %*% t(design[,-c(1,2)]))

betaAdj     <- ilogit2(mAdj)

rm(fit)
rm(mAdj)
rm(mAdj.fit)

gc()

###Interactome hotspots combined with methylation


library(FEM)

###load real data, i.e HIPPIE protein adjacency network

load("hprdAsigH-13Jun12.Rd")



##with adjusted betas

###arrange the order!!

beta.test<-betaAdj[,colnames(betaAdj) %in% pD.w12$Code]
beta.test<-beta.test[,match(pD.w12$Code,colnames(beta.test))]

statM.o <- GenStatM(beta.test,pD.w12$RESP,"EPIC")

intEpi.o<-DoIntEpi450k(statM.o,hprdAsigH.m,c=1)

source("DOEpiMod2.R")##get nominal p-vals

EpiMod.o<-DoEpiMod2(intEpi.o,nseeds=100,gamma=0.5,nMC=1000,sizeR.v=c(1,100), minsizeOUT=10,writeOUT=F,nameSTUDY="R_NR_wk12",ew.v=NULL);

save(EpiMod.o,file="wk12_Modules/Epimod_R_NRwk12.RData")
save(intEpi.o,file="wk12_Modules/Stats_Epimod_R_NRwk12.RData")


##plots

 
library(marray)
library(corrplot)


for(m in 1:length(names(EpiMod.o$topmod))){
  FemModShow(EpiMod.o$topmod[[m]],
 name=paste0("wk12_Modules/",names(EpiMod.o$topmod)[m]),EpiMod.o,mode="Epi")}


##GO enrichment for each module

#####GO enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)
#library(GOstats)
library(GO.db)
library(DOSE)
library(annotate)


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

  png(paste0("wk12_Modules/","GO_Module","_","RvsNRwk12","_",names(EpiMod.o$topmod)[m],".png",sep=""),width = 1200, height = 800)
  
  p<-dotplot(ego2, title=names(EpiMod.o$topmod)[m],showCategory=35)
 
   print(p)
  
  dev.off()
 # myplots[[i]] <- p
}
#print(myplots)
