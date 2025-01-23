library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(FlowSorted.Blood.EPIC)
library(minfi)
library(DMRcate)

################################################

##load samplesheet

setwd("/mnt/md127/Analysis_Projects/FIS_Meth/P6_validation/")

samplesheet<-read.csv("samplesheet_P6_validation.csv")
###moderates are non responders

samplesheet$RESP<-ifelse(samplesheet$response=="GOOD","R","NR")
load("P6_validation.RData")

betas<-getBeta(countsEPIC$normalizedData)

betas<-betas[,colnames(betas) %in% samplesheet$Name]
betas<-betas[,match(samplesheet$Name,colnames(betas))]


cell_comp<-countsEPIC$counts
cell_comp<-cell_comp[match(samplesheet$Name,rownames(cell_comp)),]

betas<-rmSNPandCH(betas, dist = 2, mafcut = 0.05, and = TRUE, rmcrosshyb = TRUE, rmXY=TRUE)
rm(countsEPIC)
##validate Predictor

#build design matrix
library(limma)

design <- model.matrix(~samplesheet$RESP+samplesheet$week+samplesheet$Age + cell_comp)##only females

fit <- lmFit(betas, design)
mAdj.fit    <- fit$coefficients[,-c(1:2)]

mAdj<- as.matrix(betas - mAdj.fit %*% t(design[,-c(1:2)]))

betaAdj     <- ilogit2(mAdj)

rm(fit)
rm(mAdj)
rm(mAdj.fit)

gc()

targets$grptime <- paste(targets$RESP,targets$Week,sep="_")

#design <- model.matrix(~0 + targets$SEX + targets$Age + cell_comp + targets$grptime)##
cell_comp<-cell_comp[rownames(cell_comp) %in% targets$xipName,]
cell_comp<-cell_comp[match(targets$xipName,rownames(cell_comp)),]


design <- model.matrix(~0 + targets$grptime + targets$SEX + targets$Age + cell_comp)##

colnames(design)<-gsub("targets\\$","",colnames(design))

corfit<-duplicateCorrelation(betas,design,block=targets$Code)##0.2665931

fit<-lmFit(betas,design,block=targets$Code,correlation=corfit$consensus.correlation)#0.264631 without outliers
fit<-lmFit(betas,design,block=targets$Code,correlation=0.2665931)

#Now you can extract the differences you like, using appropriate contrast matrix. For example, which genes changed over time in R

##change contrast.matrix , i.e comment uncomment depending on what you want to get

cont.matrix <- makeContrasts(
  
  T0vsT1inR=grptimeR_1 - grptimeR_0,##effect time R
  # T0vsT1inNR=grptimeNR_1 - grptimeNR_0, ##effect time NR
  #  Diff=(grptimeR_1 - grptimeR_0)-(grptimeNR_1 - grptimeNR_0),#3 Rvs NR across all time
  levels=design)


fit2<- contrasts.fit(fit,cont.matrix)

fit2 <- eBayes(fit2)


results <- topTable(fit2,num=dim(fit2)[1],sort.by="P",adjust.method = "BH") ####Comp

res<-subset(results,results$P.Value<.01)


dmp <- merge(res,ann850k[,c("UCSC_RefGene_Name","Relation_to_Island","UCSC_RefGene_Group")], by="row.names")
rownames(dmp)<-dmp$Row.names
dmp<-dmp[,-1]

###remove non annotated cpgs

#dmp<-dmp[!(is.na(dmp$UCSC_RefGene_Name) | dmp$UCSC_RefGene_Name==""), ]


write.csv(dmp,"without_outliers/Cpg_diff_meth_R_NR_anytimepoint.csv")



###Heatmap

annot<-targets[,c(1,36,37)]#w0

colnames(annot)[1]<-"ID"
annot$ID<-gsub("-","",annot$ID)
annot$IF<-paste(annot$ID,annot$grptime,sep="_")


fit <- lmFit(betas, design)
mAdj.fit    <- fit$coefficients[,-c(1:4)]

mAdj<- as.matrix(betas - mAdj.fit %*% t(design[,-c(1:4)]))

betaAdj     <- ilogit2(mAdj)

rm(fit)
rm(mAdj)
rm(mAdj.fit)

gc()

###change colnames betaAdj

colnames(betaAdj)<-paste(annot$ID,annot$grptime,sep="_")

###heatmap anytime point

##order by adequate way
Rw0<-as.character(annot[annot$grptime=="R_0","IF"])
Rw12<-as.character(annot[annot$grptime=="R_1","IF"])

NRw0<-as.character(annot[annot$grptime=="NR_0","IF"])
NRw12<-as.character(annot[annot$grptime=="NR_1","IF"])


##get the betas

my.data<-as.data.frame(betaAdj[rownames(betaAdj) %in% rownames(res) ,])#for R

##change column names
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)


# normalise and scale the data
data <- t(scale(t(my.data))) # z-score normalise each row (feature)
data <- apply(data, 2, function(x) ifelse(x > 4, 4, x)) # compress data within range
data <- apply(data, 2, function(x) ifelse(x < -4, -4, x)) # compress data within range


##Try to get the genes behaving in some way or another in R and NR
##better use log values not scaled
###R

#######select sample names

data<-as.data.frame(data)

data$uRW0<-rowMeans(data[,Rw0],na.rm=T)
data$uRW12<-rowMeans(data[,Rw12],na.rm=T)

data$uNRW0<-rowMeans(data[,NRw0],na.rm=T)
data$uNRW12<-rowMeans(data[,NRw12],na.rm=T)


toclust<-data[,c(171,173,172,174)]##with  data +4/-3
# "uRW0"   "uNRW0"  "uRW12"  "uNRW12"



##now define classes

###UP in 0 and DOWN in 12

UP0D12<-toclust[(toclust[,1]>toclust[,2] & toclust[,3]>toclust[,4] ),]##Up in all resp
D0D12<-toclust[(toclust[,1]<toclust[,2] & toclust[,3]<toclust[,4]) ,]##Up in all NR

##########plot heatmap

library(ComplexHeatmap)
###first, check number of real clusters

library(RColorBrewer)

data$real_clust<-ifelse(rownames(data) %in% rownames(UP0D12),"UP_RESP",ifelse(rownames(data) %in% rownames(D0D12),"DOWN_RESP","TRANS"))

cluster_real<-data$real_clust
names(cluster_real)<-rownames(data)

ht<-Heatmap(data[,c(171,173,172,174)], split = cluster_real, row_names_gp = gpar(fontsize = 6)
            ,cluster_columns = FALSE,show_row_names = F,
            
            heatmap_legend_param = list(title = "CpG changes"))

##k means split


ht<-make_row_cluster(ht)
row_order<-unlist(ht@row_order_list)

png("Heatmap_changes_anytime_point.png",width = 1200, height = 800)

ht

draw(ht, column_title = "CpG R vs NR anytime")


dev.off()










my.data<-as.data.frame(betaAdj[rownames(betaAdj) %in% rownames(res) ,colnames(betaAdj) %in% annot[annot$RESP=="R","IF"]])#for R
my.data<-as.data.frame(betaAdj[rownames(betaAdj) %in% rownames(res) ,colnames(betaAdj) %in% annot[annot$RESP=="NR","IF"]])#for NR

annot<-annot[annot$RESP=="R",]
annot<-annot[annot$RESP=="NR",]


annot$IF<-gsub("_","",annot$IF)
rownames(annot)<-annot$IF
#I used the coefficient of variation to remove more less variable genes

CV <- apply(my.data,1,sd)/(rowMeans(my.data))
names <- names(CV)[CV>mean(CV)] # 0.175, 0.17

##same order as annot
data3 <- subset(my.data, row.names(my.data) %in% names)
#data3<-my.data
colnames(data3)<-gsub("_","",colnames(data3))
##avoid problems with annotation

colnames(annot)<-c("Code","RESP","grptime","ID")

res.m3c <- M3C(data3, cores=7,des=annot,seed = 123, removeplots = F)

# get the data out of the results list (by using $ - dollar sign), use 2 clusters (see RCSI plot)
data <- res.m3c$realdataresults[[2]]$ordered_data # this is the data
annon <- res.m3c$realdataresults[[2]]$ordered_annotation # this is the annotation
ccmatrix <- res.m3c$realdataresults[[2]]$consensus_matrix # this is the consensus matrix


# normalise and scale the data
data <- t(scale(t(data))) # z-score normalise each row (feature)
data <- apply(data, 2, function(x) ifelse(x > 4, 4, x)) # compress data within range
data <- apply(data, 2, function(x) ifelse(x < -4, -4, x)) # compress data within range



library(heatmap.2x)## approach


samples<-ifelse(annon$grptime=="NR_1","red","blue1")

cons<-ggsci::pal_futurama()(max(levels(annon$consensuscluster)))[as.factor(annon$consensuscluster)]#5 clusters


spcol <- rbind(cons, samples)

cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)


heatmap.2x(data, col=rev(cols), Colv = NA,Rowv=TRUE, scale="none", ColSideColors=spcol,
           trace="none", dendrogram="both", 
           cexRow=1, cexCol=1.4,
           main="NR Time effect",
           labCol=NA,labRow=NA, density.info="none",
           hclust=function(x) hclust(x,method="complete"),
           distfun=function(x) as.dist((1-cor(t(x)))/2)
           
)

legend("topright",c("W0","W12"),pch=20:20,col=c("blue1","red"))



##GO enrichment
#####GO enrichment analysis

library(ggplot2)

#####GO enrichment analysis
library(methylGSA)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

###enrichment

cpgval<-dmp$P.Value
names(cpgval)<-rownames(dmp)

res1 <- methylglm(cpg.pval = cpgval, minsize = 10, maxsize = 1000, GS.type = "GO",array.type = "EPIC")

res1<-res1[res1$pvalue<.05,]

write.csv(res1,paste0("without_outliers/GO_R_NR_.","EffecTime",".csv",sep=""),row.names=T) 


###Forest_plot GOmeth
fp<-ggplot(data=gometh.res[1:55,],aes(x=GeneRatio,y=Term))+
  geom_point(aes(size=DE,fill=P.DE), colour="black",shape=21)+
  theme(axis.text.y = element_text(margin=margin(.1,.5,.8,.5,"pt"),size=7))+
  ggtitle("GO R vs NR Anytime point")+
  ylab(NULL) +
  theme(plot.margin=unit(c(.01,1,.1,1.2),"cm"))


write.csv(gometh.res,"GO_UP_NR_unfiltered.csv")
###DMR


library(DMRcate)

##use Ms

myMs.noSNPs <- logit2(betas)


nrow(myMs.noSNPs)


#design <- model.matrix(~targets$grptime + targets$SEX + targets$Age + cell_comp)##

#colnames(design)<-gsub("targets\\$","",colnames(design))
myannotation <- cpg.annotate("array", myMs.noSNPs, what="M", arraytype = "EPIC",analysis.type="differential", design=design,contrasts = TRUE, cont.matrix = cont.matrix,
                             coef="T0vsT1inR",fdr = 0.05)## 3 for time in R


##use DMP info from previous part

pos.dmps<-which(myannotation$ID %in%  rownames(dmp))

myannotation$is.sig[pos.dmps]<-"TRUE"
myannotation$is.sig[-pos.dmps]<-"FALSE"

myannotation$is.sig<-as.logical(myannotation$is.sig)

pvals.dmp<-dmp$P.Value

myannotation$indfdr[pos.dmps]<-pvals.dmp
myannotation$indfdr[-pos.dmps]<-0.49

## 2 for time in NR
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2,pcutoff = 0.05,p.adjust.method="none")

results.ranges <- extractRanges(dmrcoutput, genome = "hg38")

results.dmr<-as.data.frame(results.ranges)
results.dmr<-results.dmr[results.dmr$minfdr<.05,]

#results.dmr<-results.dmr[!(is.na(results.dmr$overlapping.promoters)),]

write.csv(results.dmr,"DMR_R_w0vsw12_fixed.csv")


for(m in 1:length(results.ranges)){
  
  groups <- c(R="blue1", NR="red")
  
  cols <- groups[as.character(phenoData3$RESP)]
  samps <- c(grep("^R",names(cols))[1:5], grep("^NR",names(cols))[1:5])
  
  
  png(paste0("Results_FIS_May_PACTABA_approach/DMR","_",m,"_",unique(unlist(strsplit(results.ranges[m]$overlapping.promoters,"\\-|\\,")))[1],".png",sep=""),width = 800, height = 800)
  
  p<-DMR.plot(ranges=results.ranges, dmr=m, CpGs=myMs.noSNPs, what="M", arraytype = "EPIC",phen.col=cols, genome="hg38",samps = samps)
  
  print(p)
  
  dev.off()
  
}

###Interactome hotspots combined with methylation


library(FEM)
data(Realdata)

##subset by samples

##use betaAdj because R is in the colnames and it is in the same orders as betas

betas.fis<-betaAdj[,grepl("_R_",colnames(betaAdj),perl=T)]
betas.fis<-betaAdj[,grepl("_NR_",colnames(betaAdj),perl=T)]


statM.o <- GenStatM(betas.fis,targets[grepl("^R_",targets$grptime,perl=T),"grptime"],"EPIC")##time effect
statM.o <- GenStatM(betas.fis,targets[grepl("^NR_",targets$grptime,perl=T),"grptime"],"EPIC")##time effect
statM.o <- GenStatM(betaAdj,targets$grptime,"EPIC")##anytime point

##with adjusted betas

intEpi.o<-DoIntEpi450k(statM.o,Realdata$adjacency,c=1)

source("DOEpiMod2.R")##get nominal p-vals

EpiMod.o<-DoEpiMod2(intEpi.o,nseeds=100,gamma=0.5,nMC=1000,sizeR.v=c(1,100), minsizeOUT=10,writeOUT=TRUE,nameSTUDY="FIS_Long",ew.v=NULL);

save(EpiMod.o2,file="Results_RESP_NRESP_fix/Epimod_NR_Time_fixed.RData")

EpiMod.o2<-EpiMod.o
load("Epimod_Long_Time_anytime.RData")



##plots

library(marray)
library(corrplot)

for(m in 1:length(names(EpiMod.o2$topmod))){
  FemModShow(EpiMod.o2$topmod[[m]],
             name=names(EpiMod.o2$topmod)[m],EpiMod.o2,mode="Epi")}


##GO enrichment for each module

#####GO enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)
#library(GOstats)
library(GO.db)
library(DOSE)
library(annotate)


#plot.list<-list()

for(m in 1:length(names(EpiMod.o2$topmod))){
  
  
  genesid<-as.character(EpiMod.o2$topmod[[m]][,"Symbol"])
  
  genesid<-genesid[ genesid != "" ]##remove empty elements from a vector
  genesid <- genesid[!is.na(genesid)]
  
  eg<-bitr(genesid, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  ego2 <- enrichGO(gene         = eg$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable=T
  )
  
  png(paste0("Results_RESP_NRESP_fix/GO_Module","_","Long_anytime","_",names(EpiMod.o$topmod)[m],".png",sep=""),width = 1200, height = 800)
  
  p<-dotplot(ego2, title=names(EpiMod.o$topmod)[m],showCategory=35)
  
  print(p)
  
  dev.off()
  # myplots[[i]] <- p
}
#print(myplots)

####

##get HCK betas

HCK<-"cg20547606"  ##all
