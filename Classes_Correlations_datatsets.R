library(doParallel)
library(dplyr)
library(edgeR)
library(limma)
library(ggplot2)

setwd("/Users/agomez/IMIDOMICS/CLASS/")
####load counts

sampleTable<-read.csv("Final_Table_Curated_Sex_FINAL.csv")
rownames(sampleTable)<-sampleTable$celgIdAll

##load cts
load("Raw_Counts_RNAseq_OK.RData")


#load original, our data and test genes in each class to get the most differentially expressed if they are in the rank

##load data

##load pheno


all_cts<-all_cts[,colnames(all_cts) %in% rownames(sampleTable)]
sampleTable<-sampleTable[rownames(sampleTable) %in% colnames(all_cts),]

all_cts<-all_cts[,match(rownames(sampleTable),colnames(all_cts))]

##select only disease

A <- as.character(sampleTable[sampleTable$imid == "SLE","celgIdAll"])

samples2test<-sampleTable[rownames(sampleTable) %in% A,]

##Build DG object
cts<-all_cts[,colnames(all_cts) %in% A]

##load classes as phenodata

classes<-read.csv("Classes_SLE.csv")

rownames(classes)<-classes$ID

samples2test<-merge(samples2test,classes, by="row.names")

##put cts same order as phenoclasses

cts<-cts[,match(samples2test$ID,colnames(cts))]

group<-samples2test$consensuscluster
group<-droplevels.factor(group)


y <- DGEList(cts,group=group)

y$samples$lib.size <- colSums(y$counts)

##filter genes with low counts

keep <- rowSums(cpm(y)>1) >= 150#### change group depending replicates
y <- y[keep, ]

dim(y)
##Normalize
y <- calcNormFactors(y, method="TMM")

##Test for each class against others

samples2test$Sex<-ifelse(samples2test$sex=="FEMALE","F","M")
##BATCH

my.batch<-paste("B", samples2test$new.plate,sep="")

##add extra column when compare each class: i.e 1 is A and rest is B

samples2test$COMP1<-ifelse(samples2test$consensuscluster=="1","B","A")
samples2test$COMP2<-ifelse(samples2test$consensuscluster=="2","B","A")
samples2test$COMP3<-ifelse(samples2test$consensuscluster=="3","B","A")###more than 4000
samples2test$COMP4<-ifelse(samples2test$consensuscluster=="4","B","A")
samples2test$COMP5<-ifelse(samples2test$consensuscluster=="5","B","A")

design <- model.matrix(~samples2test$Sex + my.batch + samples2test$age+ samples2test$COMP5)

colnames(design)[16]<-"Comp"
rownames(design)<-rownames(samples2test)

###add design to the object
y2 <- estimateDisp(y, design)
#Now estimate gene-specific dispersions:

plotBCV(y2)

####radical solution

fit <- glmQLFit(y2, design,robust=T)
res <- glmQLFTest(fit, coef="Comp")##or  contrast=B.LvsP

results<-topTags(res,n= dim(res$table)[1],adjust.method="bonferroni", sort.by="PValue")
results<-results$table

####establish high significance threshold to select genes
res.s<-results[results$FWER<.001,]

##build the rank


class1.genes<-as.character(rownames(res.s))
class2.genes<-as.character(rownames(res.s))
class3.genes<-as.character(rownames(res.s))
class4.genes<-as.character(rownames(res.s))
class5.genes<-as.character(rownames(res.s))

##select the all the genes

genes.list <- unique(c(class1.genes,class2.genes,class3.genes,class4.genes,class5.genes))
save(genes.list,file = "Gene_DEG_5classes.RData")
######TEST

###First find classes in test.data

###RNAseq Science
data.to.test<-read.delim("GSE72509_SLE_RPKMs.csv")
data.to.test<-data.to.test[!(duplicated(data.to.test$SYMBOL)),]
rownames(data.to.test)<-data.to.test$SYMBOL
data.to.test<-data.to.test[,-1]


##get only SLE

##RNAseq
data.to.test<-data.to.test[,grepl("SLE",colnames(data.to.test))]




##filter those ones with low number of counts
keep <- rowSums(data.to.test>1) >= 36#### change group depending replicates
data.to.test <- data.to.test[keep, ]

##??elect 2510 genes

data.to.test<-data.to.test[rownames(data.to.test) %in% genes.list,]
data.to.test<-log2(data.to.test+1)

dat <-data.to.test

#remvoe outlier S02

PCA1 <- M3C::pca(dat)

#dat<-dat[,-1]
#dat<-dat[!(rownames(dat) %in% "MIR342"),]
###Heatmap

library(gplots)
library(RColorBrewer)

#Use M3C approach

library(M3C)
library(NMF) # loading for aheatmap plotting function
library(ggsci) # more cool colours
#we run the algorithm using the default settings (100x monte carlo iterations and 100x inner replications).
#We have found the results generally stable using these parameters, 
#although they may be increased or the algorithm may simply be run a few extra times to test stability. 
#Plots from the tool and an .csv file with the numerical outputs may be printed into the working directory 
#(by adding printres = TRUE). We will set the seed in this example, incase you wish to repeat our results exactly (seed = 123). 
#We will add the annotation file for streamlined downstream analyses (des = annot).


res.m3c <- M3C(dat, cores=7,seed = 123, removeplots = F)
stats<-res.m3c$scores
write.csv(stats,"STATS_SLE_GSE72509_5classes_DEG.csv")
# get the data out of the results list (by using $ - dollar sign), use 2 clusters (see RCSI plot)
data <- res.m3c$realdataresults[[5]]$ordered_data # this is the data
annon <- res.m3c$realdataresults[[5]]$ordered_annotation # this is the annotation
ccmatrix <- res.m3c$realdataresults[[5]]$consensus_matrix # this is the consensus matrix


# normalise and scale the data
data <- t(scale(t(data))) # z-score normalise each row (feature)
data <- apply(data, 2, function(x) ifelse(x > 4, 4, x)) # compress data within range
data <- apply(data, 2, function(x) ifelse(x < -4, -4, x)) # compress data within range


library(heatmap.2x)## approach

library(RColorBrewer)
#samples<-ifelse(annon$cond=="RA_HI","red","black")
#annon$cons<-ifelse(annon$consensuscluster=="1","1","2")
#annon$cons<-as.factor(annon$cons)

cons<-ggsci::pal_futurama()(max(levels(annon$consensuscluster)))[as.factor(annon$consensuscluster)]#5 clusters
#cons<-ggsci::pal_futurama()(max(levels(annon$cons)))[as.factor(annon$cons)]#5 clusters

spcol<-cons
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)


###order data

data<-data[complete.cases(data),]

hp<-heatmap.2x(data, col=rev(cols), Colv = NA,Rowv=TRUE, scale="none", ColSideColors=spcol,
               trace="none", dendrogram="both", 
               cexRow=.1, cexCol=1.4,
               main="SLE Classes GSE72509 5 classes DEG ",
               labCol=NA,
               labRow=NA, 
               density.info="none",
               hclust=function(x) hclust(x,method="complete"),
               distfun=function(x) as.dist((1-cor(t(x)))/2)
               
)


#####reorder data according heatmap


annon$ID<-rownames(annon)

write.csv(annon,"SLE_classes_GSE72509_5_classes_DEG.csv")

####correspondence

##find the correlation between our classes and GSE classes: Hint:correlate the same metric genes


group<-ifelse(annon$consensuscluster==3, "B","A")

# Specify a design matrix without an intercept term
design <- model.matrix(~ 0 + group)
rownames(design)<-rownames(annon)

##get the complete matrix again


data.to.test<-data.to.test[,match(rownames(annon),colnames(data.to.test))]

###add design to the object
fit <- lmFit(data.to.test,design = design)


cont.matrix <- makeContrasts(Comp=groupB - groupA,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)

res.gse<-topTable(fit.cont,coef=1,dim(data.to.test)[1],sort.by="p")

res.gse$fcSign<-sign(res.gse$logFC)
res.gse$logP<-(-log10(res.gse$adj.P.Val))
res.gse$metric<-res.gse$logP/res.gse$fcSign
res.gse<-res.gse[order(-res.gse$metric),]

#gse.sel.1<-res.gse
#gse.sel.2<-res.gse
#gse.sel.3<-res.gse

###select same in GSE

gse.sel.1<-res.gse[rownames(res.gse) %in% genes.list,]

gse.sel.2<-res.gse[rownames(res.gse) %in% genes.list,]

gse.sel.3<-res.gse[rownames(res.gse) %in% genes.list,]

gse.sel.4<-res.gse[rownames(res.gse) %in% genes.list,]

gse.sel.5<-res.gse[rownames(res.gse) %in% genes.list,]

save(gse.sel.1,gse.sel.2,gse.sel.3,gse.sel.4,gse.sel.5,file="Ranks_LAST_DEG_GSE72509.RData")

##correlate

##correlations


##try automatically in lists

##correlate

##correlations
##chnage depending what you want to correlate

gse.sel<-list(gse.sel.1,gse.sel.2,gse.sel.3,gse.sel.4,gse.sel.5)
gse.sel<-list(gse.sel.1,gse.sel.2,gse.sel.3)##Carlucci's data

#load("Ranks_original_classes.RData")
ori<-list(class1.rank,class2.rank,class3.rank,class4.rank,class5.rank)
save(gse.sel,file="Ranks_GSE72509.RData")
save(gse.sel,file="Ranks_GSE110865.RData")

results.corrs<-data.frame()
out<-data.frame()


for (j in 1:length(ori)){
  #  print(j)
  orig.test<-ori[[j]]
  
  for (i in 1:length(gse.sel)){ 
    #print(i)
    
    to.test<-gse.sel[[i]][match(rownames(ori[[j]]),rownames(gse.sel[[i]])),]
   # orig.test<-ori[[j]]
    
    r2 <- cor.test(as.numeric(to.test$metric),as.numeric(orig.test$metric),method="spearman",exact = F)$estimate
    r2<-signif(r2, digits=3)
    my.p<-cor.test(as.numeric(to.test$metric),as.numeric(orig.test$metric),method="spearman",exact = F)$p.value
    my.p<-signif(my.p, digits=3)
    
    
    results.corrs[j,1]<-j
    results.corrs[j,2]<-r2
    results.corrs[j,3]<-my.p
    results.corrs[j,4]<-i
    
    
    out<-rbind(out,results.corrs)
    out<-unique(out,MARGIN=1)

    
  }
}

colnames(out)<-c("Ori_class","Corr","pval","External_data_class")
write.csv(out,"Correlations_GSE72509_Ori_classes.csv")


##chnage depending what you want to correlate

#gse.sel<-gse.sel.2
#ori<-class5.rank

#ori<-ori[match(rownames(gse.sel),rownames(ori)),]

#GSE72509_Class1_Class1

#mod1 <-lm(as.numeric(gse.sel$metric)~as.numeric(ori$metric))## 
#modsum<-summary(mod1)

#r2 <- cor.test(as.numeric(gse.sel$metric),as.numeric(ori$metric),method="spearman",exact = F)$estimate
#my.p<-cor.test(as.numeric(gse.sel$metric),as.numeric(ori$metric),method="spearman",exact = F)$p.value
#my.p<-signif(my.p, digits=3)


#mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))

#plot(as.numeric(gse.sel$metric)~as.numeric(ori$metric), main=paste("Correlation","GSE72509 Class 2 vs Class2", "p-value", my.p, sep=" "), xlab="GSE72509", ylab="SLE Data", pch=20,col="grey40")

#abline(mod1, col="red")

#legend('topleft', legend = mylabel, bty = 'n')


