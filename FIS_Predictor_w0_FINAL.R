library(limma)
library(minfi)
library(DMRcate)

setwd("/mnt/md127/Analysis_Projects/FIS_Meth/")
library(gplots)
library(RColorBrewer)
library(M3C)
library(NMF) # loading for heatmap plotting function
library(ggsci)
library(heatmap.2x)## approach
library(randomForest)
library(parallel)
library(ggplot2)
library(rmarkdown)
library(glmnet)
library(doParallel)
library(HandTill2001)
library(pROC)
library(caret)
library(Pomona)
##loading accessory scripts

source("../Classes_RNAseq_2019/Predictor_scripts/multi_brier.R")
source("../Classes_RNAseq_2019/Predictor_scripts/multi_logloss.R")
source("../Classes_RNAseq_2019/Predictor_scripts/RF_utilities.R")
source("../Classes_RNAseq_2019/Predictor_scripts/makefolds.R")





######################################end of function definition
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(doParallel)

setwd("/mnt/md127/Analysis_Projects/FIS_Meth/FIS_Dat/")
####load all samples FIS and P6

idat.folder <- "/mnt/md127/Analysis_Projects/FIS_Meth/FIS_Dat/predictor/" ###change for wour path dir
targets <- read.metharray.sheet(base=idat.folder)

registerDoParallel(cores = 8)

###loading data

#rgset <- read.metharray.exp(targets = targets,force=TRUE)
#library(FlowSorted.Blood.EPIC)

#countsEPIC<-estimateCellCounts2(rgset, compositeCellType = "Blood",
 #                                                                 processMethod = "preprocessNoob",
  #                                                                probeSelect = "IDOL",
   #                                                              cellTypes = c("CD8T", "CD4T", "NK", "Bcell",
       #                                                           "Mono", "Neu"),
    #                                                              referencePlatform =
        #                                                        "IlluminaHumanMethylationEPIC",
     #                                                            referenceset = NULL,
         #                                                        IDOLOptimizedCpGs =IDOLOptimizedCpGs,
        #                                                         returnAll = T)



#save(countsEPIC,targets,file="Predictor_ALLCelltypes.RData")
load("predictor/Predictor_ALLCelltypes.RData")
betas<-getBeta(countsEPIC$normalizedData)
targets<-targets[!(is.na(targets$Age)),]

betas<-betas[,colnames(betas) %in% targets$Name]
betas<-betas[,match(targets$Name,colnames(betas))]


cell_comp<-countsEPIC$counts
cell_comp<-cell_comp[rownames(cell_comp) %in% targets$Name,]

cell_comp<-cell_comp[match(targets$Name,rownames(cell_comp)),]

betas<-rmSNPandCH(betas, dist = 2, mafcut = 0.05, and = TRUE, rmcrosshyb = TRUE, rmXY=TRUE)
rm(countsEPIC)

gc()

targets$SEX<-ifelse(targets$Sex=="FEMALE" | targets$Sex=="F","F","M")


##validate Predictor

#01 build design matrix


design <- model.matrix(~targets$RESP+targets$Age + targets$SEX  + cell_comp)##only wk0



####generate adjusted beta matrix
fit <- lmFit(betas, design)
mAdj.fit    <- fit$coefficients[,-c(1:2)]

mAdj<- as.matrix(betas - mAdj.fit %*% t(design[,-c(1:2)]))

betaAdj     <- ilogit2(mAdj)

rm(fit)
rm(mAdj)
rm(mAdj.fit)
#rm(betas)
gc()

##try to find batches in data
##same without batch
pca <- prcomp(t(betaAdj))
pca_data_perc<-round(100*((pca$sdev)^2/sum(pca$sdev^2)),digits=2)



df_pca_data<-data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], batch=targets$Project)

ggplot(df_pca_data, aes(PC1,PC2, color = batch))+
  geom_point(size=3)+
  labs(x=paste0("PC1 (",pca_data_perc[1],")"), y=paste0("PC2 (",pca_data_perc[2],")"))


##split in test set and training set

training.set<-targets[targets$Project=="FIS",]
test.set<-targets[targets$Project=="P6",]


#################Previous filtering

##get only Promoter + Island
##DMP Annotation
ann850k<-getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann850k<-as.data.frame(ann850k)
dmp <- ann850k[,c("UCSC_RefGene_Name","Relation_to_Island","UCSC_RefGene_Group")]

## Outputs
proms<-dmp[grep("TSS1500|TSS200|5'UTR|1stExon",dmp$UCSC_RefGene_Group),]#

proms.island<-proms[grep("Island",proms$Relation_to_Island),]#


###remove non annotated cpgs

proms.island<-proms.island[!(is.na(proms.island$UCSC_RefGene_Name) | proms.island$UCSC_RefGene_Name==""), ]

##proms

proms<-proms[!(is.na(proms$UCSC_RefGene_Name) | proms$UCSC_RefGene_Name==""), ]
proms<-rownames(proms)

rm(ann850k)
rm(dmp)
gc()

proms.island<-rownames(proms.island)

#get only betas to train

to.train<-betaAdj[,colnames(betaAdj) %in% training.set$Name]##adjust cellcomp
#order them
to.train<-to.train[,match(training.set$Name,colnames(to.train))]

#####calculate sd

#res.sd<-apply(to.train,1,sd)

#sel<-res.sd[res.sd>0.01]
##random Forest Procedure

###########RF procedure to test classes

##using selected genes diff express

to.train<-to.train[rownames(to.train) %in% proms.island,]##proms.islands
to.train<-to.train[,match(training.set$Name,colnames(to.train))]

##LIMMA

###implement limma best
y <-as.factor(training.set[,"RESP"])

design.m<-model.matrix(~y)
fit<-lmFit(to.train,design.m)
fit<-eBayes(fit)
res<-topTable(fit,n= dim(fit)[1],adjust.method="BH", sort.by="P")  
res<-res[res$P.Value<.01,]

pv.sig<-rownames(res)



pca <- prcomp(t(to.train[pv.sig,]))
pca_data_perc<-round(100*((pca$sdev)^2/sum(pca$sdev^2)),digits=2)



df_pca_data<-data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], batch=training.set$RESP)

ggplot(df_pca_data, aes(PC1,PC2, color = batch))+
  geom_point(size=3)+
  labs(x=paste0("PC1 (",pca_data_perc[1],")"), y=paste0("PC2 (",pca_data_perc[2],")"))


####SVM-RF filtering
y <-as.factor(training.set[,"RESP"])
library(sigFeature)

sigfeatureRankedList <- sigFeature(t(to.train), y)
save(sigfeatureRankedList,file="SVM_RF_filtering.RData")

sels<-rownames(to.train)sigfeatureRankedList[1:100]
##03 Build RF model

ntrees<-999
##define response
y <-as.factor(training.set[,"RESP"])


rf.pred <- randomForest(
 # t(to.train[pv.sig,]),
  t(to.train),
  y,
  #mc=8,
  ntree=ntrees,
#  strata=y,
  mtry=sqrt(ncol(t(to.train))),
 # mtry=100,
  sampsize=rep(min(table(y)),length(table(y))),#Balancing by sampling stratification
#cutoff=c(1-rare.class.prevalence,rare.class.prevalence),##Balancing by voting rule
  proximity=TRUE,
  oob.prox=TRUE,
  importance=TRUE,
  keep.inbag=TRUE,
  do.trace=FALSE,
  seed=seed
)


##test on training


################make Folds

folds<-3
nfolds <- makenestedfolds(y,folds)
cores<-8


# check minimal class sizes for inner training loops
minclasssize <- matrix(0,ncol=length(nfolds),nrow=length(nfolds))

for(i in 1:length(nfolds)){
  for(j in 1:length(nfolds))
    minclasssize[i,j]  <- min(table(y[nfolds[[i]][[2]][[j]]$train]))
}
colnames(minclasssize) <- paste0("innfold",1:folds)
rownames(minclasssize) <- paste0("fold",1:folds)

print(minclasssize)





####select genes in each external fold, select the one with less error rate and construct RF


setwd("../Predictor_selinside_fold_antiTNF/")
ntrees <- 999
p <- 50 ##select 1st 1000 with High Gini Scores
y <-as.factor(training.set[,"RESP"])

#-----------------------------------------------------------------------------------
# script to train the classifier in each CV fold, including feature selection.
#
# Note, in this example we reduced the number of features to 20k probes by sd filtering before applying the random forest for 
# feature selection. To perform feature selection as described in the paper remove line 25,
# this will increase the computation time significantly.
# 
# Martin Sill
# m.sill@dkfz.de                                                                  
# 
# 2018-03-14 UTC
#------------------------------------------------------------------------------------                   
#source("../../Classes_RNAseq_2019/Predictor_scripts/feature_selection_in_each_fold.R")


##now run CV


for(K in 1:folds){
  
  for(k in 0:folds){
    
    if(k>0){  message("calculating fold ",K,".",k,"  ...",Sys.time())
      fold <- nfolds[[K]][[2]][[k]]
    }else{
      message("calculating outer fold ",K,"  ...",Sys.time())
      fold <- nfolds[[K]][[1]][[1]]
    }
    
    rf.scores <- calcultateCVfold(t(to.train),y,fold,p,cores,ntrees)###original

    save(rf.scores,file=paste("CVfold",K,k,"RData",sep="."))
    
    gc()
  }
}

#this script fits for each CV fold a multinomial L2-penalized glmnet model to calibrate RF scores
# and one final calibration model using RF scores generated in the out loop of the CV 
# After calibration CV results are calculated in CVresults.Rmd
#


######calibration


cores <- 8

registerDoParallel(cores)

#message("loading data ...",Sys.time())
anno<-y

for(i in 1:length(nfolds)){
  scores <- list() 
  idx <- list()
  for(j in 1:length(nfolds)){
    load(paste0("CVfold.",i,".",j,".RData"))
    scores[[j]] <- rf.scores
    idx[[j]] <- nfolds[[i]][[2]][[j]]$test
  }
  scores <- do.call(rbind,scores)
  idx <- unlist(idx)
  y <- anno[idx]         
  
  message("fitting calibration model fold ",i," ...",Sys.time())
  # fit multinomial logistic ridge regression model
  suppressWarnings(cv.calfit <- cv.glmnet(y=y,x=scores,family="multinomial",type.measure="mse",
                                          alpha=0,nlambda=100,lambda.min.ratio=10^-6,parallel=TRUE))
  
  
  load(paste0("CVfold.",i,".",0,".RData"))
  
  message("calibrating raw scores fold ",i," ...",Sys.time())
  probs <- predict(cv.calfit$glmnet.fit,newx=rf.scores,type="response"
                   ,s=cv.calfit$lambda.1se)[,,1] # use lambda estimated by 10fold CVlambda
  
  
  err <- sum(colnames(probs)[apply(probs,1,which.max)] != anno[nfolds[[i]][[1]][[1]]$test])/length(nfolds[[i]][[1]][[1]]$test)
  
  message("misclassification error: ",err)
  
  save(probs,file=paste0("probsCVfoldB.",i,".",0,".RData"))
}

scores <- list()
idx <- list()

for(i in 1:length(nfolds)){
  load(paste0("CVfold.",i,".",0,".RData"))
  scores[[i]] <- rf.scores
  idx[[i]] <- nfolds[[i]][[1]][[1]]$test
}
scores <- do.call(rbind,scores)

probl <- list()
for(i in 1:length(nfolds)){
  load(paste0("probsCVfoldB.",i,".",0,".RData"))
  probl[[i]] <- probs
}
probs <- do.call(rbind,probl)


idx <- unlist(idx)
y <- anno[idx] 

ys <- colnames(scores)[apply(scores,1,which.max)]
yp <- colnames(probs)[apply(probs,1,which.max)]


errs <- sum(y!=ys)/length(y)
errp <- sum(y!=yp)/length(y)

message("overall misclassification error scores: ",errs)
message("overall misclassification error calibrated: ",errp)

suppressWarnings(cv.calfit <- cv.glmnet(y=y,x=scores,family="multinomial",type.measure="mse",
                                        alpha=0,nlambda=100,lambda.min.ratio=10^-6,parallel=TRUE))

save(cv.calfit,file="calfit_FIS_adj.RData")

save(rf.pred,scores,probs,y,ys,yp,file="CVresults_Pred_FIS_adj.RData")


###try to plot results



AUCscores <- HandTill2001::auc(multcap(response=as.factor(y),predicted=scores))
AUCprobs <-  HandTill2001::auc(multcap(response=as.factor(y),predicted=probs))

errs <- sum(y!=ys)/length(y)
errp <- sum(y!=yp)/length(y)

briers <- brier(scores,y)
brierp <- brier(probs,y) 

logls <- mlogloss(scores,y)
loglp <- mlogloss(probs,y) 

out <- cbind(c(AUCscores,AUCprobs),c(errs,errp),c(briers,brierp),c(logls,loglp))
colnames(out) <- c("auc","misclassification","brier score","logloss")
rownames(out) <- c("raw scores","calibrated scores")
out
write.csv(out, "results_Pred_FIS_adj_w0.csv")


y.true <- y
y.score <- ys
y.calibrated <- yp
table(y.true,y.score)
table(y.true,y.calibrated)

par(mfrow=c(1,2))
for(i in 1:2){
  plot(x=scores[,i],y=probs[,i],ylim=c(0,1),pch=16,xlim=c(0,1),col=as.factor((y==colnames(scores)[i])),main=colnames(scores)[i] ,xlab="scores",ylab="calibrated scores")
  abline(a=0,b=1,col="green",lty=2)
}

predictions<-as.data.frame(scores)

predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$predict <-ys
predictions$observed <- y
head(predictions)

dev.off()
save(predictions,file="Predictions_FIS_Diff_adj.RData")
#  ROC curve,
library(plotROC)
rocplot <- ggplot(predictions, aes(m = R, d = observed))+ geom_roc(n.cuts=20,labels=FALSE)
rocplot + style_roc(theme = theme_grey) + geom_rocci(fill="pink") +
  annotate("text", x = .75, y = .25, 
           label = paste("AUC =", round(calc_auc(rocplot)$AUC, 2))) +
  scale_x_continuous("1 - Specificity", breaks = seq(0, 1, by = .1))


rocplot

load("/mnt/md127/Analysis_Projects/FIS_Meth/Classes_Results/FIS_P7/EPIC_annotation.RData")
##DMP Annotation


##prune those ones with no importance

prune.var<-as.data.frame(rf.pred$importance)

prune.var<-prune.var[!(prune.var$MeanDecreaseGini==0),]
prune.var<-prune.var[prune.var$MeanDecreaseGini > mean(prune.var$MeanDecreaseGini),]

pv.sig<-rownames(prune.var)


dmp <- merge(pv.sig,ann850k[,c("UCSC_RefGene_Name","Relation_to_Island","UCSC_RefGene_Group")], by="row.names")
rownames(dmp)<-dmp$Row.names
dmp<-dmp[,-1]
write.csv(dmp,"DMPS_Predictor_Adj.csv")
#####GO
library(methylGSA)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)


pv.vals<- 0.005

cpgval<-pv.vals[rownames(rf.pred$importance)]

res1 <- methylglm(cpg.pval = cpgval, minsize = 10, maxsize = 1000, GS.type = "GO",array.type = "EPIC")

res1<-res1[res1$padj<.05,]

write.csv(res1,"GO_Predictor_Selected_cpgs.csv")

################################LOAD prediction P6 Data


##calibrate with glmnet and give response

#probs <- predict(cv.calfit$glmnet.fit,newx=RF_predictions_probs,type="response",s=cv.calfit$lambda.1se)[,,1] # use lambda estimated by 5fold CVlambda

##05 predict

##Select P6 samples to predict

predictor_data<-t(betaAdj[rownames(betaAdj) %in% rownames(rf.pred$importance),colnames(betaAdj) %in% test.set$Name])##
#predictor_data<-t(betaAdj[rownames(betaAdj) %in% rownames(rf.pred$importance),colnames(betaAdj) %in% training.set$Name])##

RF_predictions_probs<-predict(rf.pred, predictor_data, type="p")

probs <- predict(cv.calfit$glmnet.fit,newx=RF_predictions_probs,type="response",s=cv.calfit$lambda.1se)[,,1] # use lambda estimated by 5fold CVlambda

RF_predictions_responses<-predict(rf.pred, predictor_data, type="response")

RF_p<-as.data.frame(probs)
RF_p$Name<-rownames(RF_p)
RF_a<-as.data.frame(RF_predictions_responses)
RF_a$Name<-rownames(RF_a)

sols<-merge(RF_a,targets,by="Name")
sols2<-merge(RF_a,targets,by="Name")



predictions<-sols[,c("RF_predictions_responses","RESP","Name")]
predictions<-merge(predictions,RF_p,by="Name")

##select only those ones above 0.75


####calibrate predictions
predictions$RESP<-as.factor(predictions$RESP)

stats.preds<-confusionMatrix(predictions$RF_predictions_responses,  
                             as.factor(predictions$RESP))

#####print results
write.csv(RF_p,"Probs_prediction_Test.csv")  



######TEST

rocplot <- ggplot(predictions, aes(m = R, d = RESP))+ geom_roc(n.cuts=20,labels=FALSE)
rocplot + style_roc(theme = theme_grey) + geom_rocci(fill="pink") +
  annotate("text", x = .75, y = .25, 
           label = paste("AUC =", round(calc_auc(rocplot)$AUC, 2))) +
  scale_x_continuous("1 - Specificity", breaks = seq(0, 1, by = .1))


rocplot
