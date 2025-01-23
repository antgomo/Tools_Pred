setwd("/Users/agomez/IMIDOMICS/FIS/")
library(limma)
library(minfi)
library(DMRcate)
##function to calculate CV in folds
calcultateCVfold <- function(badj,y,fold,p,cores,ntrees){
  
  
  # sd pre filtering to 20k probes, to speed up the example
  #  badj$betas.train <- badj$betas.train[,order(-apply(badj$betas.train,2,sd))[1:20000]]
  
  message("performing variable selection ...",Sys.time())
  message("cores: ",cores)
  message("ntrees: ",ntrees)  
  message("n: ",nrow(badj))
  message("p: ",ncol(badj))  
  
  rf.varsel <- rfp(badj[fold$train,],
                   y=y[fold$train],
                   mc=cores,
                   ntree=ntrees,
                   sampsize=rep(min(table(y[fold$train])),length(table(y[fold$train]))),
                   importance=TRUE)
  
  # get permutation variable importance
  imp.meandecrease <- rf.varsel$importance[,dim(rf.varsel$importance)[2]-1]
  
  # reduce data matrix
  or <- order(imp.meandecrease,decreasing=T)
  
  message("training classifier ...",Sys.time())
  message("cores: ",cores)
  message("ntrees: ",ntrees)  
  message("n: ",nrow(badj))
  message("p: ",p)  
  
  rf <- rfp(badj[fold$train,or[1:p]],y[fold$train],
            sampsize=rep(min(table(y[fold$train])),length(table(y[fold$train])))
            ,mc=cores,ntree=ntrees,importance=TRUE) 
  
  message("predicting test set ...",Sys.time())
  
  rf.scores <- predict(rf,badj[fold$test,match(rownames(rf$importance),
                                               colnames(badj[fold$test,]))],type="prob")
  
  err <- sum(colnames(rf.scores)[apply(rf.scores,1,which.max)]!=y[fold$test])/length(fold$test)
  message("misclassification error: ",err)
  
  return(rf.scores)
}


######################################end of function definition

####load pheno
load("phenoData.RData")


idat.folder <- "/Users/agomez/IMIDOMICS/FIS/" ###change for wour path dir

targets <- read.metharray.sheet(base=idat.folder)
targets$Basename<-paste(idat.folder,targets$xipName, sep="")

load("FIS_norm_cellcomp_data.RData")
###begin here after normalized beta pipeline is done
####load pheno
#load("eular.RData") # load eular 
#eular.df<-as.data.frame(eular)
####get only basal week

##only betas in blood target

targets<-targets[targets$tissue=="blood",]
betas <-betas[,colnames(betas) %in% targets$xipName]
targets <-targets[targets$xipName %in% colnames(betas),]

betas<-betas[,match(targets$xipName,colnames(betas))]

pD.w0<-targets[grep("S0",targets$donation,perl=T),]

##remove sexual chromosomes and SNPs probes

betas<-rmSNPandCH(betas, dist = 2, mafcut = 0.05, and = TRUE, rmcrosshyb = TRUE, rmXY=TRUE)



##select betas depending week
betas.fis<-betas[,colnames(betas) %in% pD.w0$xipName]

##put in same order as phenodata to get clear Sample names
pD.w0$Code<-gsub("-S0","",pD.w0$donation)
betas.fis<-betas.fis[,match(as.character(pD.w0$xipName),colnames(betas.fis))]#w0
colnames(betas.fis)<-pD.w0$Code


phenoData<-read.csv("phenodata_FIS_fixed.csv")
phenoData$delta<-phenoData$S0_DAS28-phenoData$SR12_DAS28
phenoData$EULAR2<-ifelse(phenoData$SR12_DAS28<3.2 & phenoData$delta>1.2,"GOOD",
                         ifelse(phenoData$SR12_DAS28<3.2 & phenoData$delta>0.6 & phenoData$delta<1.2 ,"MOD",
                                ifelse(phenoData$SR12_DAS28<3.2 & phenoData$delta<0.6,"NON_RESP",
                                       ifelse(phenoData$SR12_DAS28>3.2 & phenoData$SR12_DAS28<5.1 & phenoData$delta>1.2,"MOD",
                                              ifelse(phenoData$SR12_DAS28>3.2 &  phenoData$SR12_DAS28<5.1 & phenoData$delta>0.6 & phenoData$delta<1.2,"MOD",
                                                     ifelse(phenoData$SR12_DAS28>3.2 & phenoData$SR12_DAS28<5.1 & phenoData$delta<0.6,"NON_RESP",
                                                            ifelse(phenoData$SR12_DAS28>5.1 & phenoData$delta>1.2,"MOD",
                                                                   ifelse(phenoData$SR12_DAS28> 5.1 & phenoData$delta>0.6 & phenoData$delta<1.2 ,"NON_RESP",
                                                                          ifelse(phenoData$SR12_DAS28> 5.1 & phenoData$delta<0.6,"NON_RESP","NA")))))))))

phenoData$RESP<-ifelse(phenoData$EULAR2=="NON_RESP","NR",ifelse(phenoData$EULAR2=="NA","NA","R"))


library(limma)

##simple limma approach, build design matrix with Cell composition

phenoData3<-phenoData[!(is.na(phenoData$Age)),]
phenoData3<-phenoData3[!(is.na(phenoData3$RESP)),]

##in case to test w12

phenoData3<-phenoData3[phenoData3$Code %in% colnames(betas.fis),]

betas.fis<-betas.fis[,colnames(betas.fis) %in% phenoData3$Code]

betas.fis<-betas.fis[,match(phenoData3$Code,colnames(betas.fis))]

###adjust betas
design <- model.matrix(~phenoData3$RESP+ phenoData3$SEX + phenoData3$Age + cell_comp[rownames(cell_comp) %in% phenoData3$xipName,])##only sex

fit <- lmFit(betas.fis, design)
mAdj.fit    <- fit$coefficients[,-c(1,2)]

mAdj<- as.matrix(betas.fis - mAdj.fit %*% t(design[,-c(1,2)]))

betaAdj     <- ilogit2(mAdj)

rm(fit)
rm(mAdj)
rm(mAdj.fit)
rm(betas)
gc()



###select to train
#-----------------------------------------------------------------------------------
# parallelized version of random Forest
#                                                                     
# Martin Sill
# m.sill@dkfz.de                                                                  
# 
# 2018-03-14 UTC
#------------------------------------------------------------------------------------                   
#options(max.print = 1000)
#options(stringsAsFactors = FALSE)
#options(scipen = 999)

library(randomForest)
library(parallel)

rfp <- function(xx, ..., ntree = ntree, mc = mc, seed = NULL)
{
  if(!is.null(seed)) set.seed(seed, "L'Ecuyer")
  rfwrap <- function(ntree, xx, ...) randomForest::randomForest(x=xx,ntree=ntree,norm.votes=FALSE, ...)
  rfpar <- mclapply(rep(ceiling(ntree/mc),mc),mc.cores=mc, rfwrap, xx=xx, ...)
  do.call(randomForest::combine, rfpar)
}

ntrees <- 10000  # 10000 in the paper, here 500 to speed up the example
cores <- 8
##first select variables
y <- as.factor(phenoData3$RESP)

##add cell comp

cell_comp<-cell_comp[rownames(cell_comp) %in% phenoData3$xipName,]

cell_comp<-cell_comp[match(phenoData3$xipName,rownames(cell_comp)),]
rownames(cell_comp)<-phenoData3$Code
##transpose to add to betAdj
cell_comp<-t(as.data.frame(cell_comp))


betaAdj<-rbind(betaAdj,cell_comp)

library(varSelRF)

#rfsel<- varSelRF(t(betaAdj),y,ntree=10000, ntreeIterat=2000, vars.drop.frac=0.2,whole.range = FALSE,keep.forest = TRUE) 

#plot(rfsel)
rfsel<-read.csv("Selected_vars_Predictor_s0.csv",head=T)
#rfsel<-read.csv("R_NR_w0_sex_age.csv",head=T)

rfsel<-rfsel$x###8 sel
#rfsel<-rfsel$X
rfsel<-c(as.character(rfsel),"CD8T","CD4T","NK","Bcell","Mono","Neu")


##save(rfsel, file="Vars_sel_FIS.RData")
#write.csv(rfsel$selected.vars, "Selected_vars_Predictor_s0.csv")
to.train<-t(betaAdj[rownames(betaAdj) %in% rfsel,])



print(colnames(to.train))


rf.pred <- randomForest(
  to.train,
  y,
  #mc=8,
  ntree=ntrees,
  #strata=y,
  mtry=sqrt(ncol(to.train)),
  #mtry=100,
  sampsize=rep(min(table(y)),length(table(y))),#Balancing by sampling stratification
  proximity=TRUE,
  oob.prox=TRUE,
  importance=TRUE,
  keep.inbag=TRUE,
  do.trace=FALSE,
  seed=seed
)

save(rf.pred,file="Predictor_8cpgs.RData")

# get permutation variable importance
imp.meandecrease <- rf.pred$importance[,dim(rf.pred$importance)[2]-1]

# save selection forest

#######################From here, same pocedure as Nature, 2018

library(limma)

ntrees <- 10000
cores <- 3
seed <- 180314
p <- dim(to.train)[2]
folds <- 5
y <- as.factor(phenoData3$RESP)

source("makefolds.R")
#source("calculateCVfold.R")


nfolds <- makenestedfolds(y,folds)



# check minimal class sizes for inner training loops
minclasssize <- matrix(0,ncol=length(nfolds),nrow=length(nfolds))

for(i in 1:length(nfolds)){
  for(j in 1:length(nfolds))
    minclasssize[i,j]  <- min(table(y[nfolds[[i]][[2]][[j]]$train]))
}
colnames(minclasssize) <- paste0("innfold",1:folds)
rownames(minclasssize) <- paste0("fold",1:folds)
print(minclasssize)


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


##now run CV


for(K in 1:folds){
  
  for(k in 0:folds){
    
    if(k>0){  message("calculating fold ",K,".",k,"  ...",Sys.time())
      fold <- nfolds[[K]][[2]][[k]]
    }else{
      message("calculating outer fold ",K,"  ...",Sys.time())
      fold <- nfolds[[K]][[1]][[1]]
    }
    
    rf.scores <- calcultateCVfold(to.train,y,fold,p,cores,ntrees)
    
    save(rf.scores,file=paste("CVfold",K,k,"RData",sep="."))
    
    gc()
  }
}

#this script fits for each CV fold a multinomial L2-penalized glmnet model to calibrate RF scores
# and one final calibration model using RF scores generated in the out loop of the CV 
# After calibration CV results are calculated in CVresults.Rmd which is compiled to an html report
#


######calibration
library(rmarkdown)
library(glmnet)
library(doParallel)
library(HandTill2001)

cores <- 3

registerDoParallel(cores)

#message("loading data ...",Sys.time())
#load("./results/Mset_filtered.RData")
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
  
  save(probs,file=paste0("probsCVfold.",i,".",0,".RData"))
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
  load(paste0("probsCVfold.",i,".",0,".RData"))
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

save(cv.calfit,file="calfit_FIS_Diff.RData")

save(scores,probs,y,ys,yp,file="CVresults_Pred_FIS_Diff.RData")


###try to plot results

#rm(list=ls())
library(pROC)
library(HandTill2001)
source("multi_brier.R")
source("multi_logloss.R")

#load("CVresults_Pred_FIS.RData")

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
write.csv(out, "results_Pred_FIS_8_cpgs_cellcomp_w0.csv")


######TEST

# predict class and then attach test class
predictions<-as.data.frame(scores)
predictions<-as.data.frame(probs)

predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$predict <-ys
predictions$observed <- y
head(predictions)

save(predictions,file="Predictions_FIS_Diff.RData")
# 1 ROC curve, mock vs non mock
roc.mock <- roc(ifelse(predictions$observed=="R", "R", "NR"), as.numeric(predictions$R))
plot(roc.mock, col = "red",main="R vs NR Week 0 Cell Comp")
text(0.3,0.9,paste("AUC = ",format(AUCprobs, digits=3, scientific=FALSE)))
 

#####GO enrichment analysis
library(ggplot2)

#####GO enrichment analysis

library(missMethyl)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
gometh.res<-gometh(colnames(to.train), all.cpg = NULL, collection ="GO", array.type ="EPIC", plot.bias = T, prior.prob = T)


gometh.res<-gometh.res[gometh.res$Ont=="BP",]
gometh.res<-gometh.res[gometh.res$P.DE<.05,]
gometh.res$logFDR<-(-log(gometh.res$P.DE))
gometh.res$GeneRatio<-gometh.res$DE/gometh.res$N



gometh.res<-gometh.res[!(grepl("development",gometh.res$Term)),]
gometh.res<-gometh.res[!(grepl("transcription",gometh.res$Term)),]
gometh.res<-gometh.res[!(grepl("kinase",gometh.res$Term)),]
gometh.res<-gometh.res[!(grepl("transport",gometh.res$Term)),]
gometh.res<-gometh.res[!(grepl("axon",gometh.res$Term)),]

gometh.res<-gometh.res[order(gometh.res$GeneRatio),]


gometh.res$Term <- factor(gometh.res$Term, levels=rev(gometh.res$Term))###ordering from top to last



###Forest_plot GOmeth
fp<-ggplot(data=gometh.res[1:55,],aes(x=GeneRatio,y=Term))+
  geom_point(aes(size=DE,fill=P.DE), colour="black",shape=21)+
  theme(axis.text.y = element_text(margin=margin(.1,.5,.8,.5,"pt"),size=7))+
  ggtitle("GO DAS28 Anytime point")+
  ylab(NULL) +
  theme(plot.margin=unit(c(.01,1,.1,1.2),"cm"))

