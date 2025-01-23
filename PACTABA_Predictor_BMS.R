library(tximport)
library(sva)
library(RUVSeq)
library(splines)
library(limma)


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

setwd("/mnt/md127/Analysis_Projects/Pactaba/")

sampleTable<-read.csv("PACTABA_RNAseq_Mastersheet.csv")
rownames(sampleTable)<-sampleTable$GSL.ID

codes<-read.csv("Correspondencia_codis.csv")
sampleTable<-merge(sampleTable,codes,by="rna.code")
###pheno

pheno<-read.csv("Pheno_sex_age.csv")
colnames(pheno)[3]<-"Code"
##2 is female
sampleTable<-merge(sampleTable,pheno,by="Code",all.x=T)
rownames(sampleTable)<-sampleTable$GSL.ID

week0<-read.csv("Week0_BMS.csv")
week12<-read.csv("Week12_BMS.csv")

####add 0 and 12 week
sampleTable<-merge(sampleTable,week0,by="Code",all.x=T)
sampleTable<-merge(sampleTable,week12,by="Code",all.x=T)
sampleTable<-sampleTable[!(duplicated(sampleTable$GSL.ID)),]

rownames(sampleTable)<-sampleTable$GSL.ID


##########establish EULAR

sampleTable$delta<-sampleTable$das28esr_w0-sampleTable$das28esr_w12
sampleTable$EULAR<-ifelse(sampleTable$das28esr_w12<=3.2 & sampleTable$delta>1.2,"GOOD",
                          ifelse(sampleTable$das28esr_w12<=3.2 & sampleTable$delta>0.6 & sampleTable$delta<=1.2 ,"MOD",
                                 ifelse(sampleTable$das28esr_w12<=3.2 & sampleTable$delta<=0.6,"NON_RESP",
                                        ifelse(sampleTable$das28esr_w12>3.2 & sampleTable$das28esr_w12<=5.1 & sampleTable$delta>1.2,"MOD",
                                               ifelse(sampleTable$das28esr_w12>3.2 &  sampleTable$das28esr_w12<=5.1 & sampleTable$delta>0.6 & sampleTable$delta<=1.2,"MOD",
                                                      ifelse(sampleTable$das28esr_w12>3.2 & sampleTable$das28esr_w12<=5.1 & sampleTable$delta<=0.6,"NON_RESP",
                                                             ifelse(sampleTable$das28esr_w12>5.1 & sampleTable$delta>1.2,"MOD",
                                                                    ifelse(sampleTable$das28esr_w12> 5.1 & sampleTable$delta>0.6 & sampleTable$delta<=1.2 ,"NON_RESP",
                                                                           ifelse(sampleTable$das28esr_w12> 5.1 & sampleTable$delta<=0.6,"NON_RESP","NA")))))))))
##write out oheno PACTABA to check
pheno_PACTABA<-sampleTable[,c("Code","rna.code","donor.code.x","donationId","das28esr_w0","das28esr_w12","delta","EULAR")]
pheno_PACTABA<-pheno_PACTABA[!(duplicated(pheno_PACTABA)),]

library(edgeR)


all_cts <- read.delim("Counts_PACTABA.txt",head=T)
colnames(all_cts)<-gsub("X","",colnames(all_cts))
colnames(all_cts)<-gsub(".bam","",colnames(all_cts))
colnames(all_cts)<-gsub("\\.","-",colnames(all_cts))
rownames(all_cts)<-all_cts$Geneid
all_cts<-all_cts[,-1]
all_cts<-all_cts[,-1]

colnames(all_cts)<-substr(colnames(all_cts),1,13)

sampleTable<-sampleTable[rownames(sampleTable) %in% colnames(all_cts),]
sampleTable$Sex<-ifelse(sampleTable$sex==2,"F",ifelse(sampleTable$sex==1,"M","NA"))

all_cts<-all_cts[,colnames(all_cts) %in% rownames(sampleTable)]
#colnames(all_cts)<-sampleTable$rna.code


###CONDITION 2 : R= R+MOD and NR
sampleTable$COND<-ifelse(sampleTable$EULAR=="GOOD"|sampleTable$EULAR=="MOD","G","NR")

###get only in S0

#samples2test<-sampleTable[grep("S0|S12",sampleTable$donationId,perl = T),]
samples2test<-sampleTable[grep("S0",sampleTable$donationId,perl = T),]

samples2test<-samples2test[!(is.na(samples2test$COND)),]


####Analysis

cts<-all_cts[,colnames(all_cts) %in% rownames(samples2test)]
cts<-cts[,match(rownames(samples2test),colnames(cts))]
y <- DGEList(cts,group=samples2test$COND)## for COND2

# filter out very lowly expressed tags, keeping genes that are expressed at a reasonable level in at least
#one treatment condition. Since the smallest group size is 30, we keep genes that achieve at least
#one count per million (cpm) in at least 30 samples:

keep <- rowSums(cpm(y)>2) >= 15#### change group depending replicates
y <- y[keep, ]

dim(y)

y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y, method="TMM")


logCPM <- cpm(y, log=TRUE, prior.count=.1) 


####Limma
my.batch<-as.factor(paste("B", samples2test$seq.plate,sep=""))

##Batch removal
samples2test$Code<-gsub("-","",samples2test$Code)


design <- model.matrix(~samples2test$COND)##simple


# Remove batch effect. 
logCPM_no_batch <- removeBatchEffect(logCPM, batch=my.batch, design=design) 


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
y <- as.factor(samples2test$COND)

set.seed(seed,kind ="L'Ecuyer-CMRG")
rf.varsel <- rfp(t(logCPM_no_batch),
                 y,
                 mc=cores,
                 ntree=ntrees,
                 sampsize=rep(min(table(y)),length(table(y))),
                 importance=TRUE)

# get permutation variable importance
imp.meandecrease <- rf.varsel$importance[,dim(rf.varsel$importance)[2]-1]


# reduce data matrix to the first 10000 Cpgs depending on its importance
or <- order(imp.meandecrease,decreasing=T)

to.train<-t(logCPM_no_batch[or[1:5000],])
library(varSelRF)
rfsel<- varSelRF(t(logCPM_no_batch),y,ntree=10000, ntreeIterat=2000, vars.drop.frac=0.2,whole.range = FALSE,keep.forest = TRUE) 

plot(rfsel)

write.csv(rfsel$selected.vars, "Results_BMS/Selected_vars_Predictor_s0.csv")
to.train<-t(logCPM_no_batch[rownames(logCPM_no_batch) %in% rfsel$selected.vars,])

#######################From here, same pocedure as Nature, 2018

library(limma)

ntrees <- 10000
cores <- 8
seed <- 180314
p <- dim(to.train)[1]
folds <- 3
y <- as.factor(samples2test$COND)

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

cores <- 8

registerDoParallel(cores)

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

save(cv.calfit,file="calfit_P6.RData")

save(scores,probs,y,ys,yp,file="CVresults_Pred_P6.RData")


###try to plot results

rm(list=ls())
library(pROC)
library(HandTill2001)
source("multi_brier.R")
source("multi_logloss.R")

load("CVresults_Pred_P6.RData")

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
write.csv(out, "Results_BMS/results_Pred_PACTABA_w0_12.csv")


######TEST

# predict class and then attach test class
predictions<-as.data.frame(scores)
predictions<-as.data.frame(probs)

predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$predict <-ys
predictions$observed <- y
head(predictions)
# 1 ROC curve, mock vs non mock
roc.mock <- roc(ifelse(predictions$observed=="G", "G", "NR"), as.numeric(predictions$G))
plot(roc.mock, col = "red",main="PACTABA Gvs NR Week 0 EULAR w0_w12")
text(0.3,0.8,paste("AUC = ",format(roc.mock$auc, digits=3, scientific=FALSE)))
 

#####GO enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)
#library(GOstats)
library(GO.db)
library(DOSE)
library(annotate)

genesid<-colnames(to.train)

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


dotplot(ego2, title="PACTABA Predictor selected vars",showCategory=50)
