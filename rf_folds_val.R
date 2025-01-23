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
setwd("/mnt/md127/Analysis_Projects/Classes_RNAseq_2019/")

library(randomForest)
library(parallel)
library(ggplot2)

rfp <- function(xx, ..., ntree = ntree, mc = mc, seed = NULL)
{
  if(!is.null(seed)) set.seed(seed, "L'Ecuyer")
  rfwrap <- function(ntree, xx, ...) randomForest::randomForest(x=xx,ntree=ntree,norm.votes=FALSE, ...)
  rfpar <- mclapply(rep(ceiling(ntree/mc),mc),mc.cores=mc, rfwrap, xx=xx, ...)
  do.call(randomForest::combine, rfpar)
}

load("SLE_HI_adj.RData")
load("SLE_LO_adj.RData")


annon<-data.pred[[2]]

classes <- paste("C",annon$consensuscluster,sep="_")
classes<-as.factor(classes)

data.tp<-hi.adj
###load data to test
data.tp<-data.tp[,match(rownames(annon),colnames(data.tp))]

ntrees <- 10000  # 10000 in the paper, here 500 to speed up the example
cores <- 8


##first step, define the differentially expressed genes among all classes
##AOV to identify genes that are variable among and between groups.

pv <- sapply(1:dim(data.tp)[1], function(i) {
  mydataframe <- data.frame(y=data.tp[i,], ig=classes)
  fit <- aov(y ~ ig, data=mydataframe)
  summary(fit)[[1]][["Pr(>F)"]][1]
})
names(pv) <- rownames(data.tp)
pv.sig <- names(pv)[pv < 0.05/1:dim(data.tp)[1]/length(classes)] ## bonferonni

##using selected genes diff express


to.train<-data.tp[rownames(data.tp) %in% pv.sig,]

#heatmap(data.tp[pv.sig,], Rowv=NA, Colv=NA, col=colorRampPalette(c('blue', 'white', 'red'))(100), scale="none", ColSideColors=rainbow(length(classes))[classes], labCol=FALSE, labRow=FALSE)# save selection forest

#######################From here, same pocedure as Nature, 2018, 5 fold internal cross validation

library(limma)

ntrees <- 999
cores <- 8
seed <- 180314
p <-250 ##number of maximum genes to select in each fold
folds <- 5
y <- classes
source("Predictor_scripts/makefolds.R")

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
# Martin Sill
# m.sill@dkfz.de                                                                  
# 
# 2018-03-14 UTC
#------------------------------------------------------------------------------------                   


##now run CV

####select genes in each external fold, select the one with less error rate and construct RF

genes<-list()
errors<-c()
rf.list<-list()


folds.genes <- makefolds(y,5)


for(K in 1:folds){
  
  message("calculating outer fold for selecting genes",K,"  ...",Sys.time())

  fold <- folds.genes[[K]]
  
  rf.scores <- calcultateCVfold.genes(t(to.train),as.factor(classes),fold,p,cores,ntrees)
  genes[[K]]<-rf.scores[[2]]
  errors[K]<-rf.scores[[3]]
  rf.list[[K]]<-rf.scores[[4]]

}

##nested fold to calculate predictions and validate


for(K in 1:folds){
  
  for(k in 0:folds){
    
   if(k>0){  message("calculating fold ",K,".",k,"  ...",Sys.time())
     fold <- nfolds[[K]][[2]][[k]]
    }else{
      message("calculating outer fold ",K,"  ...",Sys.time())
      fold <- nfolds[[K]][[1]][[1]]
    }

    rf.scores <- calcultateCVfold(t(to.train),as.factor(classes),fold,p,cores,ntrees)
    save(rf.scores,file=paste("CVfold",K,k,"RData",sep="."))
    
   gc()
  }
}


#this script fits for each CV fold a multinomial L2-penalized glmnet model to calibrate RF scores
# and one final calibration model using RF scores generated in the out loop of the CV 
# After calibration CV results are calculated in CVresults.Rmd which is compiled to an html report
#



#In our paper we fitted an additional calibration model that maps the raw random forest scores to more suitable probabilities and this was validated by a nested cross validation.
#To get ‘indepedent’ random forest scores that are needed to fit such a calibration model we used the RF scores generated in the outer CV loop and 
#the inner loops are used to fit calibration models to validate the calibration step.
#The cv.calfit object in the calibration.R script in line 91 is the model fit on the outer CV loop RF scores that is later used to predict scores from the final classifier 
#that was trained on the complete data set.
#In line 42 the inner loop calibration models are fitted, that are used for validation.


######calibration
library(rmarkdown)
library(glmnet)
library(doParallel)
library(HandTill2001)

cores <- 8

registerDoParallel(cores)

anno<-classes

for(i in 1:length(nfolds)){
  scores <- list() 
  idx <- list()
  for(j in 1:length(nfolds)){
    load(paste0("CVfold.",i,".",j,".RData"))
    
    scores[[j]] <- rf.scores[[1]]
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
  probs <- predict(cv.calfit$glmnet.fit,newx=rf.scores[[1]],type="response"
                   ,s=cv.calfit$lambda.1se)[,,1] # use lambda estimated by 10fold CVlambda
  
  
  err <- sum(colnames(probs)[apply(probs,1,which.max)] != anno[nfolds[[i]][[1]][[1]]$test])/length(nfolds[[i]][[1]][[1]]$test)
  
  message("misclassification error: ",err)
  
  save(probs,file=paste0("probsCVfold.",i,".",0,".RData"))
}

###Final calibration on FR build with best genes


scores <- list()
idx <- list()

for(i in 1:length(nfolds)){
  load(paste0("CVfold.",i,".",0,".RData"))
  scores[[i]] <- rf.scores[[1]]
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

save(cv.calfit,file="calfit_FINAL.RData")


save(scores,probs,y,ys,yp,file="CVresults_Pred_FIS_Boruta.RData")##calibrate final predictions


###try to plot results

#rm(list=ls())
library(pROC)
library(HandTill2001)
source("Predictor_scripts/multi_brier.R")
source("Predictor_scripts/multi_logloss.R")

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
write.csv(out, "Predictor_balanced/results_Pred_FIS_Balanced.csv")

y.true <- y
y.score <- ys
y.calibrated <- yp
table(y.true,y.score)
table(y.true,y.calibrated)

par(mfrow=c(2,3))
for(i in 1:5){
  plot(x=scores[,i],y=probs[,i],ylim=c(0,1),pch=16,xlim=c(0,1),col=as.factor((y==colnames(scores)[i])),main=colnames(scores)[i] ,xlab="scores",ylab="calibrated scores")
  abline(a=0,b=1,col="green",lty=2)
}


####select 

######TEST

# predict class and then attach test class
#predictions<-as.data.frame(scores)
predictions<-as.data.frame(probs)

predictions$predict <- names(predictions)[1:5][apply(predictions[,1:5], 1, which.max)]
predictions$predict <-ys
predictions$observed <- y
head(predictions)

save(predictions,file="Predictions_FIS_Diff.RData")

#########################
##########################
###########################
# 1 ROC curve, mock vs non mock
roc.mock <- roc(ifelse(predictions$observed=="C_1", "C_1", "C_2"), as.numeric(predictions$C_1))
plot(roc.mock, col = "gray60",main="AUC internal predictions")

# others
roc.lethal <- roc(ifelse(predictions$observed=="C_2", "C_2", "C_1"), as.numeric(predictions$C_2))
roc.resist <- roc(ifelse(predictions$observed=="C_3", "C_3", "C_1"), as.numeric(predictions$C_3))


lines(roc.lethal, col = "blue")
lines(roc.resist, col = "red")


text(0.1,0.7,col="gray60",paste("AUC_class_1 = ",format(0.7337, digits=5, scientific=FALSE)))
text(0.1,0.65,col="blue",paste("AUC_class_2 = ",format(0.5696, digits=5, scientific=FALSE)))
text(0.1,0.6,col="red",paste("AUC_class_3 = ",format(0.7131, digits=5, scientific=FALSE)))

##get accuracy for each fold RF

#######build RF using best set of genes with low error

##loook at less error rate

which(errors==min(errors))


to.build<-to.train[rownames(to.train) %in% names(genes[[which(errors==min(errors))]]),]

rf.pred <- randomForest(
  t(to.build),
  classes,
  #mc=8,
  ntree=ntrees,
  #strata=y,
  mtry=sqrt(ncol(t(to.train))),
  #mtry=100,
  sampsize=rep(min(table(classes)),length(table(classes))),#Balancing by sampling stratification
  proximity=TRUE,
  oob.prox=TRUE,
  importance=TRUE,
  keep.inbag=TRUE,
  do.trace=FALSE,
  seed=seed
)


###get sens and spec



confusion<-rf.list[[3]]$confusion
sensitivity=(confusion[2,2]/(confusion[2,2]+confusion[2,1]))*100####for C_2
specificity=(confusion[1,1]/(confusion[1,1]+confusion[1,2]))*100
overall_error=rf.pred$err.rate[length(rf.pred$err.rate[,1]),1]*100
overall_accuracy=1-overall_error
class1_error=paste(rownames(confusion)[1]," error rate= ",confusion[1,6], sep="")
class2_error=paste(rownames(confusion)[2]," error rate= ",confusion[2,6], sep="")
overall_accuracy=100-overall_error


##Select LOW samples to predict


predictor_data<-t(lo.adj[rownames(lo.adj) %in% rownames(rf.list[[5]]$importance),])##



RF_predictions_probs<-predict(rf.list[[5]], predictor_data, type="p")





RF_predictions_responses<-predict(rf.list[[3]], predictor_data, type="response")



RF_p<-as.data.frame(RF_predictions_probs)
RF_a<-as.data.frame(RF_predictions_responses)


##calibrate with glmnet and give response

probs<-as.data.frame(probs)
write.csv(probs,"Predictor_balanced/Probs_calibrated_Balanced.csv")





##function to calculate CV in folds

##for gene selection


calcultateCVfold.genes <- function(badj,y,fold,p,cores,ntrees){
  
  
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
           ,mc=cores,ntree=ntrees,importance=TRUE,proximity=TRUE,
          oob.prox=TRUE,
         keep.inbag=TRUE,
        do.trace=FALSE) 
 
  message("predicting test set ...",Sys.time())
  
  rf.scores <- predict(rf,badj[fold$test,match(rownames(rf$importance),
                                               colnames(badj[fold$test,]))],type="prob")
  
  err <- sum(colnames(rf.scores)[apply(rf.scores,1,which.max)]!=y[fold$test])/length(fold$test)
  message("misclassification error: ",err)
  
  return(list(rf.scores,imp.meandecrease[or[1:p]],err,rf))
}





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
  
 # rf <- rfp(badj[fold$train,or[1:p]],y[fold$train],
  #          sampsize=rep(min(table(y[fold$train])),length(table(y[fold$train])))
   #         ,mc=cores,ntree=ntrees,importance=TRUE,proximity=TRUE,
    #        oob.prox=TRUE,
     #       keep.inbag=TRUE,
      #      do.trace=FALSE) 
  rf <- randomForest(badj[fold$train,or[1:p]],y[fold$train],
            sampsize=rep(min(table(y[fold$train])),length(table(y[fold$train])))
            ,mc=cores,ntree=ntrees,importance=TRUE,proximity=TRUE,
            oob.prox=TRUE,
            keep.inbag=TRUE,
            do.trace=FALSE) 
  
  
  message("predicting test set ...",Sys.time())
  
  rf.scores <- predict(rf,badj[fold$test,match(rownames(rf$importance),
                                               colnames(badj[fold$test,]))],type="prob")
  
  err <- sum(colnames(rf.scores)[apply(rf.scores,1,which.max)]!=y[fold$test])/length(fold$test)
  message("misclassification error: ",err)
  
  return(list(rf.scores,imp.meandecrease[or[1:p]],err))
}


###FUNCTIONS FOR VARIABLE SELECTION

#' Wrapper function to call random forests function.
#'
#' Provides an interface to different parallel implementations of the random
#' forest algorithm. Currently, only the \code{ranger} package is
#' supported.
#'
#' @param x matrix or data.frame of predictor variables with variables in
#'   columns and samples in rows (Note: missing values are not allowed).
#' @param y vector with values of phenotype variable (Note: will be converted to factor if
#'   classification mode is used).
#' @param ntree number of trees.
#' @param mtry.prop proportion of variables that should be used at each split.
#' @param nodesize.prop proportion of minimal number of samples in terminal
#'   nodes.
#' @param no.threads number of threads used for parallel execution.
#' @param method implementation to be used ("ranger").
#' @param type mode of prediction ("regression", "classification" or "probability").
#' @param ... further arguments needed for \code{\link[relVarId]{holdout.rf}} function only.
#'
#' @return An object of class \code{\link[ranger]{ranger}}.
#'
#' @import methods stats
#'
#' @export
#'
#' @examples
#' # simulate toy data set
#' data = simulation.data.cor(no.samples = 100, group.size = rep(10, 6), no.var.total = 200)
#'
#' # regression
#' wrapper.rf(x = data[, -1], y = data[, 1],
#'            type = "regression", method = "ranger")

wrapper.rf <- function(x, y, ntree = 500, mtry.prop = 0.2, nodesize.prop = 0.1, no.threads = 1,
                       method = "ranger", type = "regression", ...) {
  
  ## check data
  if (length(y) != nrow(x)) {
    stop("length of y and number of rows in x are different")
  }
  
  if (any(is.na(x))) {
    stop("missing values are not allowed")
  }
  
  if (type %in% c("probability", "regression") & (is.character(y) | is.factor(y))) {
    stop("only numeric y allowed for probability or regression mode")
  }
  
  ## set global parameters
  nodesize = floor(nodesize.prop * nrow(x))
  mtry = floor(mtry.prop * ncol(x))
  if (mtry == 0) mtry = 1
  
  if (type == "classification") {
    #    print("in classification")
    y = as.factor(y)
  }
  
  ## run RF
  if (method == "ranger") {
    if (type == "probability") {
      y = as.factor(y)
      prob = TRUE
    } else {
      prob = FALSE
    }
    
    rf = ranger::ranger(data = data.frame(y, x),
                        dependent.variable.name = "y",
                        probability = prob,
                        importance = "permutation", scale.permutation.importance = FALSE,
                        num.trees = ntree,
                        mtry = mtry,
                        min.node.size = nodesize,
                        num.threads = no.threads,
                        write.forest = TRUE,
                        ...)
  } else {
    stop(paste("method", method, "undefined. Use 'ranger'."))
  }
  
  return(rf)
}


#' Error calculation.
#'
#' Calculates errors by comparing predictions with the true values. For
#' regression and probability mode, it will give root mean squared error (rmse) and
#' pseudo R-squared (rsq). For classification mode, overall accuracy (acc), overall
#' error (err), Matthews correlation coefficient (mcc), sensitivity (sens) and
#' specificity (spec) are returned.
#'
#' @param rf Object of class \code{\link[ranger]{ranger}}
#' @param true vector with true value for each sample
#' @param test.set matrix or data.frame of predictor variables for test set with variables in
#'   columns and samples in rows (Note: missing values are not allowed)
#' @inheritParams wrapper.rf
#'
#' @return numeric vector with two elements for regression and probability estimation (rmse, rsq) and
#' five elements for classification (acc, err, mcc, sens, spec)
#'
#' @export
#'
#' @examples
#' # simulate toy data set
#' data = simulation.data.cor(no.samples = 100, group.size = rep(10, 6), no.var.total = 200)
#'
#' # random forest
#' rf = wrapper.rf(x = data[, -1], y = data[, 1],
#'                 type = "regression")
#'
#' # error
#' calculate.error(rf = rf, true = data[, 1])

calculate.error <- function(rf, true, test.set = NULL) {
  
  if (is(rf, "ranger")) {
    if (!is.null(test.set)) {
      pred = predict(rf, data = test.set)$predictions
    } else {
      pred = rf$predictions
    }
    if (rf$treetype == "Probability estimation") {
      pred = pred[, 2]
    }
  } else {
    stop(paste("rf needs to be of class ranger"))
  }
  
  if ((is(rf, "randomForest") && rf$type == "classification") |
      (is(rf, "ranger") && rf$treetype == "Classification")) {
    conf.matrix = table(pred = pred, true = true)
    tp = conf.matrix[2, 2]
    tn = conf.matrix[1, 1]
    fn = conf.matrix[2, 1]
    fp = conf.matrix[1, 2]
    
    ## accuracy
    acc = (tp + tn) / sum(conf.matrix)
    
    ## Matthews correlation coefficient
    mcc = (tp * tn - fp * fn) /
      sqrt( (tp + fn) * (tn + fp) * (tp + fp) * (tn + fn))
    
    ## sensitivity
    sens = tp / (tp + fn)
    
    ## specificity
    spec = tn / (fp + tn)
    
    error = c(err = 1 - acc, acc = acc, mcc = mcc, sens = sens, spec = spec)
  } else {
    mse = sum((pred - true)^2, na.rm = TRUE) / sum(!is.na(pred))
    
    ## pseudo R-squared uses sum of squared differences divided by n instead of variance!
    v = sum((true - mean(true))^2) / length(true)
    rsq = 1 - mse/v
    error = c(rmse = sqrt(mse), rsq = rsq)
  }
  
  return(error)
}


#' Get variable importance.
#'
#' Extracts variable importance depending on class of random forest object.
#'
#' @param rf Object of class \code{\link[ranger]{ranger}}
#'
#' @return numeric vector with importance value for each variable (in original order)
#'
#' @export

get.vim <- function(rf) {
  if (is(rf, "ranger")) {
    vim = ranger::importance(rf)
  } else {
    stop(paste("rf needs to be of class ranger"))
  }
  return(vim)
}


#' Variable selection using Boruta function.
#'
#' Variable selection using the Boruta function in the R package \code{\link[Boruta]{Boruta}}.
#'
#' This function selects only variables that are confirmed based on Boruta implementation.
#' For more details see \code{\link[Boruta]{Boruta}}.
#' Note that this function uses the ranger implementation for variable selection.
#'

#' @inheritParams wrapper.rf
#' @param pValue confidence level (default: 0.01 based on Boruta package)
#' @param maxRuns maximal number of importance source runs (default: 100 based on Boruta package)
#'
#' @return List with the following components:
#'   \itemize{
#'   \item \code{info} data.frame
#'   with information of each variable
#'   \itemize{
#'   \item run.x = original variable importance (VIM) in run x
#'   (includes min, mean and max of VIM of shadow variables)
#'   \item decision = Boruta decision (Confirmed, Rejected or Tentative)
#'   \item selected = variable has been selected
#'   }
#'   \item \code{var} vector of selected variables
#'   \item \code{info.shadow.var} data.frame with information about
#'   minimal, mean and maximal shadow variables of each run
#'   }
#'
#'  @examples
#' # simulate toy data set
#' data = simulation.data.cor(no.samples = 100, group.size = rep(10, 6), no.var.total = 200)
#'
#' # select variables
#' res = var.sel.boruta(x = data[, -1], y = data[, 1])
#' res$var
#'
#' @export

var.sel.boruta <- function(x, y, pValue = 0.01, maxRuns = 100,
                           ntree = 500, mtry.prop = 0.2, nodesize.prop = 0.1,
                           no.threads = 1, method = "ranger", type = "regression") {
  
  ## variable selection using Boruta function
  ## ----------------------------------------
  ## mtry.prop not used
  res.boruta = Boruta::Boruta(x = x, y = y,
                              pValue = pValue, maxRuns = maxRuns,
                              ntree = ntree, min.node.size = floor(nodesize.prop * nrow(x)),
                              num.threads = no.threads)
  
  ## select variables
  dec = res.boruta$finalDecision
  ind.sel = rep(0, ncol(x))
  ind.sel[dec == "Confirmed"] = 1
  info.sel = data.frame(decision = dec, selected = ind.sel)
  
  #   ## info about variables
  #   info.var = t(res.boruta$ImpHistory)
  #   colnames(info.var) = paste("run", 1:ncol(info.var), sep = ".")
  #   info.var = merge(info.var, info.sel, all.x = TRUE, by.x = "row.names", by.y = "row.names")
  #   rownames(info.var) = info.var[, 1]
  #   info.var = info.var[, -1]
  #   info.var = info.var[c(colnames(x), "shadowMax", "shadowMean", "shadowMin"), ]
  #
  #   ind.shadow = grep("shadow", rownames(info.var))
  #   return(list(info = info.var[-ind.shadow, ],
  #               var = sort(rownames(info.var)[info.var$selected == 1]),
  #               info.shadow.var = info.var[ind.shadow, -which(colnames(info.var) %in% c("decision", "selected"))]))
  
  ## info about variables
  info.var = t(res.boruta$ImpHistory)
  colnames(info.var) = paste("run", 1:ncol(info.var), sep = ".")
  info.shadow.var = info.var[grep("shadow", rownames(info.var)),]
  info.var = info.var[-grep("shadow", rownames(info.var)),]
  if (all.equal(rownames(info.var), rownames(info.sel))) {
    info.var = cbind(info.var, info.sel)
  } else {
    info.var = merge(info.var, info.sel, by.x = "row.names", by.y = "row.names")
  }
  
  return(list(info = info.var,
              var = sort(rownames(info.var)[info.var$selected == 1]),
              info.shadow.var = info.shadow.var))
  
}

#' Variable selection using Vita approach.
#'
#' This function calculates p-values based on the empirical null distribution from non-positive VIMs as
#' described in Janitza et al. (2015). Note that this function uses the \code{importance_pvalues} function in the R package
#' \code{\link[ranger]{ranger}}.

var.sel.vita <- function(x, y, p.t = 0.05,
                         ntree = 500, mtry.prop = 0.2, nodesize.prop = 0.1,
                         no.threads = 1, method = "ranger", type = "regression") {
  
  ## train holdout RFs
  res.holdout = holdout.rf(x = x, y = y,
                           ntree = ntree, mtry.prop = mtry.prop, nodesize.prop = nodesize.prop,
                           no.threads = no.threads, type = type)
  
  ## variable selection using importance_pvalues function
  res.janitza = ranger::importance_pvalues(x = res.holdout,
                                           method = "janitza",
                                           conf.level = 0.95)
  res.janitza = as.data.frame(res.janitza)
  colnames(res.janitza)[1] = "vim"
  
  ## select variables
  ind.sel = as.numeric(res.janitza$pvalue == 0 | res.janitza$pvalue < p.t)
  
  ## info about variables
  info = data.frame(res.janitza, selected = ind.sel)
  return(list(info = info, var = sort(rownames(info)[info$selected == 1])))
}

#' Helper function for variable selection using Vita approach.
#'
#' This function calculates a modified version of the permutation importance using two cross-validation folds (holdout folds)
#' as described in Janitza et al. (2015). Note that this function is a reimplementation of the \code{holdoutRF} function in the
#' R package \code{\link[ranger]{ranger}}.
#'
#' @param x matrix or data.frame of predictor variables with variables in
#'   columns and samples in rows (Note: missing values are not allowed).
#' @param y vector with values of phenotype variable (Note: will be converted to factor if
#'   classification mode is used).
#' @param ntree number of trees.
#' @param mtry.prop proportion of variables that should be used at each split.
#' @param nodesize.prop proportion of minimal number of samples in terminal
#'   nodes.
#' @param no.threads number of threads used for parallel execution.
#' @param type mode of prediction ("regression", "classification" or "probability").
#'
#' @return Hold-out random forests with variable importance
#'
#' @references
#' Janitza, S., Celik, E. & Boulesteix, A.-L., (2015). A computationally fast variable importance test for random forest for high dimensional data, Technical Report 185, University of Munich, https://epub.ub.uni-muenchen.de/25587.

holdout.rf <- function(x, y, ntree = 500, mtry.prop = 0.2, nodesize.prop = 0.1, no.threads = 1,
                       type = "regression") {
  
  ## define two cross-validation folds
  n = nrow(x)
  weights = rbinom(n, 1, 0.5)
  
  ## train two RFs
  res = list(rf1 = wrapper.rf(x = x, y = y,
                              ntree = ntree, mtry.prop = mtry.prop,
                              nodesize.prop = nodesize.prop, no.threads = no.threads,
                              method = "ranger", type = type,
                              case.weights = weights, replace = FALSE,
                              holdout = TRUE),
             rf2 = wrapper.rf(x = x, y = y,
                              ntree = ntree, mtry.prop = mtry.prop,
                              nodesize.prop = nodesize.prop, no.threads = no.threads,
                              method = "ranger", type = type,
                              case.weights = 1 - weights, replace = FALSE,
                              holdout = TRUE))
  
  ## calculate mean VIM
  res$variable.importance = (res$rf1$variable.importance +
                               res$rf2$variable.importance)/2
  res$treetype = res$rf1$treetype
  res$importance.mode = res$rf1$importance.mode
  class(res) = "holdoutRF"
  return(res)
}

######################################end of function definition

