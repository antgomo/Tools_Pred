library(randomForest)
library(parallel)
library(ggplot2)


##parallelized version of RF

rfp <- function(xx, ..., ntree = ntree, mc = mc, seed = NULL)
{
  if(!is.null(seed)) set.seed(seed, "L'Ecuyer")
  rfwrap <- function(ntree, xx, ...) randomForest::randomForest(x=xx,ntree=ntree,norm.votes=FALSE, ...)
  rfpar <- mclapply(rep(ceiling(ntree/mc),mc),mc.cores=mc, rfwrap, xx=xx, ...)
  do.call(randomForest::combine, rfpar)
}


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
  
  return(rf.scores)
}