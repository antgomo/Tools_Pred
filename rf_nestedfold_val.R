
#######################From here, same pocedure as Nature, 2018, 5 fold internal cross validation

#nestedCVfold(to.train,classes,250,8,999,5,my.options[i])

badj<-to.train
y<-classes
p<-250
cores<-8
ntrees<-999
folds<-5
m<-"DEG"

nestedCVfold <- function(badj,y,p,cores,ntrees,folds,m){
  
  

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
  
  
  load(paste0("CVfold.",m,i,".",0,".RData"))
  
  message("calibrating raw scores fold ",i," ...",Sys.time())
  probs <- predict(cv.calfit$glmnet.fit,newx=rf.scores,type="response"
                   ,s=cv.calfit$lambda.1se)[,,1] # use lambda estimated by 10fold CVlambda
  
  
  err <- sum(colnames(probs)[apply(probs,1,which.max)] != anno[nfolds[[i]][[1]][[1]]$test])/length(nfolds[[i]][[1]][[1]]$test)
  
  message("misclassification error: ",err)
  
  save(probs,file=paste0("probsCVfold.",m,i,".",0,".RData"))
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

save(paste0(cv.calfit,file="calfit.",m,".RData"))


save(paste0(scores,probs,y,ys,yp,file="CVresults.",m,".RData"))##calibrate final predictions


###try to plot results

#rm(list=ls())
#gc()

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



######################################end of function definition

