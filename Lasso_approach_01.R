library(msgl)
setwd("/Users/agomez/IMIDOMICS/Predictor/")
library(limma)
library(minfi)
library(DMRcate)

setwd("/Users/agomez/IMIDOMICS/Predictor/")
library(gplots)
library(RColorBrewer)
library(M3C)
library(NMF) # loading for heatmap plotting function
library(ggsci)
library(heatmap.2x)## approach
library(parallel)
library(ggplot2)
library(glmnet)
library(doParallel)
library(HandTill2001)
library(pROC)


######################################end of function definition
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(doParallel)

####load all samples FIS and P6

idat.folder <- "/Users/agomez/IMIDOMICS/Predictor/" ###change for wour path dir
targets <- read.metharray.sheet(base=idat.folder)

registerDoParallel(cores = 8)

###loading data

#rgset <- read.metharray.exp(targets = targets,force=TRUE)

load("Predictor_ALLCelltypes.RData")
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
design <- model.matrix(~targets$RESP+targets$Age + targets$SEX)##only wk0



####generate adjusted beta matrix
fit <- lmFit(betas, design)
mAdj.fit    <- fit$coefficients[,-c(1:2)]

mAdj<- as.matrix(betas - mAdj.fit %*% t(design[,-c(1:2)]))

betaAdj     <- ilogit2(mAdj)

rm(fit)
rm(mAdj)
rm(mAdj.fit)
rm(betas)
gc()

colnames(betaAdj)<-gsub("-","",targets$Code)

##try to find batches in data
##same without batch
pca <- prcomp(t(betaAdj))
pca_data_perc<-round(100*((pca$sdev)^2/sum(pca$sdev^2)),digits=2)



df_pca_data<-data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], batch=targets$Project)

ggplot(df_pca_data, aes(PC1,PC2, color = batch,label=rownames(df_pca_data)))+
  geom_point(size=.11)+geom_text(label=rownames(df_pca_data))
  labs(x=paste0("PC1 (",pca_data_perc[1],")"), y=paste0("PC2 (",pca_data_perc[2],")"))

  
  ##remove outliers
  
  betaAdj<-betaAdj[,-c(46,133,106,32,5)]
  betas<-betas[,-5]
  
 targets<-targets[-5,]
 targets$COD<-gsub("-","",targets$Code)
##split in test set and training set

training.set<-targets[targets$Project=="FIS",]
test.set<-targets[targets$Project=="P6",]


rm(pca)
rm(df_pca_data)
gc()
##filter outliers

##filter by annotation


ann850k<-getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann850k<-as.data.frame(ann850k)
## Outputs
proms<-ann850k[grep("TSS1500|TSS200|5'UTR|1stExon",ann850k$UCSC_RefGene_Group),]
proms.island<-proms[proms$Relation_to_Island=="Island",]
proms.island<-proms.island[!(proms.island$UCSC_RefGene_Name==" "|is.na(proms.island$UCSC_RefGene_Name)),]

rm(ann850k)
gc()





###########RF procedure to test classes

##using selected genes diff express

to.train<-betaAdj[rownames(betaAdj) %in% rownames(proms.island),colnames(betaAdj) %in% training.set$COD]##adjust cellcomp

to.train<-to.train[,match(training.set$COD,colnames(to.train))]
#to.train<-betaAdj[,colnames(betaAdj) %in% training.set$COD]##adjust cellcomp

CV<- apply(to.train,1,sd)/(rowMeans(to.train) )###
sel<-CV[CV>.1]




##Load data containing N samples and p features (covariates):
  
  #x <- # load design matrix (of size N x p)
  #classes <- # load class labels (a vector of size N)
    
classes.test<-training.set$RESP
names(classes.test)<-training.set$COD

    
    ##3. Estimate error using cross validation
    
   # Choose lambda (fraction of lambda.max) and alpha, with alpha = 1 for lasso, alpha = 0 for group lasso and alpha in the range (0,1) for sparse group lasso.
    
    #Use msgl::cv to estimate the error for each lambda in a sequence decreasing from the data derived lambda.max to lambda * lambda.max.
    #Lambda.max is the lambda at which the first penalized parameter becomes non-zero. A smaller lambda will take longer to fit and include more features.
    #The following code will run a 10 fold cross validation for each lambda value in the lambda sequence using 2 parallel units (using the foreach and doParallel packages.
    
    cl <- makeCluster(3)
    registerDoParallel(cl)
    
    fit.cv <- msgl::cv(t(to.train[names(sel),]), classes.test, fold = 5, alpha = 0.5, lambda = 0.1, use_parallel = TRUE)
    
    stopCluster(cl)
    
    ##We have now cross validated the models corresponding to the lambda values, one model for each lambda value. We can summarize the validation as follows.
    
    fit.cv
    
    #Hence, the best model is obtained using lambda index 95 and it has a cross validation error of 0.12. 
    #The expected number of selected features is 53.5 and the expected number of parameters is 276.9.
    
    #Fit the final model
    
    fit <- msgl::fit(t(to.train[names(sel),]), classes.test, alpha = 0.5, lambda = 0.1)
    
    ##As we saw in the previous step the model with index 79 had the best cross validation error, we may take a look at the included features using the command:
    features(fit)[[best_model(fit.cv)]] # Non-zero features in best model
   ### Hence 48 features are included in the model, this is close to the expected number based on the cross validation estimate.
    
  ##  The sparsity structure of the parameters belonging to these 48 features may be viewed using    
    
    parameters(fit)[fit.cv[[40]]]
    
   ## We may also take a look at the estimate parameters (or coefficients)
    
    coef(fit, 40) # First 5 non-zero parameters of best model
    
    ## If we count the total number of non-zero parameters in the model we get in this case 250, which is close to the expected based on the cross validation estimate.
    
    ## Use your model for predictions
    
  #  Use the final model to predict the classes of the M samples in x.test.
    predictor_data<-t(betaAdj[,colnames(betaAdj) %in% test.set$Name])##
    
    res <- predict(fit, to.train)
    
    res$classes[,best_model(fit.cv)] # Classes predicted by best model
    
    classes.test # True classes
    
    ##We may also get the estimated probabilities for each of the classes
    
    res$response[[best_model(fit.cv)]]
      