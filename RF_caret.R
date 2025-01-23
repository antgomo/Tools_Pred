#library(RandomForestUtils)
setwd("/mnt/md127/Analysis_Projects/Classes_RNAseq_2019/")

##load lo samples to predict
load("Test_Predictor.RData")
source("RF_Utilities_Multiclass.R")
##load classes 

annon<-data.pred[[2]]

classes <- paste("C",annon$consensuscluster,sep="_")
classes<-as.factor(classes)

data.tp<-data.pred[[1]]
###load data to test
data.tp<-data.tp[,match(rownames(annon),colnames(data.tp))]
input_features<-t(data.tp)
data.pred<-list(data.tp,annon,lo.adj)

##save(data.pred,file="Test_Predictor.RData")

#recode to case and control

library(forcats)


classes <- fct_recode(classes, Case="Patient",
                      Control="Healthy")
head(classes)


SAVE_PATH <- "./"
set.seed(1995)
seeds <- sample.int(10000000, 10)
RandomForestUtils::set_cores(10)

rf_results <- Run_RF_Pipeline(feature_table = input_features,
                                                 classes=classes,
                                                 metric = "ROC",
                                                 sampling=NULL,
                                                 repeats=10,
                                                 path=SAVE_PATH,
                                                 nmtry=4,
                                                 ntree=501,
                                                 nfolds=3,
                                                 ncrossrepeats = 5,
                                                 pro=0.8,
                                                 list_of_seeds = seeds)

boxplot(rf_results[[1]], rf_results[[2]], names=c("Cross Validation AUC", "Testing AUC"))


##Now we will re-run the pipeline but scrambling the assignment of case and control for each sample. This should give us a good understanding of how well our model would work by "chance".
scrambled_classes <- list()
set.seed(1880)
for(i in 1:10){
  scrambled_classes[[i]] <- sample(classes)
}

SAVE_PATH <- "./"
random_resullts <- get_random_rf_results(feature_table = input_features,
                                         list_of_scrambles = scrambled_classes,
                                         metric = "ROC",
                                         sampling = NULL,
                                         repeats = 10,
                                         path=SAVE_PATH,
                                         nmtry=4,
                                         ntree=1001,
                                         nfolds = 3,
                                         ncrossrepeats = 5,
                                         pro = 0.8,
                                         list_of_seeds = seeds)



##Now we can compare the test AUCs from the model trained on the real class labels and the model trained on scrambled train labels. This will gives us a reallly good idea about how well our model is preforming. 


boxplot(rf_results[[2]], random_resullts[[2]], names=c("Test AUROC Real", "Test AUROOC Random"))

#We can see that our model performs significantly better then a model trained on random labellings. 

#We can also use scripts in this package to generate a AUROC for each test and training data splits. 
#Note that the paramter labels is a dataframe that contains the sample names as rows and a column that 
#indicates the class of each sample as either "Case" or "Control"

```{r}
ROC_plot <- generate_ROC_curve(RF_models = rf_results[[5]], dataset=input_features, labels=clean_metadata, title = "AUROC of Enteric Disease Classifcation")
ROC_plot

##The last thing that we can look at is feature importance. This will give us an idea of what genera in the dataset are most predictive. 

Feature_importance_scores <- Calc_mean_accuray_decrease(rf_results[[4]])

knitr::kable(Feature_importance_scores[, c("Mean", "SD", "Min", "Max")])


#Multi-Class Summary Function
#Based on caret:::twoClassSummary
# From: http://moderntoolmaking.blogspot.com/2012/07/error-metrics-for-multi-class-problems.html

# RES: disable compilation for debugging
# require(compiler)
# multiClassSummary <- cmpfun(function (data, lev = NULL, model = NULL){
multiClassSummary <- function (data, lev = NULL, model = NULL){
  
  #Check data
  if (!all(levels(data[, "pred"]) == levels(data[, "obs"]))) 
    stop("levels of observed and predicted data do not match")
  
  ## Overall multinomial loss
  lloss <- mnLogLoss(data = data, lev = lev, model = model)
  
  #Calculate custom one-vs-all ROC curves for each class
  prob_stats <- lapply(levels(data[, "pred"]), 
                       function(class){
                         #Grab one-vs-all data for the class
                         pred <- ifelse(data[, "pred"] == class, 1, 0)
                         obs  <- ifelse(data[,  "obs"] == class, 1, 0)
                         prob <- data[,class]
                         
                         #Calculate one-vs-all AUC
                         prob_stats <- as.vector(ModelMetrics::auc(obs, prob))
                         names(prob_stats) <- c('ROC')
                         return(prob_stats) 
                       })
  roc_stats <- mean(unlist(prob_stats))
  
  #Calculate confusion matrix-based statistics
  CM <- confusionMatrix(data[, "pred"], data[, "obs"])
  
  #Aggregate and average class-wise stats
  #Todo: add weights
  # RES: support two classes here as well
  #browser() # Debug
  if (length(levels(data[, "pred"])) == 2) {
    class_stats <- c(CM$byClass, roc_vals)
  } else {
    class_stats <- colMeans(CM$byClass)
    names(class_stats) <- paste("Mean", names(class_stats))
  }
  
  # Aggregate overall stats
  overall_stats <- c(CM$overall, lloss, Mean_ROC = roc_stats)
  
  # Combine overall with class-wise stats and remove some stats we don't want 
  stats <- c(overall_stats, class_stats)
  stats <- stats[! names(stats) %in% c('AccuracyNull', "AccuracyLower", "AccuracyUpper",
                                       "AccuracyPValue", "McnemarPValue", 
                                       'Mean Prevalence', 'Mean Detection Prevalence')]
  
  # Clean names
  names(stats) <- gsub('[[:blank:]]+', '_', names(stats))
  
  if (length(levels(data[, "pred"]) == 2)) {
    # Change name ordering to place most useful first
    # May want to remove some of these eventually
    stats <- stats[c("logLoss", "Mean_ROC", 
                     "Mean_Sensitivity", "Mean_Specificity", "Accuracy", "Kappa", 
                     "Mean_Pos_Pred_Value", "Mean_Neg_Pred_Value", "Mean_Detection_Rate",
                     "Mean_Balanced_Accuracy")]
  }
  
  return(stats)
}