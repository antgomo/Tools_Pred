library(msgl)

##Load data containing N samples and p features (covariates):
  
  #x <- # load design matrix (of size N x p)
  #classes <- # load class labels (a vector of size N)
    
    data(PrimaryCancers)
    
    ##build test set
    
    idx <- 1:10
    x.test <- x[idx,]
    x <- x[-idx,]
    classes.test <- classes[idx]
    classes <- classes[-idx]
    
    ##3. Estimate error using cross validation
    
   # Choose lambda (fraction of lambda.max) and alpha, with alpha = 1 for lasso, alpha = 0 for group lasso and alpha in the range (0,1) for sparse group lasso.
    
    #Use msgl::cv to estimate the error for each lambda in a sequence decreasing from the data derived lambda.max to lambda * lambda.max.
    #Lambda.max is the lambda at which the first penalized parameter becomes non-zero. A smaller lambda will take longer to fit and include more features.
    #The following code will run a 10 fold cross validation for each lambda value in the lambda sequence using 2 parallel units (using the foreach and doParallel packages.
    
    cl <- makeCluster(2)
    registerDoParallel(cl)
    
    fit.cv <- msgl::cv(x, classes, fold = 10, alpha = 0.5, lambda = 0.1, use_parallel = TRUE)
    
    stopCluster(cl)
    
    ##We have now cross validated the models corresponding to the lambda values, one model for each lambda value. We can summarize the validation as follows.
    
    fit.cv
    
    #Hence, the best model is obtained using lambda index 95 and it has a cross validation error of 0.12. 
    #The expected number of selected features is 53.5 and the expected number of parameters is 276.9.
    
    #Fit the final model
    
    fit <- msgl::fit(x, classes, alpha = 0.5, lambda = 0.1)
    
    ##As we saw in the previous step the model with index 79 had the best cross validation error, we may take a look at the included features using the command:
    features(fit)[[best_model(fit.cv)]] # Non-zero features in best model
   ### Hence 48 features are included in the model, this is close to the expected number based on the cross validation estimate.
    
  ##  The sparsity structure of the parameters belonging to these 48 features may be viewed using    
    
    parameters(fit)[[best_model(fit.cv)]]
    
   ## We may also take a look at the estimate parameters (or coefficients)
    
    coef(fit, best_model(fit.cv))[,1:5] # First 5 non-zero parameters of best model
    
    ## If we count the total number of non-zero parameters in the model we get in this case 250, which is close to the expected based on the cross validation estimate.
    
    ## Use your model for predictions
    
    Use the final model to predict the classes of the M samples in x.test.
    
    res <- predict(fit, x.test)
    
    res$classes[,best_model(fit.cv)] # Classes predicted by best model
    
    classes.test # True classes
    
    ##We may also get the estimated probabilities for each of the classes
    
    res$response[[best_model(fit.cv)]]
      