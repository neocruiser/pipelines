pkgs <- c('RColorBrewer', 'pvclust', 'gplots', 'vegan',
          'dplyr', 'mRMRe', 'glmnet', 'caret', 'foreach',
          'doSNOW', 'lattice', 'ROCR', 'earth')
lapply(pkgs, require, character.only = TRUE)

## logging
source("./script/lsos.R")

### DEFINE DATA FITTING MODEL
classification=TRUE
regression=FALSE

## define feature structure
grouped=TRUE
binomial=FALSE


##########################
## Define new functions ##
##########################
## Classification
modelTune.clas <- function(dat, train, method, folds=10, rep=5, tune=10, grid=TRUE){
    ## requires caret
    ## GRID search HYPERPARAMETERS tuning
    ## Train model and test on independant dataset
    ## returns a classification error
    trainCtrl <- trainControl(method="repeatedcv",
                              number=folds,
                              repeats=rep,
#                              classProbs=T,
                              summaryFunction=defaultSummary)
    
    ## Choosing the right Hyper-parameters. GRID ANALYSIS
    if ( grid == TRUE ) {
        
        # Tune hyper-parameters
        if ( method == "svmLinear" ) {
            ## support vector machines with linear kernel
            grid_models <- expand.grid(.C=seq(1, 20, length=40))
        } else if ( method == "svmPoly") {
            ## svmRadial
            grid_models <- expand.grid(.degree=seq(0,20,.5),.scale=10^seq(-1,-3,length=15),.C=seq(1:5))
        } else if ( method == "svmRadialSigma") {
            ## svm with radial basis function kernel
            grid_models <- expand.grid(.C=seq(1:20), .sigma=10^seq(-1,-3,length=40))
        } else if ( method == "svmLinear3") {
            ## l2 regularized support vector machine (dual) with linear kernel
            grid_models <- expand.grid(.cost=seq(0.1,1,0.01), .Loss=seq(0:7))
        } else if ( method == "lda2" ) {
            ## Linear discriminant analysis
            grid_models <- expand.grid(.dimen=seq(1,15,.25))
        } else if ( method == "bagFDA" || method == "fda" ) {
            ## bagged flexible discriminant analysis
            grid_models <- expand.grid(.degree=seq(.1,2,length=10), .nprune=seq(1,50,length=15))
        } else if ( method == "pda" ) {
            ## penalized discriminant analysis
            grid_models <- expand.grid(.lambda=10^seq(-0.05,-4,length=35))
        } else if ( method == "loclda" ) {
            ## localized linear discriminant analysis
            grid_models <- expand.grid(.k=seq(1,400,length=40))
        } else if ( method == "bagFDAGCV" ) {
            ## bagged FDA using gCV pruning
            grid_models <- expand.grid(.degree=seq(.1,2,length=40))
        } else if ( method == "C5.0" ) {
            ## c5 decision trees
            grid_models <- expand.grid(.trials=seq(1,100,5), .model=c("rules","tree"), .winnow=c(TRUE,FALSE))
        } else if ( method == "LogitBoost" ) {
            ## Boosted logistic regression
            grid_models <- expand.grid(.nIter=seq(1,100,length=40))
        } else if ( method == "regLogistic" ) {
            ## regularized logistic regression
            grid_models <- expand.grid(.cost=seq(.1,5,length=10),.loss=c(0:7),.epsilon=seq(-.5,-3,length=10))
        } else if ( method == "kernelpls" ) {
            ## partial least squares
            grid_models <- expand.grid(.ncomp=seq(0.1,20,length=40))
        } else if ( method == "multinom" ) {
            ## penalized multinomial regression
            grid_models <- expand.grid(.decay=10^seq(-.5,-5,length=40))
        } else if ( method == "nnet" ) {
            ## neural networks
            grid_models <- expand.grid(.size=seq(1,20,length=40), .decay=10^seq(-.5,-5,length=40))
        } else if ( method == "monmlp" ) {
            ## monotone multi-layer perceptron neural network
            grid_models <- expand.grid(.hidden1=seq(1:5),.n.ensemble=c(10))
        } else if ( method == "dnn" ) {
            ## stacked autoencoder deep neural network
            grid_models <- expand.grid(.layer1=seq(1:10),.layer2=seq(1:5),.layer3=seq(1:2),
                                       .hidden_dropout=seq(.1,2,length=5),
                                       .visible_dropout=seq(.1,2,length=2))
        } else if ( method == "rf" ) {
            ## random forest
            grid_models <- expand.grid(.mtry=seq(1:15))
        } else if ( method == "RRF" ) {
            ## regularized random forest
            grid_models <- expand.grid(.mtry=seq(1:10),.coefreg=10^seq(-1,-3,length=10),.coefimp=10^seq(-1,-2,length=5))
        } else if ( method == "kknn" ) {
            ## weighted k nearest neighbors
            grid_models <- expand.grid(.kmax=seq(1:15),.distance=seq(1:5),.kernel=c("optimal","rank","gaussian",
                                                                                    "inv","cos","rectangular",
                                                                                    "triangular","biweight",
                                                                                    "triweight"))
        } else if ( method == "naive_bayes" ) {
            ## naive bayes
            grid_models <- expand.grid(.laplace=seq(0,2,length=5),.usekernel=c(FALSE,TRUE),.adjust=c(0,1,2))
        } else if ( method == "gbm" ) {
            ## stochastic gradient boosting
            grid_models <- expand.grid(.n.trees=seq(5,300,length=20),.interaction.depth=seq(1,5,length=10),
                                       .shrinkage=10^seq(-.1,-2,length=10),.n.minobsinnode=seq(1:10))
        }
        
        # train the model
        lapsed <- system.time(modelTrain <- train(y~., data=dat[train,],
                                                  method=method,
                                                  trControl= trainCtrl,
                                                  preProc=c("center","scale"),
                                                  tuneGrid=grid_models,
                                                  tuneLength=tune))

    } else if ( grid == FALSE ) {
        ## do not tune the hyperparameters
        ## just use the default options for each model
        lapsed <- system.time(modelTrain <- train(y~., data=dat[train,],
                                                  method=method,
                                                  trControl= trainCtrl,
                                                  preProc=c("center","scale")))
    }

    pdf("model.train.pdf")
    plot(modelTrain)
    dev.off()

    ## get performance scores from out-of-bag observations
    ## after cross validation during hyperparameter tuning
    results <- modelTrain$results

    ## get prediction accuracy
    ## confusion matrix for classification
    Predd <- predict(modelTrain, newdata=dat[-train,], type="raw")
    conf.m <- confusionMatrix(data=Predd, dat[-train,1])

    ## Compile everything for later accessibility
    ## time, best model based on training performance metrics
    ## tuning accuracy per parameter
    output <- list(timeLapsed=lapsed,
                   bestModel=modelTrain,
                   Results=results,
                   Hyperparameters=modelTrain$bestTune,
                   ConfusionMatrix=conf.m)
    return(output)
}




## Regression without model tuning
model.reg <- function(dat,train,method,folds=10,rep=5,tune){
    ## requires caret
    ## Train model and test on independant dataset
    ## returns RMSE
    trainCtrl <- trainControl(method="repeatedcv",
                              number=folds,
                              repeats=rep)
    lapsed <- system.time(modelTrain <- train(y~., data=dat[train,],
                                              method=method,
                                              trControl= trainCtrl,
                                              preProc=c("center","scale"),
                                              tuneLength=tune ))
    ploted <- plot(modelTrain)
    Predd <- predict(modelTrain, newdata=dat[-train,], type="raw")
    ## Test set MSE for regression
    rmse <- mean((Predd - y[-train])^2)
    output <- list(ploted,TimeLapsed=lapsed,Prediction.Estimates=Predd,Hyperparameters=modelTrain$bestTune, RMSE=rmse)
    return(output)
}

## Regression
modelTune.reg <- function(dat,train,method,folds=10,rep=5,tune,ctl){
    ## requires caret
    ## Uses GRID for HYPERPARAMETER tuning
    ## Train model and test on independant dataset
    ## returns RMSE
    trainCtrl <- trainControl(method="repeatedcv",number=folds, repeats=rep)	## Regression
    lapsed <- system.time(modelTrain <- train(y~., data=dat[train,],
                                              method=method,
                                              trControl= trainCtrl,
                                              preProc=c("center","scale"),
                                              tuneGrid=ctl,
                                              tuneLength=tune ))
    ploted <- plot(modelTrain)
    Predd <- predict(modelTrain, newdata=dat[-train,], type="raw")
    ## Test set MSE for regression
    rmse <- mean((Predd - y[-train])^2)
    output <- list(ploted,TimeLapsed=lapsed,Prediction.Estimates=Predd, Hyperparameters=modelTrain$bestTune, RMSE=rmse)
    return(output)
}

ensemble.mean <- function(a,b){
    ## Ensemble Methods
    ## calculates RMSE of joint predictions
    ## Weighted averaging of 2 base learners
    E.pred1 <- (a[[2]]+b[[2]])/2
    E.pred2 <- (a[[2]]*2+b[[2]])/3
    E.pred3 <- (a[[2]]+b[[2]]*2)/3
    M1 <- mean((E.pred1 - y[test])^2)		## Test set MSE for regression
    M2 <- mean((E.pred2 - y[test])^2)		## Test set MSE for regression
    M3 <- mean((E.pred3 - y[test])^2)		## Test set MSE for regression
    ploted <- plot(y=c(M1,M2,M3),x=1:3, lty=5,cex=1,pch=21:23,type="b",bg="red")
    output <- list(ploted,model1.ab=M1,model2.2ab=M2,model3.a2b=M3)
    return(output)
}



##########################
## Load expression data ##
##########################
# optimized for t-statistics microarray expressions
# rows=genes
# col=samples
cat("\n\nNormalized expression scores: Samples are columns and genes are rows\n")
means <- read.table("expressions", sep="\t", header=T, row.names=1)
dim(means)
cat("\n\nStandardized transformed scores: Genes are columns and samples are rows\n")
adj.x <- t(decostand(means, "standardize"))
dim(adj.x)


## restructure the dataset
metadata <- read.table("summary/phenodata", sep = "\t", header = T) %>%
    dplyr::select(SAMPLE_ID, Timepoint, GROUP, SITE, Score, Prediction, ABClikelihood) %>%
    filter(Timepoint != "T2") %>%
    mutate(Groups = case_when(GROUP %in% c("CNS_RELAPSE_RCHOP",
                                            "CNS_RELAPSE_CHOPorEQUIVALENT",
                                            "CNS_DIAGNOSIS") ~ "CNS",
                               GROUP %in% c("TESTICULAR_NO_CNS_RELAPSE", "NO_RELAPSE") ~ "NOREL",
                               GROUP == "SYTEMIC_RELAPSE_NO_CNS" ~ "SYST",
                               TRUE ~ "CTRL")) %>%
    mutate(ABClassify = case_when(ABClikelihood >= .9 ~ "ABC",
                                  ABClikelihood <= .1 ~ "GCB",
                                  TRUE ~ "U")) %>%
    mutate(ABCScore = case_when(Score > 2412 ~ "ABC",
                                Score <= 1900 ~ "GCB",
#                                Score == NA ~ "NA",
                                TRUE ~ "U")) %>%
    mutate(Nodes = case_when(SITE == "LN" ~ "LN",
                             SITE == "TO" ~ "LN",
                             SITE == "SP" ~ "LN",
                             TRUE ~ "EN")) %>%
    mutate(Lymphnodes = case_when(Nodes == "LN" ~ 1, TRUE ~ 0))

# factorize
metadata$Groups <- as.factor(metadata$Groups)
metadata$ABClassify <- as.factor(metadata$ABClassify)
metadata$ABCScore <- as.factor(metadata$ABCScore)
metadata$Nodes <- as.factor(metadata$Nodes)
metadata$Lymphnodes <- as.factor(metadata$Lymphnodes)


# choose samples strutcture from metadata
y <- metadata$Groups

if ( nlevels(y) > 2 ){
    response="multinomial"
} else if ( nlevels(y) == 2 ){
    response="bionomial"
} else {
    response="mgaussian"
}


## prepare testing dataset
# set seed for reproducibility
ed <- floor(abs(rnorm(1) * 10000000))
set.seed(ed)

# Split the dataset into 80% training data
training <- sample(1:nrow(adj.x), nrow(adj.x)/1.25)
tr <- length(training)


## EXPERIMENTAL
# create dummy variables or multi level contrasts
# first level (the baseline) is rolled into the intercept
# all other levels have a coefficient that differ from the baseline
# Helmert regressors compare each level with the average of the preceding ones
# first coefficient is the mean of the first two levels minus the first level
# second coefficient is the mean of all three levels minus the mean of the first two levels

associations=y
contrasts(associations) <- "contr.helmert"
contrasts(associations)

associations <- model.matrix(~0 + y)
colnames(associations) <- levels(y)




########################
## FEATURE EXTRACTION ##
########################
## Result: Lambda for Lasso feature selection with highest accuracy outcome
# feature extraction
# |-- randomize seed
# |   |-- iterate multiple seeds
# |
# |-- parameter tuning
# |   |-- grid tuning
# |   |-- nested cross validation
# |
# |-- fit linear model
# |   |-- first on training set
# |   |-- second on testing set
# |   |-- iterate multiple model fitting
# |
# |-- get prediction scores for each iteration
# |-- plot accuracy scores
# |   |-- iterate multiple accuracy tests
# |   |-- ROC curves for Cross validation test
# |   |-- ROC curves for prediction/validation set
# |
# |-- create summary of all iterations
# |
# |-- get the best lambda
# |   |-- best accuracy has the best lambda
# |   |-- use this lambda for subsequent analyses
# |   |-- plot final accuracy at this lambda cutoff


# WARNING: sometimes LOGNET throws an error of 0 or 1 observations
# this is due to the unbalanced nature of cross validation
# SOLUTION: the while condition will repeat the test until success
success=FALSE
iterations=5

while (success == FALSE) {
    pdf(paste0("cvROC.iterations",iterations,".",response,".pdf"))
    couleurs <- brewer.pal(nlevels(y), name = 'Dark2')
    dm=NULL
    df=NULL
    ## The control class is throwing an error
    # prediction cannot be done on small sample size
    # remove CTRL class
    nl <- c(1:nlevels(y))[levels(y)!="CTRL"]

    if (!require(ROCR)) {
        stop("Can't continue can't load ROCR")
    }

    for (i in nl) {
        for (e in 1:iterations) {
            # set seed for reproducibility
            ede <- floor(abs(rnorm(1) * 10000000))
            set.seed(ede)

            # Split the dataset into 80% training data
            training <- sample(1:nrow(adj.x), nrow(adj.x)/1.25)
            tr <- length(training)

            # make sure all sample categories are included
            # in the training and testing sets
            # if not, errors occur during training
            # for missing observations in certain classes
            while (
            (length(unique(y[training])) != nlevels(y))
            &
            (length(unique(y[-training])) != nlevels(y))
            ) {
                ede <- floor(abs(rnorm(1) * 10000000))
                set.seed(ede)
                # split 70% of the data
                randomize <- sample(nrow(adj.x))
                training <- sample(randomize, nrow(adj.x)/1.429)
                #table(y[-training])
                tr <- length(training)
            }

            # fit a generalized linear model via penalized maximum likelihood
            # if alpha 1 then lasso L1 penality and discard genes
            # if alpha 0 then ridge regression then L2 and rank genes
            if ( grouped == TRUE ){
                index="grouped"
                ncv=10
                setalpha=1
                # fitting a symmetric multinomial model,
                grid <- 10^seq(5, -5, length=100)
                lasso.trained <- glmnet(adj.x[training,],
                                        y[training],
                                        alpha=setalpha,
                                        lambda=grid,
                                        family = response,
                                        standardize=F,
                                        type.multinomial=index)
            } else {
                stop("Alternative feature indexes are not yet implemented")
            }

            # Cross validation for hyperparameter tuning.
            # Optimization of model selection to avoid overfitting
            # also, silence errors due to unbalanced observations
            cv.out <- try(cv.glmnet(adj.x[training,],
                                y[training],
                                alpha=setalpha,
                                family = response,
                                standardize=F,
                                nfolds = ncv,
                                type.multinomial=index),
                          silent = TRUE)

            # make sure CV did not encounter unbalanced observations
            # if unbalanced obs exist, the whole iteration will repeat
            if (class(cv.out) != "cv.glmnet") {
                success=FALSE
                ## initiate iteration
                # not working, must be set early
                #if ( e > 1 ) { e = e - 1 } 
                break
            } else {
                success=TRUE
            }

            # get best lambda
            # lowest probability to overfit
            bestlam <- cv.out$lambda.min

            # get probabilities
            lasso.response <- predict(lasso.trained, s=bestlam, newx=adj.x[-training,], type="response")

            ## create list of iteration scores for each class
            probability.scores <- as.data.frame(lasso.response)[, i]
            dummy.labels <- as.vector(model.matrix(~0 + y[-training])[, i])

            # create multi dimensional list across all iterations
            if ( e == 1 ){
                ps <- list(probability.scores)
                dl <- list(dummy.labels)
            } else {
                ps <- c(ps, list(probability.scores))
                dl <- c(dl, list(dummy.labels))
            }

            # rename listing headers
            names(ps)[e]=paste0(levels(y)[i],"-class",i,".iteration",e)
            names(dl)[e]=paste0(levels(y)[i],"-class",i,".iteration",e)
            
            # create summary of probabilities
            # from a list of genes that have a mean probability to classify
            # all cases into correct patient diagnosis
            dm <- rbind(dm, data.frame(iterations=e,
                                       class=levels(y)[i],
                                       seed=ede,
                                       lambda=bestlam,
                                       probabilityScore=mean(ps[[e]]) ))

            ## compile accuracies into summary file
            # get sample labels
            lasso.labels <- predict(lasso.trained, s=bestlam, newx=adj.x[-training,], type="class")

            # get gene coefficients at selected lambda
            lasso.coef <- predict(lasso.trained, s=bestlam, type = "coefficients")

            ## get non-zero genes
            selected.genes <- lasso.coef[[1]]@i[ lasso.coef[[1]]@i >= 1]
            original.genes <- colnames(adj.x)
            selected.final <- original.genes[selected.genes]
            len <- length(selected.genes)

            ## build classification confusion-rate matrix
            if ( classification == TRUE ) {
                tab <- table(lasso.labels, y[-training])
                freq <- as.data.frame.matrix(tab)

                for ( q in 1:nrow(freq) ) {
                    ## prepare a table summary
                    ## of regularization accuracy
                    n=names(rowSums(freq)[q])
                    f <- freq[n,n] / rowSums(freq)[[q]] * 100

                    if ( setalpha == 1 ) {me="lasso"} else {me="ridge"}
                    if ( index == "grouped") {ind=TRUE} else (ind=FALSE)

                    df <- rbind(df, data.frame(iterations=e,
                                               group=n,
                                               accuracy=f,
                                               seed=ede,
                                               classification=response,
                                               cv=ncv,
                                               method= me,
                                               grouped=ind,
                                               lambda=bestlam,
                                               totalNgenes=dim(means)[1],
                                               regNgenes=len,
                                               trainingPercent=c(tr/ncol(means)*100)))
                }
            } else if ( regression == TRUE ) {
                ## Test set MSE only for regression-type analysis
                mean((lasso.labels - y[test])^2)
            } else {
                stop("Data must be fit as either a classification or regression model")
            }
        }
        ## plot the lot of iterations
        pred <- prediction(ps, dl)
        perf <- performance(pred, 'tpr', 'fpr')
        par(new=TRUE)
        plot(perf, lty=3, col=couleurs[i],
             xlab="Specificity (1-False positive rate)",
             ylab="Sensitivity (True positive rate)")
        par(new=TRUE)
        plot(perf,col=couleurs[i],lty=1,lwd=2,avg="vertical",spread.estimate="stderror",add=TRUE)
    }
legend("bottomright", levels(y)[nl], lty=1, lwd=5, col = couleurs[nl])
dev.off()
}


## extract regularization metrics for all iterations
write.table(dm, paste0("logSummary.lambda.iterations",iterations,".",response,".probabilities.txt"), sep="\t", quote=F)
write.table(df, paste0("logSummary.lambda.iterations",iterations,".",response,".accuracies.txt"), sep="\t", quote=F)




#########################
## Get best parameters ##
#########################
# iterate across sample labels
set.seed(ed)
nl <- c(1:nlevels(y))[levels(y)!="CTRL"]
results=NULL
for ( z in nl ) {
    ## get lambda and iteration seed based on maximum accuracy
    patient_labels<- levels(y)[z]
    max_accuracy <- df %>%
        filter(group == patient_labels) %>%
        filter(accuracy == max(accuracy)) %>%
        filter(regNgenes == max(regNgenes))

    results <- rbind(results, max_accuracy)
}

 
# get non-zero genes for best parameters
# each sample group has a different lambda
# each lambda was generated from a different seed
# each lambda generates a different yet somewhat redundant gene subset
lasso.selected.genes = NULL
for (l in 1:nrow(results)) {
    ed <- results$seed[[l]]
    bestlam <- results$lambda[[l]]
    set.seed(results$seed[[l]])
    patient_labels <- results$group[[l]]

    # get lasso coefficients
    training <- sample(1:nrow(adj.x), nrow(adj.x)/1.25)
    lasso.trained <- glmnet(adj.x[training,],
                            y[training],
                            alpha=setalpha,
                            lambda=bestlam,
                            family = response,
                            standardize=F,
                            type.multinomial=index)

    # get gene coefficients at selected lambda
    lasso.coef <- predict(lasso.trained, s=bestlam, type = "coefficients")
    selected.genes <- lasso.coef[[1]]@i[ lasso.coef[[1]]@i >= 1]
    original.genes <- colnames(adj.x)
    selected.final <- original.genes[selected.genes]
    len <- length(selected.genes)

    ## extract expression of regularized genes
    if ( length(selected.final) == length(selected.genes) ) {
        lasso.select <- adj.x[, selected.final ]
        dim(lasso.select)
        write.table(t(lasso.select), paste0("expressions.",patient_labels,
                                         ".iterations",iterations,".cv",ncv,
                                         ".lambda_",sprintf("%.5f", bestlam),
                                         "_",response,".regularization",setalpha,
                                         ".",index,".genes",len,".seed",ed,".txt"),
                    quote=F,sep="\t")
    } else {
        stop("Number of selected genes do not match the original dataset")
    }

    # get all genes selected from each grouping
    lasso.selected.genes <- rbind(lasso.selected.genes, t(lasso.select))
}

# plot lambda iterations
#par(mfrow = c(3,4))
#pdf(paste0("grid.lambda.",patient_labels,response,
#           ".regularization",setalpha,".",index,
#           ".genes",len,".seed",ed,".pdf"))
#plot(lasso.trained, xvar="lambda", label=T)
## fraction deviance explained =R2
#plot(lasso.trained, xvar="dev", label=T)
#plot(cv.out)
#dev.off()


############################
## Validating classifiers ##
############################
# combine TRUE sample grouping (observed) with predicted
# lasso.select contains the genes that were selected for their
# high probability to assign samples to their correct group
lasso.selected.genes.nodup <- unique(lasso.selected.genes[, ])
dat <- data.frame(y=y, t(lasso.selected.genes.nodup))
dim(dat)
set.seed(ed)


# machine learning models used
model_types <- c("svmLinear", "svmPoly", "svmRadialSigma", "svmLinear3",
                 "lda2", "bagFDA", "fda", "pda", "loclda", "bagFDAGCV",
                 "C5.0", "LogitBoost", "regLogistic",
                 "kernelpls", "multinom",
                 "nnet", "monmlp", "dnn",
                 "rf", "RRF",
                 "kknn", "naive_bayes", "gbm")

# number of parameters per model
# the deep network used in this step is an automated model
# hence the low number of parameters to adjust
parameter_counts <- c(1,3,2,2,
                      1,2,2,1,1,1,
                      3,1,3,
                      1,1,
                      2,2,5,
                      1,3,
                      3,3,4)

######## STEP I ##########
## performance metrics for best model without hyperparameter optimization
# create an index for the iteration holder below
icc <- function(){ i=0; function(){ i <<- i + 1;  i }}
modet=ite=NULL
ite=icc()
models=icc()

for ( iterations in c(1:2) ) {
    ## as many iterations will be executed for each model
    ## iterations are done in addition to 25 resampling for each model
    ie=ite()

    for ( m in 1:length(model_types) ) {
        
        mods=model_types[[m]]
        param=parameter_counts[[m]]

        ## models are trained in succession
        ## output is saved
        start <- format(Sys.time(), "%X")
        cat("\nIteration", ie, "on model", mods, "started at", start)
        modnam=models()
        model.name <- paste0(mods,"|",ie,"|",param)
        
        trainCtrl <- trainControl(method="cv",
                                  number=10)

        trained.model <- train( y ~ .,
                               data = dat,
                               trControl= trainCtrl,
                               method = mods )

        # aggregate all performance metrics
        if ( modnam == 1 ){
            performance_summary <- list(trained.model)
            names(performance_summary)[modnam] <- model.name
        } else if ( modnam > 1 ) {
            performance_summary <- c(performance_summary, list(trained.model))
            names(performance_summary)[modnam] <- model.name        
        }
        end <- format(Sys.time(), "%X")
        cat(". Execution was successful at", end)
    }

}

## summary 1
## performance metrics between models with repeated iterations
sink(paste0("performance1.multianalysis.seed",ed))
performance_summary %>% resamples %>% summary
sink()

## save output
save(list=ls(pattern="perf*summary"),file="performanceSummaryNoTune.Rdata")


######## STEP II ##########
## Classification across models with hyperparameter optimization
## with hyperparameter tuning
systems.metrics <- NULL

if ( classification == TRUE & grouped == TRUE ) {
    for ( m in 1:length(model_types) ) {
        
        mods=model_types[[m]]
        param=parameter_counts[[m]]

        ## start logging
        start <- format(Sys.time(), "%X")
        cat("\nTraining on model", mods, "started at", start)

        ## models are trained in succession
        ## predicted output is saved
        model.name <- paste0(mods,"|",param)

        model.metrics <- modelTune.clas(dat,training,method=mods,folds=10,r=5,tune=30, grid=TRUE)

        ## aggregate all performance metrics
        ## only for predicted features
        ## contains results and tuned/selected parameters per model
        if ( m == 1 ){
            performance_summaryFull <- list(model.metrics$bestModel)
            names(performance_summaryFull)[m] <- model.name
        } else if ( m > 1 ) {
            performance_summaryFull <- c(performance_summaryFull, list(model.metrics$bestModel))
            names(performance_summaryFull)[m] <- model.name        
        }

        ## get iteration speed
        ## get sensitivity, specificity, precision scores ...
        systems.metrics <- rbind(systems.metrics,
                                 data.frame(model=mods,
                                            durationSeconds=model.metrics$timeLapsed[[3]],
                                            model.metrics$ConfusionMatrix$byClass,
                                            accuracy=model.metrics$ConfusionMatrix$overall[[1]],
                                            kappa=model.metrics$ConfusionMatrix$overall[[2]],
                                            accuracyPval=model.metrics$ConfusionMatrix$overall[[6]]))

        ## end logging
        end <- format(Sys.time(), "%X")
        cat(". Execution was successful at", end)
    }

} else if ( classification == TRUE & binomial == TRUE ) {

    output_summary <- modelTune.clas(dat,training,method="rf",folds=10,r=10,tune=30, grid=FALSE)

} else if ( regression == TRUE & bionomial == TRUE ) {

    output_summary <- modelTune.reg(dat,training,method="nnet",folds=10,r=10,tune=30,ctl)

}

## summary 2
## performance metrics for best model while hyperparameters tuning
sink(paste0("performance2.hyperTuning.seed",ed))
performance_summaryFull %>% resamples %>% summary
sink()

## summary 3
## performance metrics for best classifier
write.table(systems.metrics,
            paste0("performance3.full.hyperTuning.seed",ed),
            sep = "\t", quote = FALSE)

## save classification output
save(list=ls(pattern="systems.metrics"),file="systemsMetrics.Rdata")
save(list=ls(pattern="perf*Full"),file="performanceSummaryFull.Rdata")




##### debugging #####
#model.reg(dat,training,method="rf",folds=10,r=5,tune=10)
#train(y~., data = dat, method = "naive_bayes")


###### experimental ######
#compare_models(nnet0,rf1)


## print full view of variables, data frames, matrices, funcitons...
lsos()



##############################
# Ensemble Learning (bagging)
##############################
require(foreach)
require(doSNOW)
cl <- makeCluster(3)
registerDoSNOW(cl)
#ctl=expand.grid(.size=17, .decay=0)	## nnet
ctl=expand.grid(.mtry=0,.coefReg=0.89,.coefImp=0.5556)## Regularized Random forest (RRF)
#ctl=expand.grid(.C=8, .sigma=0.5556)	## svmRadial
#ctl=expand.grid(.ncomp=1)	# PCR
#ctl=expand.grid(.C=1)	## svmLinear
#ctl=expand.grid(.mstop=500,.prune="no")	## glmboost
#ctl=expand.grid(.lambda=0.007499)	## Ridge
dat <- data.frame(y=y, lasso.select); dim(dat)
set.seed(1445612321)
RFF_bag500<- baggingTune(dat[training,],dat[-training,],m=1.1,ite=500,methods="RRF",tune=10,gridZ=ctl)## For tuning the hyper-parameters
test100 <- bagging(dat[training,],dat[-training,],m=1.1,ite=100,methods="ridge",tune=10)## for testing
bagging.clas(dat[training,],dat[-training,],m=1.1,ite=100,methods="nnet",tune=10)## For Classification
## END RUN
stopCluster(cl)		## close cluster only after finishing w all modelse

## COMPUTE Ensemble RMSE
## (1) load ensembleMethods.Rdata
## (2) re-run y[test]
ensemble.mean(svmLinear_bag2K,ridge_bag200)

require(lattice)
ensRMSE <- read.table("clipboard",sep="\t",header=T); ensRMSE
xyplot(Method1+Method2+Method3~Model2|Model1,data=ensRMSE,type=c("p","a","g"),auto.key = list(columns=3,points=F,lines=T,title="Ensemble methods averaging"), ylab="Root-mean-square error (RMSE)",xlab="Tested learners",pch=21,cex=1)
## Bagging

##############################
# Unsupervised PCR and Clustering
##############################
## START PCR
pr.out <- prcomp(lasso.select, scale=T)
y <- c(rep("E",3), rep("T",3), rep("VC",3),rep("PC",3),rep("JC",3),rep("VT",3),rep("PT",3), rep("JT",1))
y <- c(rep("E",3), rep("L",9), rep("PL",3), rep("L",6),rep("PL",1))
y <- c(rep("E",3), rep("L",9), rep("PL",3), rep("LT",6),rep("PLT",1))
y <- c(rep("Healthy", 15), rep("Deficient",7))
y <- c(rep("None", 6), rep("Cocktail", 9), rep("Tiso",7))
y <- c(rep("L",9), rep("PL",6), rep("L",3),rep("PL",4))
par(mfrow = c(3,2))
plot(pr.out$x[,1:2],col=Cols(y),pch=19)
plot(pr.out$x[,c(1,3)],col=Cols(y),pch=19)
## 2D scatterplots
require(scatterplot3d)
scatterplot3d(pr.out$x[,1:3], pch=21, color=Cols(y), type="h", lty.hplot = "dashed", angle = 55, scale.y = .7)
## Plot 1 in 3D
require(rgl)
plot3d(pr.out$x[,1:3], size=15, col=Cols(y), type = "p", box=T, axes=T,top=F)
par3d(zoom = 1.1) # larger values make the image smaller
par3d(windowRect = c(50, 50, 500, 500)) # make the window large
text3d(pr.out$x[,1:3], text=y, font=1, cex=0.8, adj=c(1,1.5))
setwd("~/Downloads");rgl.postscript("PCR.eps",fmt="eps")
## plot 2 in 3D
## START recording
M <- par3d("userMatrix") # save your settings to pass to the movie
movie3d(par3dinterp(userMatrix=list(M,rotate3d(M, pi, 1, 0, 0), rotate3d(M, pi, 0, 1, 0) ) ), duration = 5, fps = 50, movie = "MyMovie", dir = ("~/Downloads/PCA3D"))
## extract multiple frames for GIF implementation
#movie3d(spin3d(axis = c(0, 0, 1), rpm = 4), duration=15, movie="TestMovie", type="gif", dir = ("~/Downloads/PCA3D"))
#rgl.snapshot('PCA3D.png', fmt="png")
## extract 2
## END
## 3D scatterplot
summary(pr.out)
pve <- 100*pr.out$sdev^2/sum(pr.out$sdev^2)
plot(pr.out)
plot(pve, type="o", ylab = "PVE",xlab = "Principal Component", col="blue")
plot(cumsum(pve), type="o",ylab = "Cumultive PVE", xlab = "Principal Component", col="brown3")
## (Unsupervised) Principal Componenent analysis (reduced dimension) (paper3)


##############################
# Tryouts
##############################



require(qmap)
qmap("Magdalen Islands", zoom=14)
## specific place map

sd.data <- scale(lasso.select)
par(mfrow = c(1,1))
data.dist <- dist(sd.data)
plot(hclust(data.dist), labels=y, main="Complete Linkage")
plot(hclust(data.dist, method="single"), labels=y, main="Complete Linkage")
plot(hclust(data.dist, method="average"), labels=y, main="Average Linkage")
abline(h=18, col="red")
## Hierarchicla clustering
hc.out <- hclust(data.dist)
hc.clusters <- cutree(hc.out, 4)
table(hc.clusters, y)
hc.out
## cut the dendrogram at the height that yield a particular number
set.seed(1)
km.out <- kmeans(sd.data, 4, nstart = 20)
km.clusters <- km.out$cluster
table(km.clusters, hc.clusters)
## K means clustering
hc.out <- hclust(dist(pr.out$x[,1:6]))
plot(hc.out, labels=y, main="Hier. Clust., on First Six Score Vectors")
table(cutree(hc.out, 4), y)
## perform clustering on the first 6 rpincipal components

set.seed(2)
train <- sample(1:nrow(x), 2*nrow(x)/2.5)
test <- -train
dat <- data.frame(x[train,1:5000], y=as.factor(y[train]))
require(e1071)
out <- svm(y~., data=dat, kernel = "radial",cost=10, gamma=1)
summary(out)
table(out$fitted, dat$y)
## Training observations
dat.te <- data.frame(x[test,1:5000], y=as.factor(y[test]))
pred.te <- predict(out, newdata=dat.te)
table(pred.te, dat.te$y)
cat("\nMisclassification Error rate",(4/22)*100,"\n")
## support vector regression



####### Save ########
lsos(pat="locus.select|mRMR|grid|bag")
save(list=ls(pattern="mRMR|grid|bag"),file="EnsembleMethods.Rdata")	## save
load("EnsembleMethods.Rdata", .GlobalEnv)
