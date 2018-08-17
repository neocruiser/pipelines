pkgs <- c('RColorBrewer', 'pvclust', 'gplots',
          'dplyr', 'glmnet', 'caret', 'foreach',
          'doSNOW', 'lattice', 'ROCR', 'earth', 'vegan',
          'reshape2', 'ggplot2', 'tidyr', 'plyr',
          'plot3D', 'ggrepel', 'ggdendro')
lapply(pkgs, require, character.only = TRUE)

## logging
source("./scripts/lsos.R")


##############
## SWITCHES ##
##############
### DEFINE DATA FITTING MODEL
classification = TRUE
regression = FALSE

## define feature structure
grouped = TRUE
binomial = FALSE

## apply normalization methods within samples
## applying a regularization transformation
## increases the number of selected genens
standardization = TRUE
lasso.std = FALSE

## how many iterations during feature selection (max 50)
## how many iterations during baseline calculation of model testing (max 10)
## how many CV folds during model validation (max 10)
## use grid tuning or not for model optimization
## length of the training data set (percent, max 85)
lambda.epochs = 5
baseline.epochs = 5
validate.folds = 10
validate.tune = FALSE
train.size = 80

## choose between methods of normalization and correlation
## variable will be extracted from ids2modules
## the latter is generated from gene networks
## which correlation, power, normalization methods used from network analysis
## choose from ids2modules colnames
clustering.strategy = 8

## choose contrasts
## choose summary file with all FDR adjusted pvalues
## file generated from expression analysis
grouping = c("systemicRelapse", "systemicRelapseNodes", "systemicRelapseCOOprediction")
pvals <- read.table("summary/summary.lmfit.all.txt", header = TRUE, fill = TRUE)
## names of each sample files
ids <- read.table("summary/sampleIDs")
 
## ids2modules is a summary file generated from the weighted nets script
## ids2description is an annotation file of the selected genes from the networks
ids2modules <- read.table("./ids2modules.summary.txt", header = TRUE)
colnames(ids2modules)
ids2description <- read.table("./ids2description.summary.txt", header = FALSE, fill = TRUE)
colnames(ids2description) <- c("genes", "chromosome", "ensembl", "symbol", "function", "site", "symbol2")



##########################
## Define new functions ##
##########################
## Stochastic mini-batch sampling using heuristic selection
miniBatch.balancedSampling <- function (meta.selected, y, adj.x, batch = 65, miniBch = 85, repl = FALSE) {
    ## when imbalanced samples are found
    ## select an equal distribution of samples across classes
    ## batch sampling from the original dataset
    ## mini-batch subsampling from the balanced batch

    ## get minimum length of smallest imabalanced class
    ## select mini batch size
    min.class <- min(as.data.frame(table(meta.selected$Groups))$Freq)
    mini.batch <- round( (batch / 100) * nrow(adj.x) / nlevels(y), 0)
    raw.index <- data.frame(index = 1:nrow(adj.x), samples = row.names(adj.x))

    ## create first random sub sampled data set
    if ( mini.batch <= min.class ) {
        sub.training=NULL
        for (i in 1:nlevels(y) ) {
            selected.group <- meta.selected %>%
                select(SAMPLE_ID, Groups) %>%
                filter(Groups == levels(y)[i] )

            mini.class.samples <- sample(selected.group$SAMPLE_ID, mini.batch, replace = repl)
            final <- raw.index[ raw.index$samples %in% mini.class.samples, 1]
            sub.training <- c(sub.training, final)
        }
    } else {
        print("Error! Batch size is too large. Set a lower number to generate a smaller stochastic mini-batch than the minimum length of the smallest imbalanced class.")
    }

    ## create a balanced random training set
    tr <- length(sub.training) * c( miniBch / 100 )
    training <- sample(sub.training, tr)
    return(training)
}



## Multi model classification
modelTune.clas <- function(dat, train, method, folds=10, grid=TRUE, confusion.metrics=NULL){
    ## requires caret
    ## GRID search HYPERPARAMETERS tuning
    ## Cross validation for parameter tuning
    ## minimum 10 epochs repeated 5 times
    ## over 20 different flexible and less flexible machine learning algorithms
    ## output are logged under different formats for different performance summarizations
    trainCtrl <- trainControl(method="repeatedcv",
                              number=folds,
                              repeats=c(folds/2),
                              summaryFunction=defaultSummary)
    
    ## Choosing the right Hyper-parameters. GRID ANALYSIS
    if ( grid == TRUE ) {

        ## Tune hyper-parameters
        if ( method == "svmLinear" ) {
            ## support vector machines with linear kernel
            grid_models <- expand.grid(.C=seq(.001,25, length=100))
        } else if ( method == "svmPoly") {
            ## svmRadial
            grid_models <- expand.grid(.degree=seq(0,10,.5),
                                       .scale=10^seq(-1,-3,length=10),
                                       .C=seq(.1,2, length=20))
        } else if ( method == "svmRadialSigma") {
            ## svm with radial basis function kernel
            grid_models <- expand.grid(.C=seq(.1,2, length=20), .sigma=10^seq(-1,-5,length=60))
        } else if ( method == "svmLinear3") {
            ## l2 regularized support vector machine (dual) with linear kernel
            grid_models <- expand.grid(.cost=10^seq(-1,2,30), .Loss=seq(0:7))
        } else if ( method == "lda2" ) {
            ## Linear discriminant analysis
            grid_models <- expand.grid(.dimen=seq(10^-2,40,length=70))
        } else if ( method == "bagFDA" || method == "fda" ) {
            ## bagged flexible discriminant analysis
            grid_models <- expand.grid(.degree=seq(.1,2,length=20), .nprune=seq(1,60,length=30))
        } else if ( method == "pda" ) {
            ## penalized discriminant analysis
            grid_models <- expand.grid(.lambda=10^seq(-1,-5,length=50))
        } else if ( method == "loclda" ) {
            ## localized linear discriminant analysis
            grid_models <- expand.grid(.k=seq(1,350,length=10))
        } else if ( method == "bagFDAGCV" ) {
            ## bagged FDA using gCV pruning
            grid_models <- expand.grid(.degree=seq(10^-3,5,length=100))
        } else if ( method == "C5.0" ) {
            ## c5 decision trees
            grid_models <- expand.grid(.trials=seq(1:200), .model=c("rules","tree"), .winnow=c(TRUE,FALSE))
        } else if ( method == "LogitBoost" ) {
            ## Boosted logistic regression
            grid_models <- expand.grid(.nIter=seq(1,100,length=100))
        } else if ( method == "regLogistic" ) {
            ## Regularized logistic regression
            ## the loss parameter triggers L1 and L2 loss regularizations
            ## its suggested the difference between primal and dual (through loss)
            ## has effects on speed of execution and little on accuracy performance
            grid_models <- expand.grid(.cost=seq(.1, 2, length=20),
                                       .loss=seq(0:7),
                                       .epsilon=10^seq(-1,-2,length=10))
        } else if ( method == "kernelpls" ) {
            ## partial least squares
            grid_models <- expand.grid(.ncomp=seq(1,20,length=40))
        } else if ( method == "multinom" ) {
            ## penalized multinomial regression
            grid_models <- expand.grid(.decay=10^seq(-1,-5,length=200))
        } else if ( method == "rf" ) {
            ## random forest
            grid_models <- expand.grid(.mtry=seq(1,40,length=100))
        } else if ( method == "RRF" ) {
            ## regularized random forest
            grid_models <- expand.grid(.mtry=seq(1:10),.coefReg=10^seq(-1,-5,length=10),.coefImp=10^seq(-1,-5,length=10))
        } else if ( method == "kknn" ) {
            ## weighted k nearest neighbors
            grid_models <- expand.grid(.kmax=seq(1,15,length=20),.distance=seq(1,5,length=10),
                                       .kernel=c("optimal","rank","gaussian",
                                                 "inv","cos","rectangular",
                                                 "triangular","biweight",
                                                 "triweight"))
        } else if ( method == "naive_bayes" ) {
            ## naive bayes
            grid_models <- expand.grid(.laplace=seq(10^-2,2,length=10),.usekernel=c(FALSE,TRUE),
                                       .adjust=c(10^-2,2,length=10))
        } else if ( method == "gbm" ) {
            ## stochastic gradient boosting
            grid_models <- expand.grid(.n.trees=seq(5,300,length=10),.interaction.depth=seq(1:3),
                                       .shrinkage=10^seq(-2,-4,length=5),.n.minobsinnode=seq(5,10,length=3))
        } else if ( method == "nnet" ) {
            ## neural networks
            grid_models <- expand.grid(.size=seq(1:5), .decay=10^seq(-1,-2,length=7))
        } else if ( method == "pcaNNet" ) {
            ## neural networks with inclusive feature extraction
            grid_models <- expand.grid(.size=seq(1:5), .decay=10^seq(-1,-2,length=7))
        } else if ( method == "monmlp" ) {
            ## monotone multi-layer perceptron neural network
            grid_models <- expand.grid(.hidden1=seq(1:5),.n.ensemble=c(5))
        } else if ( method == "mlpSGD" ) {
            ## multilayer perceptron network by stochastic gradient descent
            grid_models <- expand.grid(.size=seq(1:3),.l2reg=10^seq(-3,-5,length=2),.lambda=0,
                                       .learn_rate=10^seq(-1,-7,length=3),
                                       .gamma=10^seq(0,-1,length=3),
                                       .momentum=10^seq(0,-1,length=3),
                                       .minibatchsz=seq(1,120,length=15),    
                                       .repeats=2)
        } else if ( method == "mxnet" ) {
            ## deep neural network with GPU computing
            ## relu (rectified linear units) faster than sigmoid function
            ## because relu converge faster (faster learning)
            ## also, a reduced likelihood of vanishing gradient (weights and biases)
            ## and adds sparcity resulting in less dense network representation
            ## length 6 was the maximum threshold ie functional in 60 hours
            grid_models <- expand.grid(.layer1=seq(1:15),.layer2=seq(1:3),.layer3=seq(1:2),
                                       .learning.rate=10^seq(-1,-7,length=4),
                                       .momentum=10^seq(0,-1,length=4),
                                       .dropout=10^seq(0,-5,length=4),
                                       .activation=c("relu"))
        } else if ( method == "mxnetAdam" ) {
            ## deep neural network
            num.round=1
            grid_models <- expand.grid(.layer1=seq(1:15),.layer2=seq(1:3),.layer3=seq(1:2),
                                       .dropout=10^seq(-1,-5,length=3),
                                       .beta1=seq(0,1,length=3),
                                       .beta2=seq(0,1,length=3),                                       
                                       .learningrate=10^seq(-2,-7,length=3),
                                       .activation=c("relu"))
        } else if ( method == "dnn" ) {
            ## stacked autoencoder deep neural network
            grid_models <- expand.grid(.layer1=seq(1:10),.layer2=seq(1:3),.layer3=seq(1:2),
                                       .hidden_dropout=10^seq(-1,-7,length=10),
                                       .visible_dropout=10^seq(-1,-7,length=10))
        }

        ## train the model
        lapsed <- system.time(modelTrain <- train(y~., data=dat[train,],
                                                  method=method,
                                                  trControl= trainCtrl,
                                                  preProc=c("center","scale"),
                                                  ## weights = model_weights,
                                                  tuneGrid=grid_models,
                                                  tuneLength=c(folds*3)))

    } else if ( grid == FALSE ) {
        ## do not tune the hyperparameters
        ## just use the default options for each model
        lapsed <- system.time(modelTrain <- train(y~., data=dat[train,],
                                                  method=method,
                                                  trControl= trainCtrl,
                                                  preProc=c("center","scale")))
    }

    ## get performance scores from out-of-bag observations
    ## after cross validation during hyperparameter tuning
    for ( iterations in c(1: c(folds * 20)) ) {
        ## iterative prediction test for model performance
        edd <- floor(abs(rnorm(1) * 10000000))
        set.seed(edd)
        subtrain <- sample( 1:nrow(dat), c(nrow(dat) * c(train.size/100)) )

        cat("\n  >>> Sub-iteration", iterations, "on model", method, "started")
        ## predict over best hyperparameters
        end <- system.time(Predd <- predict(modelTrain, newdata=dat[-subtrain,], type="raw"))
        conf.m <- confusionMatrix(data=Predd, dat[-subtrain,1])
        ## aggregate all prediction metrics
        confusion.metrics <- rbind(confusion.metrics,
                                   data.frame(iteration=iterations,
                                              model=method,
                                              seed=edd,
                                              conf.m$byClass,
                                              accuracy=conf.m$overall[[1]],
                                              accLow=conf.m$overall[[3]],
                                              accHigh=conf.m$overall[[4]],
                                              kappa=conf.m$overall[[2]],
                                              accPval=conf.m$overall[[6]]))
        cat(" &", method, "validated successfully in (ms)", c(end[[1]]*1000), "@accuracy", conf.m$overall[[1]])
    }

    ## get prediction accuracy
    ## confusion matrix for classification
    Predd <- predict(modelTrain, newdata=dat[-train,], type="raw")
    conf.m <- confusionMatrix(data=Predd, dat[-train,1])

    ## Compile everything for later accessibility
    ## time, best model based on training performance metrics
    ## tuning accuracy per parameter
    output <- list(timeLapsed=lapsed,
                   bestModel=modelTrain,
                   Results=modelTrain$results,
                   Hyperparameters=modelTrain$bestTune,
                   ConfusionMatrix=conf.m,
                   metrics=confusion.metrics)
    return(output)
}


## get probabilityy scores for support vector machines (used as best predictor)
svmPoly.prob <- function(dat, train, folds = 10, risk.prob = NULL){
    ## requires previous functions
    ## requires hyperparameter tuning on best model
    ## execution will run on 100 epochs
    hp <- parameters_summaryFull[["svmPoly"]]
    trainCtrl <- trainControl(method="repeatedcv",number=folds,repeats=c(folds/2),classProbs=T)
    grid_models <- expand.grid(.degree=hp[[1]], .scale=hp[[2]], .C=hp[[3]])
    cat("\nModel training initiated [ok]\n")
    modelTrain <- train(y~.,
                        data=dat[train,],
                        method="svmPoly",
                        trControl= trainCtrl,
                        preProc=c("center","scale"),
                        tuneGrid=grid_models,
                        tuneLength=c(folds*3))
    cat("Model training finished [ok]\n")    
    for ( iterations in c(1: c(folds * 20)) ) {
        cat("  --> iteration", iterations, "[ok]\n")
        ## iterative prediction test for model performance
        edx <- floor(abs(rnorm(1) * 10000000))
        set.seed(edx)
        subtrain <- sample( 1:nrow(dat), c(nrow(dat) * c(train.size/100)) )
        Predd <- predict(modelTrain, newdata=dat[-subtrain,], type="prob")
        risk.prob <- rbind(risk.prob, data.frame(epochs = iterations, original = y[-subtrain], Predd))
    }
    cat("Model validation finished [ok]")    
    return(risk.prob)
}


##########################
## Load expression data ##
##########################
## optimized for t-statistics microarray expressions
## rows are genes & col are samples
## remove controls (imbalanced samples)
## append gene names instead of transcript numbers
cat("\n\nNormalized expression scores: Samples are columns and genes are rows\n")
means <- read.table("expressions", sep="\t", header=T, row.names=1) %>%
    select(-matches("CNR800.T1")) %>%
    select(-matches("CNR900.T1")) 
pre.annot <- dim(means)[2]

gene.names <-  ids2description[ ids2description$genes %in% rownames(means), c(1:2,4:5)]
gene.names$chromosome <- as.factor(gsub("_.*","",gene.names$chromosome))
write.table(gene.names, "log.gene.names_allnetworks.txt", sep = '\t', row.names = FALSE, quote = FALSE)

idx <- gene.names[, -4] %>%
    mutate(id = paste0(symbol,"-",gsub(".hg.1","",genes)))

means$genes <- rownames(means)
new.means <- full_join(means, idx, by = "genes")
rownames(new.means) <- new.means$id
means <- new.means %>%
    select(-genes, -id, -symbol, -chromosome)
post.annot <- dim(means)[2]

if ( pre.annot != post.annot ) { stop("Can't continue, newly indexed data differ from the original") }


## normalization of gene expression
cat("\n\nStandardized transformed scores: Genes are columns and samples are rows\n")
if ( standardization == TRUE ) {
    adj.x <- t(decostand(means, "standardize"))
    print("Data have been transformed as requested for better interpretability and overall fitting under flexible models")
} else {
    adj.x <- t(means)
    print("Data were not transformed and kept as is.")
}
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
    filter(Groups != "CTRL") %>%
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

## reorder sample names
metadata <- arrange(metadata, factor(SAMPLE_ID, levels = ids$V1))

meta.selected <- metadata %>%
    mutate(Contrast1 = as.factor(paste0(Groups, ".", Prediction))) %>%
    mutate(Contrast2 = as.factor(paste0(Groups, ".", Nodes)))


## choose samples strutcture from metadata
sample.classes = "Groups"
y <- meta.selected$Groups

## designate a multivariate analysis
if ( nlevels(y) >= 2 ) {
    response="multinomial"
} else if ( nlevels(y) == 2 ) {
    ## experimental
    response="binomial"
} else {
    ## experimental
    response="mgaussian"
}

## prepare testing dataset
# set seed for reproducibility
ed <- floor(abs(rnorm(1) * 10000000))
set.seed(ed)

## Split the dataset
training <- sample(1:nrow(adj.x), c(nrow(adj.x) * c(train.size/100)) )
tr <- length(training)

##  proportions (in percentages) of imbalanced samples
print(round((table(y[training])/tr) * 100),2)
sink("1-expression.data.loaded.ok")
sink()

########################
## FEATURE EXTRACTION ##
########################

     ##############################################################################
     ## Result: Lambda for Lasso feature selection with highest accuracy outcome ##
     ## Reduce collinearity or excessive correlation among genes                 ##
     ## improve identification of optimal set of variables                       ##
     ##  feature extraction                                                      ##
     ##  |-- randomize seed                                                      ##
     ##  |   |-- iterate multiple seeds                                          ##
     ##  |   |-- mini-batch sampling                                             ##
     ##  |   |-- imbalanced sample weighting (optional)                          ##
     ##  |                                                                       ##
     ##  |-- parameter lambda tuning                                             ##
     ##  |   |-- grid search                                                     ##
     ##  |   |-- nested cross validation                                         ##
     ##  |                                                                       ##
     ##  |-- fit linear model                                                    ##
     ##  |   |-- first on training set                                           ##
     ##  |   |-- second on testing set                                           ##
     ##  |   |-- iterate multiple model fitting over multiple epochs             ##
     ##  |                                                                       ##
     ##  |-- get prediction scores for each iteration                            ##
     ##  |-- plot accuracy scores                                                ##
     ##  |   |-- iterate multiple accuracy tests                                 ##
     ##  |   |-- ROC curves for Cross validation test                            ##
     ##  |   |-- ROC curves for prediction/validation set                        ##
     ##  |                                                                       ##
     ##  |-- create summary of all iterations                                    ##
     ##  |                                                                       ##
     ##  |-- get the best lambda                                                 ##
     ##  |   |-- best accuracy has the best lambda                               ##
     ##  |   |-- use this lambda for subsequent analyses                         ##
     ##  |   |-- plot final accuracy at this lambda cutoff                       ##
     ##############################################################################

## WARNING: sometimes LOGNET throws an error of 0 or 1 observations
## this is due to the unbalanced nature of cross validation
## SOLUTION: the while condition will repeat the test until success
success=FALSE

while (success == FALSE) {
    pdf(paste0("cvROC.shrinking.epochs",lambda.epochs,".",response,".pdf"))
    couleurs <- brewer.pal(nlevels(y), name = 'Dark2')
    dm=NULL
    df=NULL
    gpd=NULL

    ##Remove levels without observations
    ##obsLevels <- levels(droplevels(allObs))
    
    ## The control class is throwing an error
    ## prediction cannot be done on small sample size
    ## remove CTRL class
    if ( sample.classes == "Contrast1" ) {
        nl <- c(1:nlevels(y))[ levels(y) != c("CTRL.ABC", "CTRL.GCB", "CTRL.NA", "CTRL.U")]
    } else if ( sample.classes == "Contrast2" ) {
        nl <- c(1:nlevels(y))[ levels(y) != c("CTRL.EN", "CTRL.LN")]
    } else if ( sample.classes == "Groups" ) {
        nl <- c(1:nlevels(y))[ levels(y) != "CTRL"]
    }

    if (!require(ROCR)) {    stop("Can't continue can't load ROCR")   }

    for (nid in nl) {
        gc()
        
        for (e in 1:lambda.epochs) {
            ## set seed for reproducibility
            ede <- floor(abs(rnorm(1) * 10000000))
            set.seed(ede)

            ## Split the dataset into 65% training data
            training <- miniBatch.balancedSampling(meta.selected, y, adj.x, batch = 70, miniBch = train.size, rep = FALSE)
            tr <- length(training)
            print(round((table(y[training])/tr) * 100),2)

            ## fit a generalized linear model via penalized maximum likelihood
            ## if alpha 1 then lasso L1 penality and discard genes
            ## if alpha 0 then ridge regression then L2 and rank genes
            if ( grouped == TRUE ){
                index="grouped"
                ncv=10
                setalpha=1
                nfold <- 10

                ## fitting a symmetric multinomial model,
                grid <- 10^seq(5, -5, length=100)
                lasso.trained <- glmnet(adj.x[training,],
                                        y[training],
                                        alpha = setalpha,
                                        lambda = grid,
                                        family = response,
                                        standardize = lasso.std,
                                        type.multinomial=index)

                ## Cross validation for hyperparameter tuning.
                ## Optimization of model selection to avoid overfitting
                ## also, silence errors due to unbalanced observations
                cv.out <- try(cv.glmnet(adj.x[training,],
                                        y[training],
                                        alpha=setalpha,
                                        family = response,
##                                        foldid = foldid.custom,
##                                        weights = weights.class,
                                        standardize = lasso.std,
                                        nfolds = ncv,
                                        type.multinomial=index),
                              silent = TRUE)
            } else {
                stop("Alternative feature indexes are not yet implemented")                
            }

            ## make sure CV did not encounter unbalanced observations
            ## if unbalanced obs exist, the whole iteration will repeat
            if (class(cv.out) != "cv.glmnet") {
                success=FALSE
                ## initiate iteration
                ## not working, must be set early
                ##if ( e > 1 ) { e = e - 1 } 
                break
            } else {
                success=TRUE
            }

            ## get best lambda
            ## lowest probability to overfit
            bestlam <- cv.out$lambda.min

            ## get probabilities
            lasso.response <- predict(lasso.trained, s=bestlam, newx=adj.x[-training,], type="response")

            ## create list of iteration scores for each class
            probability.scores <- as.data.frame(lasso.response)[, nid]
            dummy.labels <- as.vector(model.matrix(~0 + y[-training])[, nid])

            ## create multi dimensional list across all iterations
            if ( e == 1 ){
                ps <- list(probability.scores)
                dl <- list(dummy.labels)
            } else {
                ps <- c(ps, list(probability.scores))
                dl <- c(dl, list(dummy.labels))
            }

            ## rename listing headers
            names(ps)[e]=paste0(levels(y)[nid],"-class",nid,".epochs",e)
            names(dl)[e]=paste0(levels(y)[nid],"-class",nid,".epochs",e)
            
            ## create summary of probabilities
            ## from a list of genes that have a mean probability to classify
            ## all cases into correct patient diagnosis
            dm <- rbind(dm, data.frame(epochs=e,
                                       class=levels(y)[nid],
                                       seed=ede,
                                       lambda=bestlam,
                                       probabilityScore=mean(ps[[e]]) ))

            ## extract probability of classification per sample
            ## during training
            ## to get general distribution of the whole classification
            gpd <- rbind(gpd, data.frame(probabilities = ps[[e]],
                                         classes = as.factor(y[-training]),
                                         genes = as.factor(levels(y)[nid])))


            ## compile accuracies into summary file
            ## get sample labels
            lasso.labels <- predict(lasso.trained, s=bestlam, newx=adj.x[-training,], type="class")

            ## get gene coefficients at selected best lambda
            lasso.coef <- predict(lasso.trained, s=bestlam, type = "coefficients")

            ## get non-zero genes
            selected.genes <- lasso.coef[[nid]]@i[ lasso.coef[[nid]]@i >= 1]
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
                    nas = names(rowSums(freq)[q])
                    f <- freq[nas, nas] / rowSums(freq)[[q]] * 100

                    if ( setalpha == 1 ) {me="lasso"} else {me="ridge/elastic"}
                    if ( index == "grouped") {ind=TRUE} else (ind=FALSE)

                    ## contains predicted probablities
                    df <- rbind(df, data.frame(epochs=e,
                                               patients=nas,
                                               featureGrouping=levels(y)[nid],
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

            cat(paste(">> Iteration", e, "finished for", levels(y)[nid] ))
        }
        ## plot the lot of iterations
        pred <- prediction(ps, dl)
        perf <- performance(pred, 'tpr', 'fpr')
        par(new=TRUE)
        plot(perf, lty=3, col=couleurs[nid],
             xlab="Specificity (1-False positive rate)",
             ylab="Sensitivity (True positive rate)")
        par(new=TRUE)
        plot(perf,col=couleurs[nid],lty=1,lwd=2,avg="vertical",spread.estimate="stderror",add=TRUE)
    }
    legend("bottomright", levels(y)[nl], lty=1, lwd=5, col = couleurs[nl])
    dev.off()
}


if ( exists('dm') && exists('df') ) {

    ## extract regularization metrics for all iterations
    write.table(dm, paste0("logSummary.lambda.epochs",lambda.epochs,".",response,".probabilities.txt"), sep="\t", quote=F)
    write.table(df, paste0("logSummary.lambda.epochs",lambda.epochs,".",response,".accuracies.txt"), sep="\t", quote=F)
    write.table(gpd, paste0("logSummary.lambda.epochs",lambda.epochs,".",response,".densities.txt"), sep="\t", quote=F)
    
    sink("2-feature.extraction.ok")
    sink()
}




## CHART 0
pdf(paste0("density.shrinking.epochs",lambda.epochs,".",response,".pdf"))
gpd %>%
##    filter(classes == c("CNS", "SYST", "NOREL")) %>%
    ggplot(aes(x = probabilities,
               fill = classes)) +
    geom_density(alpha = 0.3) +
##    geom_density(position = "stack") +    
    theme_minimal() +
    facet_wrap( ~ genes,
               ncol = 1,
               scales = "free") +
    theme(legend.position = "top") +
    scale_fill_brewer(palette = "Dark2") +
    labs(x = "Probability of selected genes to predict different patient classes",
         y = "Density (fitted probabilities with lambda = 1)")
dev.off()
try(dev.off(), silent = TRUE)



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
        filter(featureGrouping == patient_labels) %>%
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
    patient_labels <- results$featureGrouping[[l]]

    # get lasso coefficients
    index="grouped"
    training <- miniBatch.balancedSampling(meta.selected, y, adj.x, batch = 70, miniBch = train.size, rep = FALSE)
    lasso.trained <- glmnet(adj.x[training,],
                            y[training],
                            alpha=setalpha,
                            lambda=bestlam,
                            family = response,
                            standardize = lasso.std,
                            type.multinomial=index)

    # get gene coefficients at selected lambda
    lasso.coef <- predict(lasso.trained, s=bestlam, type = "coefficients")
    selected.genes <- lasso.coef[[1]]@i[ lasso.coef[[1]]@i >= 1]
    original.genes <- colnames(adj.x)
    selected.final <- original.genes[selected.genes]
    len <- length(selected.genes)

    ## extract expression of regularized genes
    if ( length(selected.final) == length(selected.genes) ) {
        lasso.select <- adj.x[ , selected.final ]
        print(dim(lasso.select))
        write.table(t(lasso.select), paste0("expressions.",patient_labels,
                                         ".epochs",lambda.epochs,".cv",ncv,
                                         ".lambda_",sprintf("%.5f", bestlam),
                                         "_",response,".regularization",setalpha,
                                         ".",index,".genes",len,".seed",ed,".txt"),
                    quote=F,sep="\t")


        ## create multi dimensional list of genes selected by class
        ## used below for venn diagram 
        if ( l == 1 ){
            selgenes <- list(rownames(t(lasso.select)))
            names(selgenes)[l] <- as.character(results$featureGrouping[[l]])
        } else {
            selgenes <- c(selgenes, list(rownames(t(lasso.select))))
            names(selgenes)[l] <- as.character(results$featureGrouping[[l]])
        }

    } else {
        stop("Number of selected genes do not match the original dataset")
    }

    # get all genes selected from each grouping
    lasso.selected.genes <- rbind(lasso.selected.genes, t(lasso.select))
}


# combine TRUE sample grouping (observed) with predicted
# lasso.select contains the genes that were selected for their
# high probability to assign samples to their correct group
lasso.selected.genes.nodup <- unique(lasso.selected.genes[, ])
dat <- data.frame(y=y, t(lasso.selected.genes.nodup))
dim(dat)
cat("\nGroup imabalances across available all samples: ")
round(table(dat$y)/length(dat$y)*100, 2)

## get gene names
gene.names[ gene.names$genes %in% colnames(dat), ] %>%
    write.table("log.gene.names_bestLambda.txt", sep = '\t', row.names = FALSE, quote = FALSE)

## # plot lambda iterations
## par(mfrow = c(3,4))
## pdf(paste0("grid.lambda.",patient_labels,response,
##            ".regularization",setalpha,".",index,
##            ".genes",len,".seed",ed,".pdf"))
## plot(lasso.trained, xvar="lambda", label=T)
## # fraction deviance explained =R2
## plot(lasso.trained, xvar="dev", label=T)
## plot(cv.out)
## dev.off()

#####################################
## Visualization of selected genes ##
#####################################
## summarize variation between genes (response variables)
## variation explained by the sample grouping (explanatory variables)

## create dummy variables or multi level contrasts
## first level (the baseline) is rolled into the intercept
## all other levels have a coefficient that differ from the baseline
## Helmert regressors compare each level with the average of the preceding ones
## first coefficient is the mean of the first two levels minus the first level
## coefficient is the mean of all three levels minus the mean of the first two levels
set.seed(ed)
associations <- y
contrasts(associations) <- "contr.helmert"
contrasts(associations)

associations <- as.data.frame(model.matrix(~0 + y))
colnames(associations) <- levels(y)

## redundancy analysis variant of canonical correspondence analysis
## implementation based on Legendre & Legendre's (1998) algorithm
## chi-squared data fitted into unweighted linear regression
## on constraining variables (sample groups)
## fitted values are then transformed by singular value decomposition
##rda.results <- vegan::rda( dat[, -1] ~ ., associations[, nl], scale = TRUE )
rda.results <- vegan::rda( adj.x ~ ., associations[, nl], scale = TRUE )
rda.scores <- vegan::scores(rda.results)$species

## variance inflation factor to identify collinearity
## higher the value of the R squared value of the regression between variables, higher the collinearity
vs <- vif.cca(rda.results)

## Tests of significance
## Total variance
tv <- rda.results$CCA$tot.chi / rda.results$tot.chi
ana <- anova(rda.results, step=2000)
ena <- envfit(rda.results ~ ., associations[, nl], perm=2000, dis="lc", scaling=2)

## Akaike information criterion
## AIC estimates relative quality of models
## quality of each model relative to each other is used to select the best one
extended <- rda(adj.x ~ ., associations[, nl])
reduce0 <- rda(adj.x ~ 1, associations[, nl])
reduce <- step(reduce0, scope=formula(extended), test="perm")
reduce
max.aic <- round(max(reduce$anova$AIC),2)
min.aic <- round(min(reduce$anova$AIC),2)

## chart 1
## for scaling check ?biplot.rda {vegan}
## here the scores are scaled symmetrically by square root of eigenvalues
## create color palette
xcol <- length(unique(ids2modules[, clustering.strategy]))
available.colors <- brewer.pal(12, name = "Paired")
selected.colors <- colorRampPalette(available.colors)(n = xcol)
selected.forms <- c(15:19)

pdf(paste0("rda.bestLassoLambda.seed.pdf"))
par(mar=c(1,1,1,1), fig = c(.05,1,.05,1))
plot(rda.results, dis=c("cn","sp"), yaxt="n", scaling=3, type="n", bty = "n")

## project genes from network
axis(2, las = 1)
points(rda.scores, pch=1, col="grey", cex=.5)

## add vectors
## add centroids of factor constraints
##text(rda.results, dis="cn", col="chocolate", font=4)
plot(ena, add = TRUE, col = "chocolate")

## ## show each gene association to a class (form and color)
## for ( sp in nl ) {
##     genes2modules <- ids2modules[ ids2modules$ids %in% selgenes[[sp]], c(1, clustering.strategy)]
##     ## add +1 to offset module 0
##     vcol <- (unique(genes2modules[, 2]) + 1)

##     points(rda.scores[rownames(rda.scores) %in% selgenes[[levels(y)[sp]]],],
##            pch = selected.forms[[sp]],
##            col = selected.colors[vcol], cex=c(2-(sp * 0.5)))
## }

## gene association to classes (colors only)
## distinguish genes selected by best lambda by minimizing the least squares
sel.mods <- colnames(dat[, -1])
genes2modules <- ids2modules[ ids2modules$ids %in% sel.mods, c(1, clustering.strategy)]
colnames(genes2modules) <- c("handle", "mods")
## add +1 to offset module 0
vcol <- (unique(genes2modules[, 2]) + 1)

for (co in vcol) {
    sg <- genes2modules %>%
        filter(mods == co)
    points(rda.scores[rownames(rda.scores) %in% sg$handle, ],
           pch = 16,
           col = selected.colors[vcol])
}

## content of the legend from tests of significance
##legend("bottomleft", colnames(associations), pch = selected.forms[nl], cex = .8, bty = "n")
legend("topleft",
       inset = -.01,
       title = "Sample-wise gene expression statistics",
       c("+ Variance inflation factor (collinearity)",
         paste0(names(vs),": ",round(vs, 2)),
         paste0("+ Total variance: ", round(tv, 2)),
         paste0("+ F-stats: ",round(ana$F[1],2), ", pval ",ana[[4]][1]),
         "+ Akaike information criterion (best AIC): ",
         paste0(max.aic," (",round(100-((min.aic/max.aic)*100),2),"% improvement)"),
         "+ Eigenvalues for constrained axes",
         paste0(names(rda.results$CCA$eig),": ",round(rda.results$CCA$eig, 2), "%"),
         "+ Permutation over RDA1-2 (n=2000)",
         paste0(names(ena$vectors$r),": R2 ",round(ena$vectors$r,2), ", pval ", round(ena$vectors$pvals,5))),
       cex=.7,
       bty = "n",
       horiz=FALSE)

par(fig = c(.7,1,.7,1), new = TRUE, cex = .6)
## Venn diagram showing redundancy between genes assigned to different subsets
## each specifically designed (categorized) to predict a class
try(venn(selgenes), silent = TRUE)
dev.off()
try(dev.off(), silent = TRUE)



## Chart 2
## principal component analysis of gene distribution
## get singular value deomposition which retains
## a reduced number of orthogonal covariates
## that explain as much variance as possibl
genes.pca <- prcomp(t(dat[, -1]), scale = TRUE, center = TRUE)
genes.scores <- as.data.frame(genes.pca$x)[,c(1:3)]
genes.scores$handle <- rownames(genes.scores)

## get pcs
pcs <- t(data.frame(summary(genes.pca)$importance)[2, 1:3])

## plot of observations from networks output
sel.mods <- colnames(dat[, -1])
genes2modules <- ids2modules[ ids2modules$ids %in% sel.mods, c(1, clustering.strategy)]
colnames(genes2modules) <- c("handle", "mods")
genes2modules$handle <- as.character(genes2modules$handle)

pdf("pca.genes2D.bestLassoLambda.pdf")
full_join(genes.scores, genes2modules, by = "handle") %>% 
    ggplot(aes(x = PC1, y = PC2, colour = as.factor(mods))) +
    geom_hline(yintercept = 0, colour = "gray65") +
    geom_vline(xintercept = 0, colour = "gray65") +
    geom_point() +
    xlab(paste0("PC1 (",round(pcs[1]*100,2),"%)")) +
    ylab(paste0("PC2 (",round(pcs[2]*100,2),"%)")) +
    theme(legend.position = "none") +
    theme_minimal()
dev.off()
try(dev.off(), silent = TRUE)


## Chart 2.1
## 3D PCAp
dpf <- full_join(genes.scores, genes2modules, by = "handle")
fit <- lm(dpf$PC3 ~ dpf$PC1 + dpf$PC2)
# predict values on regular xy grid
grid.lines = 4
x.pred <- seq(min(dpf$PC1), max(dpf$PC1), length.out = grid.lines)
y.pred <- seq(min(dpf$PC2), max(dpf$PC2), length.out = grid.lines)
xy <- expand.grid( x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)
fitpoints <- predict(fit)

pdf("pca.genes3D.bestLassoLambda.pdf")
scatter3D(dpf$PC1, dpf$PC2, dpf$PC3,
          col.var = dpf$mods,
          col = selected.colors,
          pch = 19, cex = 1,
          bty = "g", colkey = FALSE,
          theta = 60,
          phi = 5,
          xlab = paste0("PC1 (",round(pcs[1]*100,2),"%)"),
          ylab = paste0("PC2 (",round(pcs[2]*100,2),"%)"),
          zlab = paste0("PC3 (",round(pcs[3]*100,2),"%)"),
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA, fit = fitpoints))
dev.off()
try(dev.off(), silent = TRUE)


## Chart 3
## principal component analysis of sample distribution
## get singular value deomposition
class.pca <- prcomp(dat[, -1], scale = TRUE)
class.scores <- as.data.frame(class.pca$x)[,c(1:3)]
class.scores$handle <- rownames(class.scores)
class.scores$ids <- dat[, 1]
pcs <- t(data.frame(summary(class.pca)$importance)[2, 1:3])

pdf("pca.class2D.bestLassoLambda.pdf")
class.scores %>%
    ggplot(aes(x = PC1, y = PC2, shape = as.factor(ids), color = as.factor(ids))) +
    geom_hline(yintercept = 0, colour = "gray65") +
    geom_vline(xintercept = 0, colour = "gray65") +
    geom_point() +
    xlab(paste0("PC1 (",round(pcs[1]*100,2),"%)")) +
    ylab(paste0("PC2 (",round(pcs[2]*100,2),"%)")) +
    theme(legend.position = "top") +
    theme_minimal()
dev.off()
try(dev.off(), silent = TRUE)


## CHART 4
## Correlation matrices for each group
## show correlation for selected genes from Lasso
## create heatmaps for genes specific for each sample class
vints <- attr(venn(selgenes), "intersections")
save(list=ls(pattern="vints"),file="venn.intersections.Rdata")

pdf("correlation.matrices4venn.intersection.pdf", onefile = TRUE)
for ( lev in 1:length(selgenes) ) {
    ## get the genes that intersect
    df <- adj.x[, selgenes[[lev]]]
    cordf <- round(cor(df),2)

    get_upper_tri <- function(da){
        ## function from the ggplot tutorials
        ## removes half of the correlation matrix
        da[lower.tri(da)] <- NA
        return(da)
    }

    ## cluster and arrange correlation matrix
    dd <- as.dist((1-cordf)/2)
    hc <- hclust(dd)
    cordf <-cordf[hc$order, hc$order]

    ## plot correlation matrix
    corplot <- get_upper_tri(cordf) %>%
        melt(na.rm = TRUE) %>%
        ggplot(aes(x = Var2,
                   y = Var1,
                   fill = value)) +
        geom_tile() +
        scale_fill_gradient2(low = "#67a9cf", high = "#ef8a62", mid = "#f7f7f7", 
                             midpoint = 0, limit = c(-1,1), space = "Lab", 
                             name="Pearson\nCorrelation") +
        theme_minimal()+ 
        coord_fixed() +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.ticks = element_blank(),
            legend.justification = c(1, 0),
            legend.position = c(0.6, 0.7),
            legend.direction = "horizontal") +
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                     title.position = "top", title.hjust = 0.5)) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                         size = 5, hjust = 1)) +
        ggtitle(paste0("Classifier genes for ",names(selgenes)[lev]))

    ## create dendrograms
    dendro.data <-  as.dendrogram(hc) %>%
        dendro_data 
    dendro.plot <- ggplot(segment(dendro.data)) + 
        geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
        theme_minimal()

    ## add dendrograms to the correlation matrices
    grid.newpage()
    print(corplot, vp=viewport(0.8, 0.8, x=0.4, y=0.4))
    print(dendro.plot, vp=viewport(0.665, 0.17, x=0.488, y=0.85))
    print(dendro.plot + coord_flip(), vp=viewport(0.17, 0.667, x=0.85, y=0.435))
}
dev.off()
try(dev.off(), silent = TRUE)




pdf("test2.pdf")
    grid.newpage()
    print(corplot, vp=viewport(0.8, 0.8, x=0.4, y=0.4))
    print(dendro.plot, vp=viewport(0.665, 0.17, x=0.488, y=0.85))
    print(dendro.plot + coord_flip(), vp=viewport(0.17, 0.667, x=0.85, y=0.435))
dev.off()

## CHART 5
## box plot all selected genes
## grouped by module from hierarchical analysis
## before inferring gene network associations
pdf("boxplots.groups2genes.intersection.pdf", onefile = TRUE)
for ( lev in 1:length(selgenes) ) {
    ## get the genes that intersect
    ## add module association from hierarchical clustering
    ## add gene function to gene ids from annotated files generated from PBS script
    genes2modules <- ids2modules[ ids2modules$ids %in% selgenes[[lev]], c(1, clustering.strategy)]
    colnames(genes2modules) <- c("genes", "modules")
    genes2description <- ids2description[ ids2description$genes %in% selgenes[[lev]], ]    

    full.list <- as.data.frame(adj.x[ , selgenes[[lev]] ]) %>%
        mutate(groups = metadata$Groups) %>%
        mutate(nodes = metadata$Nodes) %>%
        mutate(coo = metadata$Prediction) %>%            
        mutate(samples = rownames(adj.x)) %>%
        mutate(category1 = paste0(coo,"-",nodes)) %>%
        gather("genes", "expressions", 1:length(selgenes[[lev]])) %>%
        left_join(genes2modules, by = "genes") %>%
        left_join(genes2description, by = "genes") %>%
        ggplot(aes(x = reorder(paste0(genes,"-",symbol), expressions),
                   y = expressions,
                   color = groups)) +
        geom_jitter(aes(color = factor(groups)),
                    shape=16,
                    position=position_jitterdodge(dodge.width=.8),
                    cex = .15) +
        geom_boxplot(outlier.colour = NA, lwd = .1) +
        coord_flip() +
        facet_wrap(~ coo + nodes,
                   ncol = 6) +
        scale_color_brewer(palette="Dark2") +
        theme_minimal() +
        theme(legend.position = "top",
              text = element_text(size = 5),
              axis.text.y = element_text(size = rel(.5))) +
        ggtitle(paste0("Classifier genes for ",names(selgenes)[lev])) +
        xlab("") +
        ylab("Log2 scaling of expression after RMA quantile normalization (2 is 4-fold up)")

    print(full.list)

}
dev.off()
try(dev.off(), silent = TRUE)



## CHART 6
## gene expression from lmfit data from limma eBayes empirical analysis
pdf("bars.classifier.genes.pdf", onefile = TRUE)
for ( g in grouping ) {
    for ( lev in 1:length(selgenes) ) {
    gene.pvals <- pvals[ pvals$ID %in% selgenes[[lev]],  ]
    print(dim(gene.pvals))

    gene.bars <- gene.pvals %>%
        filter(Contrast == g) %>%
        ggplot(aes(x = reorder(paste0(ID,"-",Symbol), LogFC),
                   y = LogFC,
                   fill = factor(Comparison))) +
        geom_bar(stat = "identity",
                 position = "dodge") +
        coord_flip() +
        scale_fill_brewer(palette="Set1") +
        theme_minimal() +
        theme(legend.position = "top",
              text = element_text(size = 6),
              axis.text.y = element_text(size = rel(1.5)),
              axis.text.x = element_text(size = rel(1.5)),
              strip.background = element_rect(linetype = "blank",
                                              fill = "white"),
              panel.border = element_rect(linetype = "blank",
                                          fill = NA),
              panel.grid.major = element_line(linetype = "blank")) +
        ggtitle(paste0("Differential genes for\n",g, " with classifier genes for ", names(selgenes)[lev])) +
        xlab("") +
        ylab("Log2 scaling of expression after RMA quantile normalization (2 is 4-fold up)")

    print(gene.bars)
    }
}
dev.off()
try(dev.off(), silent = TRUE)



## CHART 7
##  RMA expressions
pdf("boxplots.classifier.genes.pdf")
for ( lev in 1:length(selgenes) ) {
    ## get gene annotations
    genes2description <- ids2description[ ids2description$genes %in% selgenes[[lev]], ]    
    colnames(genes2description) <- c("ID", "chromosome", "ensembl", "symbol", "function", "site", "symbol2")
    ## merge expression and annotations
    rma.selected <- means[ rownames(means) %in% selgenes[[lev]], ]
    rma.selected$ID <- rownames(rma.selected)
    dm <- left_join(rma.selected, genes2description, by = "ID")
    rownames(rma.selected) <- paste0(dm$ID,"-", dm$symbol)
    ## plot
    selgenes.box <- data.frame(y, t(rma.selected[, -c(dim(means)[2]+1)])) %>%
        gather("id", "expression", 2:dim(rma.selected)[1]+1) %>%
        ggplot(aes(x = reorder(paste0(id), expression),
                   y = expression,
                   fill = y)) +
        geom_boxplot(outlier.colour = NA, lwd = .1) +
        coord_flip() +
        scale_color_brewer(palette="Dark2") +
        theme_minimal() +
        theme(legend.position = "top",
              text = element_text(size = 7),
              axis.text.y = element_text(size = rel(.5))) +
        ggtitle(paste0("Classifier genes for ", names(selgenes[lev]))) +
        xlab("") +
        ylab("Log2 scaling of expression after RMA quantile normalization (2 is 4-fold up)")

    print(selgenes.box)
}
dev.off()
try(dev.off(), silent = TRUE)




## CHART 8
## RMA expression without contrast grouping
pdf("boxplots.genes2groups.intersection.pdf", onefile = TRUE)
for ( lev in 1:length(vints) ) {
    genes2description <- ids2description[ ids2description$genes %in% vints[[lev]], ]    
    colnames(genes2description) <- c("genes", "chromosome", "ensembl", "symbol", "function", "site", "symbol2")

    full.list <- as.data.frame(adj.x[ , vints[[lev]] ]) %>%
        mutate(groups = metadata$Groups) %>%
        mutate(nodes = metadata$Nodes) %>%
        mutate(coo = metadata$Prediction) %>%            
        mutate(samples = rownames(adj.x)) %>%
        mutate(category1 = paste0(coo,"-",nodes)) %>%
        gather("genes", "expressions", 1:length(vints[[lev]])) %>%
        left_join(genes2modules, by = "genes") %>%
        left_join(genes2description, by = "genes") %>%
        ggplot(aes(x = reorder(category1, expressions),
                   y = expressions)) +
        geom_boxplot(aes(fill = groups), outlier.colour = NA, lwd = .1) +
        facet_wrap(~ paste0(genes,"-",symbol),
                   ncol = 10) +
        coord_flip() +
        theme_minimal() +
        scale_color_brewer(palette = "Dark2") +            
        theme(legend.position = "top",
              text = element_text(size = 4),
              axis.text.y = element_text(size = rel(.5))) +
        ggtitle(paste0("Intersection between ",names(vints)[lev]," genes")) +
        xlab("") +
        ylab("Log2 scaling of expression after RMA quantile normalization (2 is 4-fold up)")

    print(full.list)
}
dev.off()
try(dev.off(), silent = TRUE)


## CHART 9
## clustering and bootstrap of lasso selected genes
## get adjusted pvalues of similar genes
## analysis will generate heatmaps for each contrast
cluster.samples <- c("pearson")
cluster.genes <- c("pearson")
normalize.genes <- c("standardize")
dissimilar.genes <- c("ward.D2")
n.bootstraps=2
p.alpha=0.95
tree.cut=1.8

pdf("heatmaps.modules4venn.intersection.pdf", onefile = TRUE)
for ( lev in 1:length(selgenes) ) {
    if ( length(selgenes[[lev]]) >= 10 ) {
        gc()
        for ( i in normalize.genes ) {
            for ( it in dissimilar.genes ) {
                cat(">> ",names(selgenes[lev]),"expression normalized with",i,"method &",it,"clustering >> ")
                
                df <- adj.x[ , selgenes[[lev]] ]

                ## add gene function to gene ids
                genes2description <- ids2description[ ids2description$genes %in% selgenes[[lev]], ]    
                colnames(genes2description) <- c("genes", "chromosome", "ensembl", "symbol", "function", "site")

                ## normalize
#                scaledata <- t(decostand(df, method = normalize.genes))
                scaledata <- t(df)
                rownames(scaledata) <- paste0(genes2description$genes, "-", genes2description$symbol)
                gct <- dim(df)[2]

                ## Dissimilarity clustering and tree cutting
                palette.hc <- brewer.pal(11, name = "RdYlBu")
                palette.hclust <- colorRampPalette(palette.hc)(n = c(gct * .05))
                hra <- hclust(as.dist(1-cor(t(scaledata), method= cluster.genes)), method= dissimilar.genes)
                hca <- hclust(as.dist(1-cor(scaledata, method= cluster.samples)), method= dissimilar.genes)
                mycl.row <- cutree(hra, h=max(hra$height) / tree.cut)
                mycl.col <- cutree(hca, h=max(hca$height) / tree.cut)

                ## Color clusters
                maxClusters.row <- length(unique(mycl.row))
                maxClusters.col <- length(unique(mycl.col))                

                if ( maxClusters.row <= 12) {
                    myrowhc <- brewer.pal(maxClusters.row, name = 'Paired')
                } else {
                    myrowhc <- colorRampPalette(brewer.pal(8, name="Dark2"))(maxClusters.row)
                }

                if ( maxClusters.col <= 12) {
                    mycolhc <- brewer.pal(maxClusters.col, name = 'Paired')
                } else {
                    mycolhc <- colorRampPalette(brewer.pal(9, name="Set1"))(maxClusters.col)
                }

                myrowhc <- myrowhc[as.vector(mycl.row)]
                mycolhc <- mycolhc[as.vector(mycl.col)]                

                ## bootstraping and pvalues calculation for each branch
                ## interval confidence (5% chance wrong clustering)
                ## retrieve members of significant clusters
                pvData.row <- pvclust(t(scaledata), method.dist="correlation", method.hclust= it,
                                      nboot= n.bootstraps, parallel=TRUE)
                pvData.col <- pvclust(scaledata, method.dist="correlation", method.hclust= it,
                                      nboot= n.bootstraps, parallel=TRUE)
                clsig.row <- unlist(pvpick(pvData.row, alpha = p.alpha,
                                           pv="au", type="geq", max.only=TRUE)$clusters)
                clsig.col <- unlist(pvpick(pvData.col, alpha = p.alpha,
                                           pv="au", type="geq", max.only=TRUE)$clusters)                

                ## plot heatmap
                ## color dendrograms
                dendroCol <- function(dend=dend, keys=keys, xPar="edgePar", bgr="red", fgr="blue", pch=20, lwd=1, ...) {
                    if(is.leaf(dend)) {
                        myattr <- attributes(dend)
                        if(length(which(keys==myattr$label))==1){
                            attr(dend, xPar) <- c(myattr$edgePar, list(lab.col=fgr, col=fgr, pch=pch, lwd=lwd))
                        } else {
                            attr(dend, xPar) <- c(myattr$edgePar, list(lab.col=bgr, col=bgr, pch=pch, lwd=lwd))
                        }
                    }
                    return(dend)
                }

                dend_colored.row <- dendrapply(as.dendrogram(pvData.row$hclust), dendroCol,
                                               keys=clsig.row, xPar="edgePar", bgr="black", fgr="red", pch=20)
                dend_colored.col <- dendrapply(as.dendrogram(pvData.col$hclust), dendroCol,
                                               keys=clsig.col, xPar="edgePar", bgr="black", fgr="red", pch=20)                

                par(cex.main=.8)
                heatmap.2(scaledata,
                          Rowv=dend_colored.row, Colv=dend_colored.col,
                          col=rev(palette.hc),
                          scale="row", trace="none",
                          RowSideColors = myrowhc, ColSideColors=mycolhc,
                          margin=c(5, 30),
                          cexRow=.2, cexCol=.15,
                          key.title = c("Log2 fold change estimates"),
                          key.ylab = NA,
                          key.xlab = c("Genes Z-scores"),
                          density.info = "none",
                          main = paste0(names(selgenes[lev]),'\n',
                                        "expression normalized with ",i,'\n',
                                        "method & ",it," clustering"),
                          keysize = .9)
            }
        }
    } else {
        print(paste0(names(selgenes[lev]), " interactions has ", length(selgenes[[lev]]),
                     " genes, which is less than the recommended gene counts for clustering."))
    }
}

dev.off()
try(dev.off(), silent = TRUE)





############################
## Validating classifiers ##
############################
# machine learning models used
model_types <- c(
    "nnet", "pcaNNet",
    "dnn",
    "kernelpls", "svmLinear", "svmPoly", "svmRadialSigma", "svmLinear3",
    "lda2", "bagFDA", "fda", "pda", "bagFDAGCV",
    "kknn", "naive_bayes", "gbm",
    "monmlp", "mlpSGD",
    "rf", "RRF",
    "multinom"
)

# number of parameters per model
# the deep network used in this step is an automated model
# hence the low number of parameters to adjust
parameter_counts <- c(
    2,2,
    5,
    1,1,3,2,2,
    1,2,2,1,1,
    3,3,4,
    2,8,
    1,3,
    1
)

## debugging single models
##model_types="svmPoly"
##parameter_counts=3


##########################
######## STEP I ##########
##########################
## setting baseline metrics
## performance metrics for best model without hyperparameter optimization
## create an index for the iteration holder below
## Split the dataset
icc <- function(){ i=0; function(){ i <<- i + 1;  i }}
modet=ite=NULL
ite=icc()
models=icc()

for ( epochs in c(1:baseline.epochs) ) {
    ## as many iterations will be executed for each model
    ## iterations are done in addition to 25 resampling for each model
    ie=ite()
    edp <- floor(abs(rnorm(1) * 10000000))
    set.seed(edp)
    base.training <- sample(1:nrow(adj.x), c(nrow(adj.x) * c(train.size/100)) )

    for ( m in 1:length(model_types) ) {
        mods=model_types[[m]]
        param=parameter_counts[[m]]

        ## models are trained in succession
        ## output is saved
        start <- format(Sys.time(), "%a-%d %X")
        cat(">> Iteration", ie, "on model", mods, "started at", start)
        modnam=models()
        model.name <- paste0(mods,"|",ie,"|",param)
        
        trainCtrl <- trainControl(method="cv",
                                  number=c(validate.folds / 2))
        trained.model <- train( y ~ .,
                               data = dat[base.training, ],
                               trControl= trainCtrl,
                               preProc=c("center","scale"),
                               method = mods )

        # aggregate all performance metrics
        if ( modnam == 1 ){
            performance_summary <- list(trained.model)
            names(performance_summary)[modnam] <- model.name
        } else if ( modnam > 1 ) {
            performance_summary <- c(performance_summary, list(trained.model))
            names(performance_summary)[modnam] <- model.name        
        }
        end <- format(Sys.time(), "%a-%d %X")
        cat(". >", mods, "execution successful at", end, "\n")
    }
}

###############
## summary 1 ##
###############
## performance metrics between models with repeated iterations
sink(paste0("performance1.multianalysis.seed",ed))
performance_summary %>% resamples %>% summary
sink()

## save output
save(list=ls(pattern="perf*summary"),file="performanceSummaryNoTune.Rdata")

if ( file.exists(paste0("performance1.multianalysis.seed",ed)) ) {
    sink("3-baseline.training.wo.tuning.ok")
    sink()
}



###########################
######## STEP II ##########
###########################
## improving the baseline metrics
## Classification across models with hyperparameter optimization
## with hyperparameter tuning
if ( classification == TRUE & grouped == TRUE ) {
    importance <- NULL
    systems.metrics <- NULL
    subtrain.metrics <- NULL

    set.seed(ed)
    training <- sample(1:nrow(adj.x), c(nrow(adj.x) * c(train.size/100)) )

    for ( m in 1:length(model_types) ) {

        mods=model_types[[m]]
        param=parameter_counts[[m]]

        ## start logging
        start <- format(Sys.time(), "%a-%d %X")
        cat("\n>> Training on model", mods, "started at", start)

        ## models are trained in succession
        ## predicted output is saved
        model.name <- paste0(mods,"|",param)

        ## multimodel analysis
        ## with summary output
        model.metrics <- modelTune.clas(dat, training, method=mods,
                                        folds=validate.folds, grid=validate.tune)

        ## get predictor importance
        ## accuracy of a predictor is calculated to detect a decrease after permutation without it
        ## reported as class-specific decreases in accuracy
        if ( mods != "bagFDA" && mods != "fda" &&
             mods != "bagFDAGCV" && mods != "gbm" &&
             mods != "mlpSGD" && mods != "rf" &&
             mods != "RRF" && mods != "multinom") {
        bestmod <- varImp(model.metrics$bestModel)
        importance <- rbind(importance, data.frame(model = mods,
                                                   round(bestmod[[1]][, levels(y)], 2)))
        }

        ## aggregate all performance metrics
        ## only for predicted features
        ## contains results and tuned/selected parameters per model
        if ( m == 1 ){
            ## extract/aggregate best model accuracy counts
            performance_summaryFull <- list(model.metrics$bestModel)
            names(performance_summaryFull)[m] <- model.name
            ## extract/aggregate best model hyperparameters
            parameters_summaryFull <- list(model.metrics$Hyperparameters)
            names(parameters_summaryFull)[m] <- mods
        } else if ( m > 1 ) {
            performance_summaryFull <- c(performance_summaryFull, list(model.metrics$bestModel))
            names(performance_summaryFull)[m] <- model.name
            parameters_summaryFull <- c(parameters_summaryFull, list(model.metrics$Hyperparameters))
            names(parameters_summaryFull)[m] <- mods
        }

        ## get iteration speed
        ## get sensitivity, specificity, precision scores
        ## for models trained and validated 
        durationMinutes <- round(model.metrics$timeLapsed[[3]]/60, 2)
        systems.metrics <- rbind(systems.metrics,
                                 data.frame(model=mods,
                                            durationMinutes=durationMinutes,
                                            model.metrics$ConfusionMatrix$byClass,
                                            accuracy=model.metrics$ConfusionMatrix$overall[[1]],
                                            kappa=model.metrics$ConfusionMatrix$overall[[2]],
                                            accuracyPval=model.metrics$ConfusionMatrix$overall[[6]]))
        ## for models validated (only) at iterated epochs
        subtrain.metrics <- rbind(subtrain.metrics, model.metrics$metrics)

        ## run best model and get probability scores
        if ( mods == "svmPoly" ) {
            class.prob.bestmod <- svmPoly.prob(dat, training)
        }
        
        ## end logging
        end <- format(Sys.time(), "%a-%d %X")
        cat(paste0("\n>> ",mods," execution successful at ", end,
                   " - Duration: ",durationMinutes,"min",
                   " (",round(durationMinutes/60,2),"H)\n"))
    }

} else if ( classification == TRUE & binomial == TRUE ) {

    output_summary <- modelTune.clas(dat,training,method="rf",folds=10, grid=FALSE)

} else if ( regression == TRUE & bionomial == TRUE ) {

    output_summary <- modelTune.reg(dat,training,method="nnet",folds=10,ctl)

}

###############
## summary 2 ##
###############
## performance metrics for best model while hyperparameters tuning
sink(paste0("performance2.hyperTuning.seed",ed))
performance_summaryFull %>% resamples %>% summary
sink()

###############
## summary 3 ##
###############
## recall, sensitivity, specificity, etc. for best parameters during training (only)
write.table(systems.metrics,
            paste0("performance3.full.hyperTuning.seed",ed),
            sep = "\t", quote = FALSE)

###############
## summary 4 ##
###############
## recall, sensitivity, specificity, etc. for best parameters during sub-training (only)
write.table(subtrain.metrics,
            paste0("performance4.mini-batch.metrics.seed",ed),
            sep = "\t", quote = FALSE)

###############
## summary 5 ##
###############
## importance of variables
write.table(importance,
            paste0("performance5.importance.seed",ed),
            sep = "\t", quote = FALSE)

###############
## summary 6 ##
###############
## probabilities for each predicted sample by class
write.table(class.prob.bestmod,
            paste0("performance6.class_probabilities.seed",ed),
            sep = "\t", quote = FALSE)

## save classification output
save(list=ls(pattern="*metrics"),file="systemsMetrics.Rdata")
save(list=ls(pattern="*summaryFull"),file="performanceSummaryFull.Rdata")


if ( file.exists(paste0("performance2.hyperTuning.seed",ed)) ) {
    sink("4-model.training.w.tuning.ok")
    sink()
}


#####################
##### debugging #####
### new functions ###
#####################
modelTune.clas.debug <- function(dat, train, method, folds=10, grid=TRUE, confusion.metrics=NULL){
    ##    trainCtrl <- trainControl(method="repeatedcv",number=folds,repeats=rep,summaryFunction=defaultSummary)
        trainCtrl <- trainControl(method="repeatedcv",number=folds,repeats=c(folds/2),classProbs=T)
    if ( grid == TRUE ) {
        if ( method == "gbm_h2o" ) {
            grid_models <- expand.grid(.ntrees=seq(2,30,length=2),
                                       .max_depth=5,
                                       .min_rows=0.7,
                                       .learn_rate=0.5,
                                       .col_sample_rate=0.8)
        } else if ( method == "nnet" ) {
            ## neural networks
            grid_models <- expand.grid(.size=seq(1,5,length=10), .decay=10^seq(-1,-2,length=10))
        } else if ( method == "svmPoly") {
            ## svmRadial
            grid_models <- expand.grid(.degree=seq(0,10,.5),
                                       .scale=10^seq(-1,-3,length=10),
                                       .C=seq(.1,2, length=20))
        } else if ( method == "LogitBoost" ) {
            ## Boosted logistic regression
            grid_models <- expand.grid(.nIter=seq(1,100,length=2))
        }

        lapsed <- system.time(modelTrain <- train(y~.,
                                                  data=dat[train,],
                                                  method=method,
                                                  trControl= trainCtrl,
                                                  preProc=c("center","scale"),
                                                  tuneGrid=grid_models,
##                                                  weights = model_weights,
                                                  tuneLength=c(folds*3)))
    } else if ( grid == FALSE ) {
        lapsed <- system.time(modelTrain <- train(y~., data=dat[train,],
                                                  method=method,trControl= trainCtrl,
                                                  preProc=c("center","scale")))  }

    for ( iterations in c(1: c(folds * 2)) ) {
        ## iterative prediction test for model performance
        edx <- floor(abs(rnorm(1) * 10000000))
        set.seed(edx)
        subtrain <- sample( 1:nrow(dat), c(nrow(dat) * c(train.size/100)) )
        cat("\n  >> Sub-iteration", iterations, "on model", method, "started")
        ## predict over best hyperparameters
        end <- system.time(Predd <- predict(modelTrain, newdata=dat[-subtrain,], type="raw"))
        conf.m <- confusionMatrix(data=Predd, dat[-subtrain,1])
        ## aggregate all prediction metrics
        confusion.metrics <- rbind(confusion.metrics,
                                   data.frame(iteration=iterations,
                                              seed=edx,
                                              conf.m$byClass,
                                              accuracy=conf.m$overall[[1]],
                                              accLow=conf.m$overall[[3]],
                                              accHigh=conf.m$overall[[4]],
                                              kappa=conf.m$overall[[2]],
                                              accPval=conf.m$overall[[6]]))
        
        cat(" >", method, "validated successfully in", c(end[[1]]*1000), "ms @accuracy", conf.m$overall[[1]])
    }

    Predd <- predict(modelTrain, newdata=dat[-train,], type="raw")
    conf.m <- confusionMatrix(data=Predd, dat[-train,1])
    output <- list(timeLapsed=lapsed,bestModel=modelTrain,
                   Results=modelTrain$results,Hyperparameters=modelTrain$bestTune,
                   ConfusionMatrix=conf.m,
                   metrics=confusion.metrics)
    return(output)
}

## unbalance.class <- table(y[training])/length(y[training])
## model_weights <- as.numeric(1 - unbalance.class[y[training]])


## model_weights <- ifelse(dat$y == "NOREL",
## (1/table(dat$y)[1]) * 0.35,
## (1/table(dat$y)[3]) * 0.5)

## mm.svm <- modelTune.clas.debug(dat,training,method="svmPoly",folds=2, grid=FALSE)

## ## importance
## bestmod <- varImp(mm.svm$bestModel)
## pdf("test.pdf")
## plot(bestmod, size =.5)
## dev.off()


## output <- svmPoly.prob(dat, train)

## model_list <- list(original = mm.o$bestModel,
##                   weighted.nnet = mm.nnet$bestModel,
##                   weighted.linear = mm.linear$bestModel,
##                   weighted = mm.w$bestModel)
## model_list %>% resamples %>% summary


## #### experimental ######
## compare_models(nnet0,rf1)
## model.reg(dat,training,method="rf",folds=10,r=5,tune=10)
## modelTrain <- train(y~., data = dat, method = "glmnet_h2o", classProbs=TRUE)






####### Save ########
## print full view of variables, data frames, matrices, funcitons...
lsos()
sink("log.R.sessionInfo.txt")
sessionInfo()
sink()

##load("EnsembleMethods.Rdata", .GlobalEnv)




############## add the gist to the modelTune.clas function before training the model
############## the following will assign weights to samples (method 1)
##        ## reassign weights for imbalanced samples
##        ## some model do not accept weights though
##        if ( sample.classes == "Groups" ) {
##            ## up-weighting imbalanced samples
##            model_weights <- ifelse(dat$y == "NOREL",
##            (1/table(dat$y)[1]) * 0.35,
##            (1/table(dat$y)[3]) * 0.5)
##            }


############## add the gist to the while loop before weighting
## make sure all sample categories are included
## in the training and testing sets
## if not, errors occur during training
## for missing observations in certain classes
## while (
## (length(unique(y[training])) != nlevels(y))
## &
## (length(unique(y[-training])) != nlevels(y))
## ) {
##     ede <- floor(abs(rnorm(1) * 10000000))
##     set.seed(ede)
##     # split 70% of the data
##     randomize <- sample(nrow(adj.x))
##     training <- sample(randomize, nrow(adj.x)/1.429)
##     #table(y[-training])
##     tr <- length(training)
## }



############### add it to the while loop before regularization
############### the following will assign weights to imbalanced samples (method 2)
## ## treating sample imbalances
## ## by manually assigning CV folds
## ## and by up-weithing over represented samples for shrinkage
## if ( sample.classes == "Groups" ) {
##     ## assign folds evenly using the modulus operator
##     foldid <- as.numeric(length(y[training]))
##     fold0 <- sample.int(sum(y[training] == "CNS")) %% nfold
##     fold1 <- sample.int(sum(y[training] == "SYST")) %% nfold
##     fold2 <- sample.int(sum(y[training] == "NOREL")) %% nfold
##     foldid[ y[training] == "CNS" ] <- fold0
##     foldid[ y[training] == "SYST" ] <- fold1
##     foldid[ y[training] == "NOREL" ] <- fold2
##     foldid.custom <- foldid + 1

##     ## assign weights
##     unbalance.class <- table(y[training])/length(y[training])
##     weights.class <- 1 - unbalance.class[y[training]]

## } else if ( sample.classes == "Contrast2" ) {
##     ## assign folds evenly using the modulus operator
##     foldid <- as.numeric(length(y[training]))
##     fold0 <- sample.int(sum(y[training] == "CNS.EN")) %% nfold
##     fold1 <- sample.int(sum(y[training] == "SYST.EN")) %% nfold
##     fold2 <- sample.int(sum(y[training] == "NOREL.EN")) %% nfold
##     fold3 <- sample.int(sum(y[training] == "CNS.LN")) %% nfold    
##     fold4 <- sample.int(sum(y[training] == "SYST.LN")) %% nfold    
##     fold5 <- sample.int(sum(y[training] == "NOREL.LN")) %% nfold    
##     foldid[ y[training] == "CNS.EN" ] <- fold0
##     foldid[ y[training] == "SYST.EN" ] <- fold1
##     foldid[ y[training] == "NOREL.EN" ] <- fold2
##     foldid[ y[training] == "CNS.LN" ] <- fold3
##     foldid[ y[training] == "SYST.LN" ] <- fold4
##     foldid[ y[training] == "NOREL.LN" ] <- fold5
##     foldid.custom <- foldid + 1

##     ## assign weights
##     unbalance.class <- table(y[training])/length(y[training])
##     weights.class <- 1 - unbalance.class[y[training]]
## }
