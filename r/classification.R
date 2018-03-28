pkgs <- c('RColorBrewer', 'pvclust', 'gplots', 'vegan',
          'dplyr', 'mRMRe', 'glmnet', 'caret', 'foreach',
          'doSNOW', 'lattice', 'ROCR')
lapply(pkgs, require, character.only = TRUE)

### DEFINE DATA FITTING MODEL
classification=TRUE
regression=FALSE

## define feature structure
grouped=TRUE
binomial=FALSE



## load expression data
# optimized for t-statistics microarray expressions
# rows=genes
# col=samples
cat("\n\nNormalized expression scores: Samples are columns and genes are rows\n")
means <- read.table("expressions.149444", sep="\t", header=T, row.names=1)
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


## prepare testing dataset


# make sure all sample categories are included
# in the training and testing sets
while (
(length(unique(y[training])) != nlevels(y))
&
(length(unique(y[-training])) != nlevels(y))
) {

    # set seed for reproducibility
    ed <- floor(abs(rnorm(1) * 10000000))
    set.seed(ed)

    # Split the dataset into 80% training data
    training <- sample(1:nrow(adj.x), nrow(adj.x)/1.25)
    tr <- length(training)
}


## choose samples strutcture from metadata
y <- metadata$Groups

if ( nlevels(y) > 2 ){
    response="multinomial"
} else if ( nlevels(y) == 2 ){
    response="bionomial"
} else {
    response="mgaussian"
}



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
# fit a generalized linear model via penalized maximum likelihood
# if alpha 1 then lasso L1 penality and discard genes
# if alpha 0 then ridge regression then L2 and rank genes
if ( grouped == TRUE ){

    index="grouped"
    ncv=10
    setalpha=1
    
    # fitting a symmetric multinomial model,
    grid <- 10^seq(5, -5, length=200)
    lasso.trained <- glmnet(adj.x[training,],
                            y[training],
                            alpha=setalpha,
                            lambda=grid,
                            family = response,
                            standardize=F,
                            type.multinomial=index)

    pdf(paste0("grid.lambda.",response,".regularization",setalpha,".",index,".features.",ed,".pdf"))
    par(mfrow = c(2,2))
    plot(lasso.trained, xvar="lambda", label=T)
    ## fraction deviance explained =R2
    plot(lasso.trained, xvar="dev", label=T)
    dev.off()


    cv.out <- cv.glmnet(adj.x[training,],
                        y[training],
                        alpha=setalpha,
                        family = response,
                        standardize=F,
                        nfolds = ncv,
                        type.multinomial=index)


    ## Cross validation for hyperparameter tuning.
    # Optimization of model selection to avoid overfitting
    bestlam <- cv.out$lambda.min
    bestlam
    pdf(paste0("cv",ncv,".lambda",sprintf("%.5f", bestlam),".",response,
               ".regularization",setalpha,".",index,".features.",ed,".pdf"))
    plot(cv.out)
    dev.off()

} else {

    stop("Alternative feature indexes are not yet implemented")

}


# get genes probabilities to predict sample cases
lasso.link <- predict(lasso.trained, s=bestlam, newx=adj.x[-training,], type="link")
lasso.response <- predict(lasso.trained, s=bestlam, newx=adj.x[-training,], type="response")
lp <- data.frame(y[-training], lasso.link, lasso.response)
colnames(lp) <- c("labels",
                  paste0(colnames(lasso.link), "-link"),
                  paste0(colnames(lasso.response), "-response"))




# plot multiclass ROC curves
pdf("ROC.pdf")
for ( i in 1:nlevels(y) ) {
    
    couleurs <- brewer.pal(nlevels(y), name = 'Dark2')

    # get scores
    probability.scores <- lp[, c(5+i)]
    dummy.labels <- as.vector(model.matrix(~0 + y[-training])[, i])

    # create list for average ROC
    if ( i == 1 ) {
        listofprobs <- list(probability.scores)
        listofdummies <- list(dummy.labels)        
    } else {
        listofprobs <- c(listofprobs, list(probability.scores))
        listofdummies <- c(listofdummies, list(dummy.labels))    
    }

    # rename iterations
    names(listofprobs)[i]=paste0(levels(y)[i],"-",i)
    names(listofdummies)[i]=paste0(levels(y)[i],"-",i)

    # get accuracy
    pred <- prediction(probability.scores, dummy.labels)
    perf <- performance(pred, 'tpr', 'fpr')

    plot(perf, lwd=1.5, col=couleurs[i],
         xlab="False positive rate (1-Specificity)",
         ylab="True positive rate (Sensitivity)"
         )
    par(new=TRUE)
}
# create average
# error bars for variation around the average curve
pred <- prediction(listofprobs, listofdummies)
perf <- performance(pred, 'tpr', 'fpr')
plot(perf,lty=3,lwd=2.5,avg="vertical",spread.estimate="stderror",add=TRUE)
legend("bottomright", levels(y), lty=1, lwd=5, col = couleurs[1:nlevels(y)])
dev.off()




# get sample labels
lasso.labels <- predict(lasso.trained, s=bestlam, newx=adj.x[-training,], type="class")

# get gene coefficients at selected lambda
lasso.coef <- predict(lasso.trained, s=bestlam, type = "coefficients")
str(lasso.coef)

## get non-zero genes
selected.genes <- lasso.coef[[1]]@i[ lasso.coef[[1]]@i >= 1]
original.genes <- colnames(adj.x)
selected.final <- original.genes[selected.genes]
len <- length(selected.genes)

## extract expression of regularized genes
if ( length(selected.final) == length(selected.genes) ) {
    lasso.select <- adj.x[, selected.final ]
    dim(lasso.select)
    write.table(lasso.select, paste0("expressions.cv",ncv,
                                    ".lambda",sprintf("%.5f", bestlam),
                                    ".",response,".regularization",
                                    setalpha,".",index,".features.",ed,".txt"),
                quote=F,sep="\t")
} else {
    stop("Number of selected genes do not match the original dataset")
}


## create a summary table
if ( classification == TRUE ) {
    ## build classification confusion matrix
    tab <- table(lasso.labels, y[-training])
    tab
    freq <- as.data.frame.matrix(tab)

    df=NULL
    for ( i in 1:nrow(freq) ) {
        ## prepare a table summary
        ## of regularization accuracy
        n=names(rowSums(freq)[i])
        f <- freq[n,n] / rowSums(freq)[[i]] * 100

        if ( setalpha == 1 ) {me="lasso"} else {me="ridge"}
        if ( index == "grouped") {ind=TRUE} else (ind=FALSE)

        df <- rbind(df, data.frame(group=n,
                                   accuracy=f,
                                   seed=ed,
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




##########################
## Define new functions ##
##########################

model.reg <- function(dat,train,test,method,folds=10,rep=5,tune){
    ## requires caret
    ## Regression
    ## Train model and test on independant dataset
    ## returns RMSE
    trainCtrl <- trainControl(method="repeatedcv",number=folds, repeats=rep)
    lapsed <- system.time(modelTrain <- train(y~., data=dat[train,],
                                              method=method,
                                              trControl= trainCtrl,
                                              preProc=c("center","scale"),
                                              tuneLength=tune ))
    ploted <- plot(modelTrain)
    Predd <- predict(modelTrain, newdata=dat[test,], type="raw")
    ## Test set MSE for regression
    rmse <- mean((Predd - y[test])^2)		
    output <- list(ploted,TimeLapsed=lapsed,Prediction.Estimates=Predd,Hyperparameters=modelTrain$bestTune, RMSE=rmse)
    return(output)
}

modelTune.reg <- function(dat,train,test,method,folds=10,rep=5,tune,ctl){
    ## requires caret
    ## Regression
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
    Predd <- predict(modelTrain, newdata=dat[test,], type="raw")
    ## Test set MSE for regression
    rmse <- mean((Predd - y[test])^2)
    output <- list(ploted,TimeLapsed=lapsed,Prediction.Estimates=Predd, Hyperparameters=modelTrain$bestTune, RMSE=rmse)
    return(output)
}

modelTune.clas <- function(dat, train, test, method, folds=10, rep=5, tune=10, grid=ctl){
    ## requires caret
    ## Classification
    ## GRID search HYPERPARAMETERS tuning
    ## Train model and test on independant dataset
    ## returns a classification error
    trainCtrl <- trainControl(method="repeatedcv",number=folds, repeats=rep, , classProbs=T,summaryFunction=defaultSummary)
    lapsed <- system.time(modelTrain <- train(y~., data=dat[train,],
                                              method=method,
                                              trControl= trainCtrl,
                                              preProc=c("center","scale"),
                                              tuneGrid=grid,
                                              tuneLength=tune,
                                              metric="ROC"))
    ploted <- plot(modelTrain)
    Predd <- predict(modelTrain, newdata=dat[test,], type="raw")
    ## confusion matrix for classification
    conf.m <- confusionMatrix(data=Predd, dat[test,1])
    Probs <- predict(modelTrain, newdata=dat[test,], type="prob")
    output <- list(ploted,TimeLapsed=lapsed, Hyperparameters=modelTrain$bestTune, ConfusionMatrix=conf.m,Probabilities=Probs)
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




############################
## TESTING THE CLASSIFIER ##
############################
## Choosing the right Hyper-parameters. GRID ANALYSIS
dat <- data.frame(y=associations, lasso.select); dim(dat)

ctl=expand.grid(.size=seq(1,20,length=40), .decay=10^seq(-1,-5,length=40))	## nnet
#ctl=expand.grid(.mtry=seq(1:15),.coefReg=10^seq(-1,-3,length=40),.coefImp=10^seq(-1,-2,length=40))## Regularized Random forest (RRF)
#ctl=expand.grid(.C=seq(1:20), .sigma=10^seq(-1,-3,length=40))	## svmRadial
#ctl=expand.grid(.ncomp=seq(1:15))	# PCR
#ctl=expand.grid(.mstop=seq(10:1000,length=20),.prune=no)	## glmboost
#ctl=expand.grid(.lambda=10^seq(10,-2,length=100))	## Ridge

model.reg(dat,training,-training,method="RRF",folds=10,r=5,tune=10)

# Tune hyper-parameters
nnet0 <- modelTune.reg(dat,training,-training,method="nnet",folds=10,r=5,tune=10,ctl)

nnet0 <- modelTune.clas(dat,training,-training,method="nnet",folds=10,r=5,tune, grid=ctl)


modelsRMSE <- read.table("clipboard", sep="\t", header=T);modelsRMSE
library(lattice)
xyplot(Gene1 + Gene2 + Gene3 + Gene4 + Gene5~Model, data=modelsRMSE, type=c("a","p"),cex = 1.5, pch=21,auto.key = list(columns=5,points=F,lines=T,title='Uncorrelated genes'),ylab = "Root-mean-square error (RMSE)",xlab="Tested base Learners")
## plot genes 1 to 5 from mRMR
xyplot(subset1 + subset2 + subset3 + subset4 ~Model, data=modelsRMSE, type=c("a","p"), pch =21,cex = 1.5,auto.key = list(columns=4,points=F,lines=T,title='Correlated subgenes'),ylab = "Root-mean-square error (RMSE)",xlab="Tested base Learners")
## plot subset trees to gene1 from mRMR
xyplot(iter100 + iter200 + iter600 + iter1000 + iter2000 + iter5000 ~ Model, data=modelsRMSE, type=c("a","p"),pch=20,cex=1,auto.key = list(space="top",points=F,lines=T),ylab = "Bagging Iteration
s",xlab="Tested base Learners")
## plot bagging iteration RMSEs for base learners
xyplot(iter100 + iter200 + iter600 + iter1000 + iter2000 + iter5000 ~ Model, data=modelsRMSE, type=c("a","p"),cex=1,lty=1, auto.key = list(columns=6,points=T,lines=F,title="Number of iterations"),ylab = "System iteration time (s)",xlab="Tested base Learners")
## plot bagging run time RMSEs for base learners


##############################
# Ensemble Learning (bagging)
##############################
require(caret)
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
# ANOVA for 29 genes
##############################

dat <- t(mmcDat[,-c(1,2)])
diet <- gl(3,9,27, label=c("No","Co","Ti"));diet <- diet[-c(7:9,26:27)];diet
stages <- gl(5,3,28, label=c("E", "T", "V", "P", "J")); stages <- stages[-c(16:21)]; stages
## prepare variables
resDat <- NULL
for(i in 1:ncol(dat)){
test <- data.frame(gene=dat[,i],diet,stages)
resDat <- rbind(resDat, coefficients(aov(gene~diet+stages+diet*stages,test)))
}
## compute residuals for manova
rownames(resDat) <- colnames(dat)
asetwd("C:/Workshop2014/Paper3")
write.csv(resDat, "29residuals.csv",quote=F)
## extract and save

##############################
# Circos
##############################
setwd("C:/workshop2014/Paper3")
#save(list=ls(pattern="locus|opt|setup|mmc"),file="circos_MS3.Rdata")	## save
load("circos_MS3.Rdata", .GlobalEnv)
lsos(pat="locus|opt|setup|mmc")

require(circlize)
circos.test(mmcCorr,5)
## plot correlations for 5 genes

## old
circos.test(optS3,5)
## optS3= all correlations of setup III from the MMC output file
## 5= number of genes to be ploted

##############################
# Tryouts
##############################

rats <- data.frame(id = paste0("rat",1:10),
  sex = factor(rep(c("female","male"),each=5)),
  weight = c(2,4,1,11,18,12,7,12,19,20),
  length = c(100,105,115,130,95,150,165,180,190,175))

## working
mat <- t(optall[1:3,1:5,drop=F])
#rownames(mat) = letters[1:3]
#colnames(mat) = LETTERS[1:6]
rn = rownames(mat)
cn = colnames(mat)
mat
factors = c(rn,cn)
factors = factor(factors, levels = factors)
col_sum = apply(mat, 2, sum)
row_sum = apply(mat, 1, sum)
xlim = cbind(rep(0, ncol(mat)+nrow(mat)), c(row_sum, col_sum))
par(mar = c(1, 1, 1, 1))
circos.par(cell.padding = c(0, 0, 0, 0), clock.wise = FALSE,gap.degree = c(ncol(mat)+nrow(mat)), start.degree = 5)
circos.initialize(factors = factors, xlim = xlim,sector.width = 2)
#circos.clear()
circos.trackPlotRegion(factors = factors, ylim = c(0, 1), bg.border = NA,
bg.col = c("red", "green", "blue", rep("grey", 6)), track.height = 0.05,
panel.fun = function(x, y) {
sector.name = get.cell.meta.data("sector.index")
xlim = get.cell.meta.data("xlim")
circos.text(mean(xlim), 1.5, sector.name, adj = c(0.5, 0))
if(sector.name %in% rn) {
for(i in seq_len(ncol(mat))) {
circos.lines(rep(sum(mat[sector.name, seq_len(i)]), 2), c(0, 1),
col = "white")
}
} else if(sector.name %in% cn) {
for(i in seq_len(nrow(mat))) {
circos.lines(rep(sum(mat[ seq_len(i), sector.name]), 2), c(0, 1),
col = "white")}}})
col = c("#FF000020", "#00FF0020", "#0000FF20")
for(i in seq_len(nrow(mat))) {
for(j in seq_len(ncol(mat))) {
circos.link(rn[i], c(sum(mat[i, seq_len(j-1)]), sum(mat[i, seq_len(j)])),
cn[j], c(sum(mat[seq_len(i-1), j]), sum(mat[seq_len(i), j])),
col = col[i], border = "white")}}
circos.clear()
## Build Circos



require(ggmap)
require(rworldmap)
par(mar=c(2,2,2,.5))
newWorld <- getMap(resolution="high")
plot(newWorld)
geocode("magdalen islands")
plot(newWorld, xlim=c(5,53), ylim=c(17.6,23),asp=1)
## world map


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

##############################
# SAVE
##############################

setwd("C:/Dropbox/Workshop2013/Work/R/ANN")
lsos(pat="locus.select|mRMR|grid|bag")
save(list=ls(pattern="mRMR|grid|bag"),file="EnsembleMethods.Rdata")	## save
save(list=ls(pattern="*.mRMR"),file="lassoSelected.Rdata")	## save
save(list=ls(pattern="locus.select"),file="mRMRselected.Rdata")	## save
load("EnsembleMethods.Rdata", .GlobalEnv)
