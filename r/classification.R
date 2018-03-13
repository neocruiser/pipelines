pkgs <- c('RColorBrewer', 'pvclust', 'gplots', 'vegan')
lapply(pkgs, require, character.only = TRUE)

## import custom functions
source('01funcs.R')

## load expression data
## optimized for t-statistics microarray expressions
means <- read.table("means22samples15K.txt", sep="\t", header=T, row.names=3)
x <- t(means[,-c(1,2)])
dim(x)

xs <- decostand(t(x), "standardize")

## prepare features
y <- c(rep("E",3), rep("T",3), rep("VC",3),rep("PC",3),rep("JC",3),rep("VT",3),rep("PT",3), rep("JT",1))
y <- c(rep("E",3), rep("L",9), rep("PL",3), rep("L",6),rep("PL",1))
y <- c(rep("E",3), rep("L",9), rep("PL",3), rep("Lc",6),rep("PLc",1))
y <- c(rep("Healthy", 15), rep("Deficient",7))
y <- c(rep("None", 6), rep("Cocktail", 9), rep("Tiso",7))
y <- c(rep("E",3), rep("Lc",6), rep("PLc",6), rep("Lt",3),rep("PLt",4))
y <- c(rep("L",9), rep("PL",6), rep("L",3),rep("PL",4))


##############################
# Filter subset selection
##############################
## Unsupervised gene selection based on high variance
## discard extreme variance
nset=5000

hi.x <- set.var(t(xs), 1, nset)
summary(apply(hi.x, 2, var))
summary(apply(hi.x, 2, sd))

cor.x <- cor(hi.x)
summary(cor.x[upper.tri(cor.x)])


## get number of discarded high variance genes
## iteration multiple thresholds
gv <- NULL
for (i in seq(500, 10000, 500)) {
hi.x <- set.var(t(xs), 1, i)

# get mean (dm) and maximum value of variance (dmv) of the whole dataset
dm <- summary(apply(hi.x, 2, var))[[4]]
dmvb <- summary(apply(hi.x, 2, var))[[6]]
offset = c( (dm * (log(dmv)) + (dmv * 0.1)) )

# get the mean for each gene
selected <- data.frame(locus=colnames(hi.x),var=apply(hi.x, 2, var)) %>%
    filter(var > offset)

gv <- rbind(gv, data.frame(dimension=i,
                           meanVariance=dm,
                           maxVariance=dmv,
                           discarded=dim(selected)[1]))
}


require(mRMRe)
feature.select <- new("mRMRe.Data",data=data.frame(hi.x[,1:1500, drop=F]))
## extract data
set.seed(1445612321)
locus.select <- new("mRMRe.Filter", data=feature.select, target_indices=1, levels=c(500,1), continuous_estimator="spearman")
## feature select min redundant max relevant genes
range_correlation(locus.select,n=1,t=1,hi.x,method="spearman")	# n=gene(set w target_indices); t=different gene mashups
## view range of correlated selected features
locus <- locusRMR(locus.select,feature.select,n=1,t=1)
length(locus)
## extract LOCUS names
## Parallelized mRMR ensemble Feature selection

##############################
# Training/Testing
##############################
foo2 <- hi.x[,locus];dim(foo2)
train <- sample(1:nrow(foo2), 2*nrow(foo2)/2.5); length(train)
test <- -train
y[test]; length(y[test])
## prepare training/testing sets
set.seed(1445612321)
dat <- data.frame(y=y, foo2)
y <- as.vector(model.matrix(~y,dat)[,2])
#y <- c(rep(0,6),rep(1,9),rep(-1,7))	## dummy variables for 3 levels None, Coc, Tiso

##############################
# Feature extraction
##############################

## LASSO
grid <- 10^seq(10, -2, length=100)
require(glmnet)
set.seed(1445612321)
lasso.mod <- glmnet(foo2[train,],y[train], alpha=0.4, lambda=grid, family = "gaussian", standardize=T, type.gaussian = "naive")
## (1) Model Training. Select the kernel function and associated kernel parameters
#coef(lasso.mod, s=0)
#par(mfrow = c(2,2));
plot(lasso.mod, xvar="lambda", label=T)
plot(lasso.mod, xvar="dev", label=T) 	## fraction deviance explained =R2
system.time(cv.out <- cv.glmnet(foo2[train,], y[train], alpha=0.4, family="gaussian", standardize=T, nfolds = 10, type.gaussian="naive"))
## (2) Regularization. Cross validation for hyperparameter tuning. Optimization of model selection to avoid overfitting
par(mfrow = c(1,1)); plot(cv.out)
bestlam <- cv.out$lambda.min; bestlam
#lasso.bestlam <- coef(lasso.mod, s=bestlam)
#bestlam <- cv.out$lambda.1se; bestlam
## Select the best hyperparameter
lasso.pred <- predict(lasso.mod, s=bestlam, newx=foo2[test,], type="nonzero"); str(lasso.pred)
lasso.pred <- predict(lasso.mod, s=bestlam, newx=foo2[test,], type="response"); lasso.pred
lasso.pred <- predict(lasso.mod, s=bestlam, newx=foo2[test,], type="class")
table(lasso.pred, y[test])		## Confusion matrix for classification
mean((lasso.pred - y[test])^2)		## Test set MSE for regression
lasso.coef <- predict(lasso.mod, s=bestlam, type = "coefficients")
str(lasso.coef)
gene1leaf.mRMR <- lasso.coef	## SAVE to .Rdata
## show results
ind.lasso <- lasso.coef@i
loc.lasso <- lasso.coef@Dimnames[[1]]
final.select <- loc.lasso[ind.lasso]
foo <- row.names(means[rownames(means) %in% final.select,])
means[rownames(means) %in% final.select,2]
lasso.select <- x[,foo]
dim(lasso.select)
## extract genes from model selection
mmcDat <- means[rownames(means) %in% final.select,]
#write.table(mmcDat,"mmcDat.txt",quote=F,sep="\t")
## LASSO (from the GLM package)

##############################
# Hierarchical clustering
##############################
require(pvclust)
require(gplots)
myanova = as.matrix(mmcDat[,-c(1,2)])
colnames(myanova) <- c(rep("Healthy", 15), rep("Deficient",7))
colnames(myanova) <- c(rep("L",9), rep("PL",6), rep("L",3),rep("PL",4))
colnames(myanova) <- c(rep("None", 6), rep("Cocktail", 9), rep("Tiso",7))
colnames(myanova) <- c(rep("E",3), rep("T",3), rep("VC",3),rep("PC",3),rep("JC",3),rep("VT",3),rep("PT",3), rep("JT",1))
#rownames(myanova) <- paste(seq(1:nrow(mmcDat)),sep="-",mmcDat[,2])
rownames(myanova) <- rownames(mmcDat)
head(myanova)
mydataanova <- t(scale(t(myanova)))
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R")
hra <- hclust(as.dist(1-cor(t(mydataanova), method="pearson")), method="complete") ## ROWS (genes)
hca <- hclust(as.dist(1-cor(mydataanova, method="pearson")), method="complete")	## COL (samples)
heatmap(myanova, Rowv=as.dendrogram(hra), Colv=as.dendrogram(hca), col=my.colorFct(), scale="row")
## to be continued in Heatmap.R


##############################
# Enrichment analysis
##############################

z=rownames(mmcDat)
geneList <- factor(as.integer(geneNames %in% z))
## extract locus names
## continue from "Build top Go data"


##############################
# Building classifier on subset model
##############################
## Choosing the right Hyper-parameters. GRID ANALYSIS
require(caret)
dat <- data.frame(y=y, lasso.select); dim(dat)
#ctl=expand.grid(.size=seq(1,20,length=40), .decay=10^seq(-1,-5,length=40))	## nnet
#ctl=expand.grid(.mtry=seq(1:15),.coefReg=10^seq(-1,-3,length=40),.coefImp=10^seq(-1,-2,length=40))## Regularized Random forest (RRF)
#ctl=expand.grid(.C=seq(1:20), .sigma=10^seq(-1,-3,length=40))	## svmRadial
#ctl=expand.grid(.ncomp=seq(1:15))	# PCR
#ctl=expand.grid(.mstop=seq(10:1000,length=20),.prune=no)	## glmboost
#ctl=expand.grid(.lambda=10^seq(10,-2,length=100))	## Ridge
set.seed(1445612321)
model.reg(dat,train,test,method="RRF",folds=10,r=5,tune=10)
nnet0 <- modelTune.reg(dat,train,test,method="nnet",folds=10,r=5,tune=10,ctl)	# Tune hyper-parameters
## Regression
model.clas(dat,train,test,method="pls",folds=10,r=5,tune=10)
modelTune.clas(dat,train,test,method="pls",folds=10,r=5,tune=10,ctl)
## Classification
## Hyper-parameters tuning and model optimization

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
RFF_bag500<- baggingTune(dat[train,],dat[test,],m=1.1,ite=500,methods="RRF",tune=10,gridZ=ctl)## For tuning the hyper-parameters
test100 <- bagging(dat[train,],dat[test,],m=1.1,ite=100,methods="ridge",tune=10)## for testing
bagging.clas(dat[train,],dat[test,],m=1.1,ite=100,methods="nnet",tune=10)## For Classification
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
