################################################################################
#    &&&....&&&    % Project: A multi-step machine learning method-Metabee     #
#  &&&&&&..&&&&&&  % Author: Bo Li, Huachun Yin, Jingxin Tao   
#  &&&&&&&&&&&&&&  % Date: Mar. 1st, 2022                                      #
#   &&&&&&&&&&&&   %                                                           #
#     &&&&&&&&     % Environment: R version 3.5.3;                             #
#       &&&&       % Platform: x86_64-pc-linux-gnu (64-bit)                    #
#        &         %                                                           #
################################################################################

### ****************************************************************************
### code chunk number 01: Improved R script for RNA-seq data simulation.
### ****************************************************************************

# Setting the work directory, under the Windows Operation System. 
# setwd("F:/")

### ------------------------------------------------------------------------ ###
### Step-01. Install all R packages required in this step. 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("preprocessCore")) BiocManager::install("preprocessCore")
if (!require("e1071")) install.packages("e1071")
if (!require("coin")) install.packages("coin")
 
library(preprocessCore)
library(e1071)
library(coin)
library(FSA)
library(pROC)
library(DT)
library(tibble)
library(dplyr)
library(readr)
library(ropls)
### ------------------------------------------------------------------------ ###
### Function-01. A set of SVM-RFE functions for feature selection. 
### Note: derived from johncolby/SVM-RFE/msvmRFE.R in Github. 

#++ Function-1.1 svmRFE.wrap

svmRFE.wrap <- function(test.fold, X, ...) {
  # Wrapper to run svmRFE function while omitting a given test fold
  train.data = X[-test.fold, ]
  test.data  = X[test.fold, ]
  
  # Rank the features
  features.ranked = svmRFE(train.data, ...)
  
  return(list(feature.ids=features.ranked, train.data.ids=row.names(train.data), test.data.ids=row.names(test.data)))
}

#++ Function-1.2 svmRFE

svmRFE <- function(X, k = 1, halve.above = 5000) {
  # Feature selection with Multiple SVM Recursive Feature Elimination (RFE) algorithm
  n = ncol(X) - 1
  
  # Scale data up front so it doesn't have to be redone each pass
  cat('Scaling data...')
  X[, -1] = scale(X[, -1])
  cat('Done!\n')
  flush.console()
  
  pb = txtProgressBar(1, n, 1, style=3)
  
  i.surviving = 1:n
  i.ranked    = n
  ranked.list = vector(length=n)
  
  # Recurse through all the features
  while(length(i.surviving) > 0) {
    if(k > 1) {
      # Subsample to obtain multiple weights vectors (i.e. mSVM-RFE)            
      folds = rep(1:k, len=nrow(X))[sample(nrow(X))]
      folds = lapply(1:k, function(x) which(folds == x))
      
      # Obtain weights for each training set
      w = lapply(folds, getWeights, X[, c(1, 1+i.surviving)])
      w = do.call(rbind, w)
      
      # Normalize each weights vector
      w = t(apply(w, 1, function(x) x / sqrt(sum(x^2))))
      
      # Compute ranking criteria
      v    = w * w
      vbar = apply(v, 2, mean)
      vsd  = apply(v, 2, sd)
      c    = vbar / vsd
    } else {
      # Only do 1 pass (i.e. regular SVM-RFE)
      w = getWeights(NULL, X[, c(1, 1+i.surviving)])
      c = w * w
    }
    
    # Rank the features
    ranking = sort(c, index.return=T)$ix
    if(length(i.surviving) == 1) {
      ranking = 1
    }
    
    if(length(i.surviving) > halve.above) {
      # Cut features in half until less than halve.above
      nfeat = length(i.surviving)
      ncut  = round(nfeat / 2)
      n     = nfeat - ncut
      
      cat('Features halved from', nfeat, 'to', n, '\n')
      flush.console()
      
      pb = txtProgressBar(1, n, 1, style=3)
      
    } else ncut = 1
    
    # Update feature list
    ranked.list[i.ranked:(i.ranked-ncut+1)] = i.surviving[ranking[1:ncut]]
    i.ranked    = i.ranked - ncut
    i.surviving = i.surviving[-ranking[1:ncut]]
    
    setTxtProgressBar(pb, n-length(i.surviving))
    flush.console()
  }
  
  close(pb)
  
  return (ranked.list)
}

#++ Function-1.3 getWeights

getWeights <- function(test.fold, X) {
  # Fit a linear SVM model and obtain feature weights
  train.data = X
  if(!is.null(test.fold)) train.data = X[-test.fold, ]
  
  svmModel = svm(train.data[, -1], train.data[, 1], cost=10, cachesize=500,
                 scale=F, type="C-classification", kernel="linear")
  
  t(svmModel$coefs) %*% svmModel$SV
}

#++ Function-1.4 WriteFeatures

WriteFeatures <- function(results, input, save=T, file='features_ranked.txt') {
  # Compile feature rankings across multiple folds
  featureID = sort(apply(sapply(results, function(x) sort(x$feature, index.return=T)$ix), 1, mean), index=T)$ix
  avg.rank  = sort(apply(sapply(results, function(x) sort(x$feature, index.return=T)$ix), 1, mean), index=T)$x
  feature.name = colnames(input[, -1])[featureID]
  features.ranked = data.frame(FeatureName=feature.name, FeatureID=featureID, AvgRank=avg.rank)
  if(save==T) {
    write.table(features.ranked, file=file, quote=F, row.names=F)
  } else {
    features.ranked
  }
}

#++ Function-1.5 FeatSweep.wrap

FeatSweep.wrap <- function(i, results, input) {
  # Wrapper to estimate generalization error across all hold-out folds, for a given number of top features
  svm.list = lapply(results, function(x) tune(svm,
                                              train.x      = input[x$train.data.ids, 1+x$feature.ids[1:i]],
                                              train.y      = input[x$train.data.ids, 1],
                                              validation.x = input[x$test.data.ids, 1+x$feature.ids[1:i]],
                                              validation.y = input[x$test.data.ids, 1],
                                              # Optimize SVM hyperparamters
                                              ranges       = tune(svm,
                                                                  train.x = input[x$train.data.ids, 1+x$feature.ids[1:i]],
                                                                  train.y = input[x$train.data.ids, 1],
                                                                  ranges  = list(gamma=2^(-12:0), cost=2^(-6:6)))$best.par,
                                              tunecontrol  = tune.control(sampling='fix'))$perf)
  
  error = mean(sapply(svm.list, function(x) x$error))
  return(list(svm.list=svm.list, error=error))
}

#++ Function-1.6 PlotErrors

PlotErrors <- function(errors, errors2=NULL, no.info=0.5, ylim=range(c(errors, errors2), na.rm=T), xlab='Number of Features',  ylab='10x CV Error') {
  # Makes a plot of average generalization error vs. number of top features
  AddLine <- function(x, col='black') {
    lines(which(!is.na(errors)), na.omit(x), col=col)
    points(which.min(x), min(x, na.rm=T), col='red')
    text(which.min(x), min(x, na.rm=T), paste(which.min(x), '-', format(min(x, na.rm=T), dig=3)), pos=4, col='red', cex=0.75)
  }
  
  plot(errors, type='n', ylim=ylim, xlab=xlab, ylab=ylab)
  AddLine(errors)
  if(!is.null(errors2)) AddLine(errors2, 'gray30')
  abline(h=no.info, lty=3)
}

### End of Function-01. 
### ------------------------------------------------------------------------ ###


### ------------------------------------------------------------------------ ###
### Function-02. Metabee for metabolites selection.
###  


Metabee <- function(data, 
                    label, 
                    nfold=10,
                    ranktop=100,
                    vipvalue=1,
                    Pvalue=0.05,
                    upFC=1.3,
                    downFC=0.6){
  M_matrix <- data
  
  if(!is.null(nfold)){nfold <- nfold}
  
  if(is.null(nfold)){nfold <- 5}

  ####translate the data
  M_matrix1 <- apply(M_matrix, 1, as.character)
  M_matrix2 <- apply(M_matrix1, 1, as.numeric)
  
  rownames(M_matrix2) <- rownames(M_matrix)
  colnames(M_matrix2) <- colnames(M_matrix)
  
  M_matrix <- M_matrix2
  
  sam.lab<- label #facotr, "Con" is Control groups, and "tet" is experimental groups.
  
  c.p <- which(sam.lab == "Con")
  e.p <- which(sam.lab == "tet")
  
  #########FC
  
  FC <-apply(M_matrix, 1, function(x) {mean(x[e.p])/mean(x[c.p])})
  
  FC <-as.data.frame(FC)
  
  FC$"log2(FC)"<-log2(FC$FC)
  
  FC$FeatureName<-rownames(FC)
  
  ##########log normalization
  
  M_matrix<-log2(M_matrix+0.1)
  
  M_matrix <-as.data.frame(t(M_matrix))
  
  input<- cbind(sam.lab, M_matrix)
  
  input <- as.data.frame(input)
  
  
  ###############OPLSDA
  
  set.seed(1)
  
  oplsda = opls(input[,-1], sam.lab, predI = 1, orthoI = NA, crossvalI=5)
  
  sample.score <- oplsda@scoreMN %>%as.data.frame() %>% mutate(label = sam.lab,
                                                               o1 = oplsda@orthoScoreMN[,1])
  
  vip <- getVipVn(oplsda)
  
  vip <- as.data.frame(vip)
  
  vip$FeatureName<-rownames(vip)
  
  ##########SVM-RFE#####

  nrows <- nrow(input)
  
  set.seed(2)
  
  folds <- rep(1:nfold, len=nrows)[sample(nrows)]
  
  folds <- lapply(1:nfold, function(x) which(folds == x))
  
  set.seed(2)
  
  results <- lapply(folds, svmRFE.wrap, input, k = 10, halve.above = 100)
  
  top.features <- WriteFeatures(results, input, save = FALSE)
  
  #View(top.features)
  
  top.features$FeatureID<-1:nrow(top.features)
  
  deg.svm <- top.features$FeatureName
  
  deg.svm <- as.character(deg.svm)
  
  ##########Permutation#####
  
  all.genes <- colnames(M_matrix)
  
  gene.count <- length(all.genes)
  
  deg.count <- NULL
  
  i <- 0
  
  pb <- txtProgressBar(min = 0, max = gene.count, style = 3, char = "+")
  
  for (g in all.genes) {
    
    i <- i + 1
    
    setTxtProgressBar(pb, i)
    
    # Display the progress bar!
    
    tmp <- Summarize(get(g) ~ sam.lab, data = input, digits = 3)
    
    #print(tmp)
    #boxplot(get(g) ~ sam.lab, data = input)
    
    deg.per <- try(independence_test(get(g) ~ sam.lab, 
                                     data = input,
                                     distribution = approximate(nresample = 10000)), 
                   silent = FALSE)
    
    # deg.Z <- deg.per@statistic@teststatistic
    deg.p <- deg.per@distribution@pvalue(deg.per@statistic@teststatistic)
    
    flush.console()
    
    # This variable, deg.count, stores all the differentially expressed genes.
    
    if (!is.na(deg.p) & deg.p < 1) deg.count <- c(deg.count, g,deg.p) else next
    
  }
  
  close(pb)
  
  
  deg.count1<-matrix(deg.count,ncol=2,byrow = T)
  
  colnames(deg.count1)<-c("FeatureName","Pvalue")
  
  deg.count1<-as.data.frame(deg.count1)
  
  deg.count1$Pvalue<-as.vector(deg.count1$Pvalue)
  
  deg.count1$Pvalue<-as.numeric(deg.count1$Pvalue)
  
  #################
  
  
  if(!is.null(ranktop)){ranktop <- ranktop}
  
  if(is.null(ranktop)){ranktop <- nrow(top.features)}
  
  top.features1<-top.features[1:ranktop,]
  
  top.features1<-merge(FC,top.features1)
  
  top.features1<-merge(top.features1,deg.count1)
  
  sumdata<-merge(top.features1,vip)
  
  sumdata<-sumdata[sumdata$vip>vipvalue & sumdata$Pvalue < Pvalue,]
  #
  sumdata<-sumdata[sumdata$FC>upFC | sumdata$FC<downFC,]
  
  results1 <- list()
  
  results1[[1]] <- sumdata
  results1[[2]] <- deg.svm
  
  names(results1) <- c("result","rank")
  
  return(results1)
 
}


### End of Function-02. 
### ------------------------------------------------------------------------ ###


### ------------------------------------------------------------------------ ###
### Function-03. TOPSIS:Technique for Order Preference by Similarity to an Ideal Solution
###  

z_value <- function(x){
  x / sqrt(sum(x^2))
}


dist <-function(x, std){
  res <- c()
  for ( i in 1 : nrow(x)) {
    res[i] = sqrt(sum((unlist(x[i,-1])-std)^2))
  }
  
  return(res)
}


dat_z <- dataAUC %>% dplyr::mutate(across(c(3:7), z_value))
dat_z <- dat_z[,-1]

z_max <- dat_z %>% summarise(across(c(2:6), max)) %>% unlist
z_min <- dat_z %>% summarise(across(c(2:6), min)) %>% unlist

du <- dist(dat_z, z_max)
dn <- dist(dat_z, z_min)

# calculate CI
dat_z <- dat_z %>% add_column(du = du, dn = dn) %>% 
  mutate(ci= dn/(du+dn)) %>%
  arrange(-ci)


### End of Function-03. 
### ------------------------------------------------------------------------ ###


### ------------------------------------------------------------------------ ###
### Function-04. Other function. 
###  

#####Function-4.1 save pdf

save_pdf <- function(x, filename, width=6, height=4) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#####Function-4.2  Calculating with many distance metrics.
dist_n <- function(x, mtd = "euclidean", p = NULL){
  
  if (!require("philentropy")) install.packages("philentropy")
  
  if (mtd == "maximum") {
    dist(x = x,
         method = mtd,
         diag = FALSE,
         upper = FALSE,
         p = p)
  }else{
    # library(philentropy)
    distance(x = x,
             method = mtd,
             p = p,
             test.na = TRUE,
             unit = "log",
             est.prob = NULL,
             use.row.names = TRUE,
             as.dist.obj = TRUE,
             diag = FALSE,
             upper = FALSE)
  }
}

#####Function-4.3  get annotation for heatmap.
get_anno_for_heatmap2<-function(annocol,annorow=NULL,color=NULL,only.color=F){
  require(plyr)
  require(stringr)
  if(is.null(color)){
    require(RColorBrewer)
    color=c(brewer.pal(12,"Set3"),brewer.pal(8,"Set2"),brewer.pal(9,"Set1"),brewer.pal(8,"Dark2"))
  }
  
  annocolor=do.call(as.list,list(x=annocol))
  annocolor=lapply(annocolor,function(x){if(is.factor(x)){x=levels(x);a=color[1:length(x)];names(a)=x;return(a)}else{x=unique(x);a=color[1:length(x)];names(a)=x;return(a)}})
  if(!is.null(annorow)){
    annocolor.row<-do.call(as.list,list(x=annorow))
    annocolor.row=lapply(annocolor.row,function(x){if(is.factor(x)){x=levels(x);a=color[1:length(x)];names(a)=x;return(a)}else{x=unique(x);a=color[1:length(x)];names(a)=x;return(a)}})
  }else{annocolor.row=NULL}
  annocolor=c(annocolor,annocolor.row)
  annocolor_col<-as.list(annocol)
  annocolor_row<-as.list(annorow)
  annocolor<-c(annocol,annorow)
  annocolor<-lapply(annocolor,function(x){if(is.factor(x)){x=levels(x);return(x)}else{x=unique(x);return(x)}})
  annocolor<-do.call(c,annocolor)
  annocolor<-data.frame(var_name=as.factor(stringr::str_replace(names(annocolor),"[0-9]{1,}$","")),
                        var=annocolor,
                        color=color[1:length(annocolor)])
  annocolor<-split(annocolor,annocolor$var_name)
  annocolor<-lapply(annocolor,function(x){a=x$var;b=as.character(x$color);names(b)=a;return(b)})
  
  
}

### End of Function-04. 
### ------------------------------------------------------------------------ ###


