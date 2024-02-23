################################################################################
#    &&&....&&&    % 
#  &&&&&&..&&&&&&  % Author: Bo Li, Huachun Yin, Jingxin Tao   
#  &&&&&&&&&&&&&&  % Date: Mar. 1st, 2024                                      #
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

if (!require("edgeR")) BiocManager::install("edgeR")
if (!require("preprocessCore")) BiocManager::install("preprocessCore")
if (!require("meta")) install.packages("meta")
if (!require("e1071")) install.packages("e1071")
if (!require("coin")) install.packages("coin")
if (!require("ggsci")) install.packages("ggsci")

library(preprocessCore)
library(compcodeR)
library(meta)
library(e1071)
library(coin)
library(venn)
library(ggsci)
library(preprocessCore)
library(FSA)
library(pROC)
library(simpleaffy)
library(RankProd)
library(bspec)
library(netClass)
library(tidyverse)
library(FSinR)
library(stringr)

### ------------------------------------------------------------------------ ###
### Function-01. A set of SVM-RFE functions for gene selection. 
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
  n <- ncol(X) - 1
  
  # Scale data up front so it doesn't have to be redone each pass
  cat('Scaling data...')
  X[, -1] <- scale(X[, -1])
  cat('Done!\n')
  flush.console()
  
  pb <- txtProgressBar(1, n, 1, style=3)
  
  i.surviving <- 1:n
  i.ranked <- n
  ranked.list <- vector(length=n)
  
  # Recurse through all the features
  while(length(i.surviving) > 0) {
    if(k > 1) {
      set.seed(78)
      # Subsample to obtain multiple weights vectors (i.e. mSVM-RFE)            
      folds <- rep(1:k, len=nrow(X))[sample(nrow(X))]
      folds <- lapply(1:k, function(x) which(folds == x))
      
      # Obtain weights for each training set
      w <- lapply(folds, getWeights, X[, c(1, 1+i.surviving)])
      w <- do.call(rbind, w)
      
      # Normalize each weights vector
      w <- t(apply(w, 1, function(x) x / sqrt(sum(x^2))))
      
      # Compute ranking criteria
      v <- w * w
      vbar <- apply(v, 2, mean)
      vsd <- apply(v, 2, sd)
      c <- vbar / vsd
    } else {
      # Only do 1 pass (i.e. regular SVM-RFE)
      w <- getWeights(NULL, X[, c(1, 1+i.surviving)])
      c <- w * w
    }
    
    # Rank the features
    ranking <- sort(c, index.return=T)$ix
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
    i.ranked <- i.ranked - ncut
    i.surviving <- i.surviving[-ranking[1:ncut]]
    
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


