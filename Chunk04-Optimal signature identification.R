#setwd("/Users/")

source("Chunk01-R script for searching feature.R")

k1 <- 5 #folds for outer loop
k2 <- 5 # folds for inner loop

#####################read the data

input21 <- readRDS("input21.rds")

y <- input21$`as.factor(lable)`

input21 <- input21[,-1]

dem <- colnames(input21)

test21 <- readRDS("test21.rds")

y.test1 <- rownames(test21)

#################################################################################

###########################01 combination searching

#################################################################################


rets <- c("threshold", "specificity", "sensitivity", "accuracy")

all.cmatAUC.tables <- lapply(2:7, function(tt) {
  
  cmat <- combinations(21,tt)
  
  cmat.tables <- lapply(1:nrow(cmat), function(ii) {
    
    dem.pos<-cmat[ii,]
    
    dem.ser <- dem[dem.pos]
    
    X <- as.data.frame(input21[,dem.ser]) 
    
    colnames(X) <- dem.ser
    rownames(X) <- rownames(input21)
    
    X.test<-test21[,dem.ser]
    
    colnames(X.test) <- dem.ser
    ##############
    folds <- k1
    
    ############
    
    AUC.10.tables <- lapply(1:100, function(t) {
      
      set.seed(t)
      
      test.fold <- split(sample(1:length(y)), 1:folds) 
      
      set.seed(10)
      
      all.pred.tables <- lapply(1:folds, function(si) {
        
        test.id <- test.fold[[si]]
        
        X.train <- as.data.frame(X[-test.id, ])
        
        rownames(X.train)<-rownames(X)[-test.id]
        
        colnames(X.train)<-colnames(X)
        
        y.train <- as.factor(y[-test.id])
        
        X.t<-as.data.frame(X[test.id, ])
        
        colnames(X.t)<-colnames(X)
        
        rownames(X.t)<-rownames(X)[test.id]
        
        X1.test <-rbind(X.test,X.t)
        
        set.seed(1)
        model <- svm(X.train, y.train, kernel = "radial", prob = TRUE, cross = k2) # some tuning may be needed
        
        predict.test <- predict(model, X1.test, prob = TRUE)
        
        prob.benign <- attr(predict.test, "probabilities")[, 2]
        
        data.frame(y.test =c(as.character(y.test1),as.character(y[test.id])), y.pred = prob.benign) # returning this
        
      })
      
      full.pred.table <- do.call(rbind, all.pred.tables)
      
      res.roc <- roc(full.pred.table$y.test, 
                     full.pred.table$y.pred, 
                     plot = FALSE)
      
      ##########
      
      pROC::ci(res.roc)
      
      all.ci <- pROC::ci.coords(res.roc, x="best", input = "threshold",best.policy="random",ret=rets)
      
      c(pROC::ci(res.roc),all.ci$threshold,
        all.ci$specificity,
        all.ci$sensitivity,
        all.ci$accuracy)
      
    })
    AUC.10.mean <- do.call(rbind, AUC.10.tables)
    
    colnames(AUC.10.mean) <- paste0(c(rep("AUC",3),
                                      rep("threshold",3),
                                      rep("specificity",3),
                                      rep("sensitivity",3),
                                      rep("accuracy",3)),
                                    rep(c(".low",".median",".high"),5))
    
    name.mat <- paste0(dem.ser,collapse="&")
    
    AUC.100.mean <- as.data.frame(AUC.100.mean)
    
    AUC.mean <- apply(AUC.100.mean,2,mean)
    AUC.mean <- as.data.frame(t(AUC.mean))
    AUC.mean$names <- name.mat
    AUC.mean
  })
  
  cmat.text <- do.call(rbind, cmat.tables)
  
  data.frame(cmat.text)
  
  
})

all.cmatAUC <- do.call(rbind, all.cmatAUC.tables)

all.cmatAUC <- as.data.frame(all.cmatAUC)

all.cmatAUC1 <- all.cmatAUC[all.cmatAUC$AUC.median>0.80,]
all.cmatAUC1 <- all.cmatAUC1[all.cmatAUC1$specificity.median>0.80,]
all.cmatAUC1 <- all.cmatAUC1[all.cmatAUC1$sensitivity.median>0.80,]
all.cmatAUC1 <- all.cmatAUC1[all.cmatAUC1$accuracy.median>0.80,]

all.cmatAUC1 <- all.cmatAUC1[order(all.cmatAUC1$accuracy.median,decreasing = T),]

all.cmatAUC1 <- all.cmatAUC1[1:20,]


#################################################################################

###########################02 TOPSIS

#################################################################################
#######TOPSIS:Technique for Order Preference by Similarity to an Ideal Solution
dataAUC <- all.cmatAUC1

namelist <- strsplit(dataAUC$names,"&")


for (ii in 1:nrow(dataAUC)) {
  
  dem.ser <-namelist[[ii]]
  
  n <- length(dem.ser)
  
  numscore <- 1/n
  
  dataAUC$numscore[ii] <- numscore
}

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

rownames(dataAUC) <- dataAUC$names

dataAUC <- dataAUC[,c("AUC.median",
                      "specificity.median",
                      "sensitivity.median",
                      "accuracy.median","numscore")]

dat_z <- dataAUC %>% dplyr::mutate(across(c("AUC.median",
                                            "specificity.median",
                                            "sensitivity.median",
                                            "accuracy.median","numscore"), z_value))


z_max <- dat_z %>% summarise(across(c("AUC.median",
                                      "specificity.median",
                                      "sensitivity.median",
                                      "accuracy.median","numscore"), max)) %>% unlist
z_min <- dat_z %>% summarise(across(c("AUC.median",
                                      "specificity.median",
                                      "sensitivity.median",
                                      "accuracy.median","numscore"), min)) %>% unlist


du <- dist(dat_z, z_max)
dn <- dist(dat_z, z_min)


dat_z <- dat_z %>% add_column(du = du, dn = dn) %>% 
  mutate(ci= dn/(du+dn)) %>%
  arrange(-ci)

#################################################################################

###########################03 survival analysis

#################################################################################



deg.cer<-dat_z$names[2]

