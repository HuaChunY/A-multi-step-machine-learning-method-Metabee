library(caret)
library(pROC)

source("Chunk01-R script for searching feature.R")

k1 <- 5 #folds for outer loop
k2 <- 5 # folds for inner loop

#####################read the data

input21 <- readRDS("input21.rds")

y <- input21$`as.factor(lable)`

deg.cer<-dat_z$names[2]

X <- input21[,deg.cer]

colnames(X) <- paste0("M",1:ncol(X)) #re-name for randomforest method

test21 <- readRDS("test21.rds")

X.test <- test21[,deg.cer]

colnames(X.test) <- paste0("M",1:ncol(X.test))#re-name for randomforest method

y.test1 <- rownames(test21)

##############################
folds <- k1

rets <- c("threshold", "specificity", "sensitivity", "accuracy")

####################################

##################01 SVM

####################################

AUC.100.tables <- lapply(1:100, function(t) {
  
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
    model <- svm(X.train, y.train, kernel = "radial", 
                 gamma = 0.1,cachesize=10,
                 type="C-classification", 
                 prob = TRUE, cross = k2)
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
AUC.100.mean <- do.call(rbind, AUC.100.tables)

colnames(AUC.100.mean) <- paste0(c(rep("AUC",3),
                                  rep("threshold",3),
                                  rep("specificity",3),
                                  rep("sensitivity",3),
                                  rep("accuracy",3)),
                                rep(c(".low",".median",".high"),5))



SVM_AUC <-  rbind(AUC.100.mean,apply(AUC.100.mean,2,mean))

rownames(SVM_AUC) <- c(paste0(rep("seed",100),1:100),"Mean")

View(SVM_AUC)

####################################

##################02 RF

####################################
AUC.100.tables <- lapply(1:100, function(t) {
  
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
    
    model <- randomForest::randomForest(y.train ~ ., X.train,
                                        proximity = TRUE,)
    
    predict.test <- predict(model, X1.test, type = "prob")
    
    prob.benign <- data.frame(predict.test)
    
    data.frame(y.test =c(as.character(y.test1),as.character(y[test.id])),  y.pred = prob.benign[,2]) # returning this
    
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
AUC.100.mean <- do.call(rbind, AUC.100.tables)

colnames(AUC.100.mean) <- paste0(c(rep("AUC",3),
                                  rep("threshold",3),
                                  rep("specificity",3),
                                  rep("sensitivity",3),
                                  rep("accuracy",3)),
                                rep(c(".low",".median",".high"),5))



RF_AUC <-  rbind(AUC.100.mean,apply(AUC.100.mean,2,mean))

rownames(RF_AUC) <- c(paste0(rep("seed",100),1:100),"Mean")

##################

####################################

##################03 NB

####################################

AUC.100.tables <- lapply(1:100, function(t) {
  
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
    
    model = caret::train(X.train,y.train,'nb',trControl=trainControl(method='cv',number=k2))
    
    predict.test <- predict(model$finalModel, X1.test, type = "prob")
    
    prob.benign <- data.frame(predict.test)
    
    data.frame(y.test =c(as.character(y.test1),as.character(y[test.id])), y.pred = prob.benign[,2]) # returning this
    
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
AUC.100.mean <- do.call(rbind, AUC.100.tables)

colnames(AUC.100.mean) <- paste0(c(rep("AUC",3),
                                  rep("threshold",3),
                                  rep("specificity",3),
                                  rep("sensitivity",3),
                                  rep("accuracy",3)),
                                rep(c(".low",".median",".high"),5))



NB_AUC <-  rbind(AUC.100.mean,apply(AUC.100.mean,2,mean))

rownames(NB_AUC) <- c(paste0(rep("seed",100),1:100),"Mean")
View(NB_AUC)
##################

####################################

##################04 KNN

####################################
AUC.100.tables <- lapply(1:100, function(t) {
  
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
    X.train$y.train <- y.train
    set.seed(1)
    
    model = train(y.train ~ ., 
                  X.train, 'knn',
                  trControl=trainControl(method='cv',number=k2))
    
    predict.test <- predict(model, X1.test, type = "prob")
    
    prob.benign <- data.frame(predict.test)
    
    data.frame(y.test =c(as.character(y.test1),as.character(y[test.id])), y.pred = prob.benign[,2]) # returning this
    
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
AUC.100.mean <- do.call(rbind, AUC.100.tables)

colnames(AUC.100.mean) <- paste0(c(rep("AUC",3),
                                  rep("threshold",3),
                                  rep("specificity",3),
                                  rep("sensitivity",3),
                                  rep("accuracy",3)),
                                rep(c(".low",".median",".high"),5))

KNN_AUC <-  rbind(AUC.100.mean,apply(AUC.100.mean,2,mean))

rownames(KNN_AUC) <- c(paste0(rep("seed",100),1:100),"Mean")
View(KNN_AUC)

AUClist <- list(SVM_AUC, RF_AUC,NB_AUC,KNN_AUC)

names(AUClist) <- c("SVM","Random forest","Naive Bayes","K-NN")

meanAUC <- rbind(SVM_AUC[101,],RF_AUC[101,],NB_AUC[101,],KNN_AUC[101,])

meanAUC <- as.data.frame(meanAUC)

meanAUC$name <- c("SVM","Random forest","Naive Bayes","K-NN")

meanAUC <- meanAUC[,c("name","AUC.median","specificity.median","sensitivity.median","accuracy.median")]
######################

AUClist <- readRDS("返修/二审返修/13test/results of four classification algorithms.rds")

data.list <- list()

for (i in names(AUClist)) {
  temp <- cbind(rep(i,100),AUClist[[i]][1:100,"AUC.median"])
  temp1 <- cbind(temp,AUClist[[i]][1:100,"specificity.median"])
  temp2 <- cbind(temp1,AUClist[[i]][1:100,"sensitivity.median"])
  temp3 <- cbind(temp2,AUClist[[i]][1:100,"accuracy.median"])
  colnames(temp3) <- c("name","AUC.median","specificity.median","sensitivity.median","accuracy.median")
  data.list[[i]] <- temp3
}


data.merge <- as.data.frame(matrix(nrow = 0,ncol = 5))
for (i in names(data.list)) {
  data.merge <- rbind(data.merge,data.list[[i]])
}


data.merge$AUC.median <- as.numeric(as.character(data.merge$AUC.median))
data.merge$specificity.median <- as.numeric(as.character(data.merge$specificity.median))
data.merge$sensitivity.median <- as.numeric(as.character(data.merge$sensitivity.median))
data.merge$accuracy.median <- as.numeric(as.character(data.merge$accuracy.median))

data.merge1 <- reshape2::melt(data.merge,
                             id.vars="name",
                             variable.name = "Group", 
                             value.name = "value")
  
  ggplot(data=data.merge1, aes(x= name,y= value))+
  geom_boxplot(aes(color=name), size=0.5, width=0.3,outlier.size = 0.1)+
  scale_color_manual(values = c("#336699","#CBBB7F","#D15C62","#669999"))+
  facet_wrap(~ Group, scales="free",ncol =4)+
  theme_bw()+
  ylab("Value")+
  theme(text = element_text(family="Times New Roman",size = 8), #"Times New Roman","serif"
        panel.grid = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,angle = 60, hjust = 1),
        axis.text.y = element_text(size=8),
        #legend.position ="top",
        legend.text=element_text(size=8),
        legend.title = element_text(size=8)
  )

######################################################################## 
  
##################05 Survival analysis
  
######################################################################## 
  
  DMs_Training <- readRDS("返修/归一化结果/DMs_Training.rds")
  DMs_Training <- log(DMs_Training+1)
  sur_data_DM<-DMs_Training[,which(substr(colnames(DMs_Training),1,3)%in%"tet")]
  
  Proteins_Training <- readRDS("返修/归一化结果/Proteins_Training.rds")
  Proteins_Training <- log(Proteins_Training+1)
  sur_data_P<-Proteins_Training[,which(substr(colnames(Proteins_Training),1,3)%in%"tet")]
  
  sur_time <- openxlsx::read.xlsx("survival.xlsx") #the survival status of patients with germinoma  shown in Tabel S1
  sur_time[which(sur_time$status %in%0),1] <- 2
  sur_time[which(sur_time$status %in%1),1] <- 0
  sur_time[which(sur_time$status %in%2),1] <- 1
  #View(sur_data)
  
  #summary(sur_data)
  library(dplyr)
  colnames(sur_data_DM)
  colnames(sur_data_P)
  sur_data <- rbind(sur_data_DM,sur_data_P)
  sur_data <- as.data.frame(t(sur_data))
  sur_data <- sur_data[,c("Inosine-","N-Acetylaspartylglutamate (NAAG)+","β-HCG")]
  sur_data1<-sur_data %>% 
    mutate(NAAG=case_when(`N-Acetylaspartylglutamate (NAAG)+` < quantile(sur_data$`N-Acetylaspartylglutamate (NAAG)+`)[3] ~"L",
                          `N-Acetylaspartylglutamate (NAAG)+` > quantile(sur_data$`N-Acetylaspartylglutamate (NAAG)+`)[3] ~"H",
                          TRUE ~"M"),
           Inosine=case_when(
             `Inosine-` < quantile(sur_data$`Inosine-`)[3] ~"L",
             `Inosine-` > quantile(sur_data$`Inosine-`)[3] ~"H",
             TRUE ~"M"),
           
           HCGG=case_when(
             `β-HCG` < quantile(sur_data$`β-HCG`)[3] ~"L",
             `β-HCG` > quantile(sur_data$`β-HCG`)[3]~"H",
             TRUE ~"M")
    )
  
  sur_data2 <- sur_data1[,c("NAAG","Inosine","HCGG")]
  sur_data2$sample <- rownames(sur_data2)
  
  sur_time2 <- sur_time[which(sur_time$sample%in%rownames(sur_data2)),]
  
  sur_data3 <- merge(sur_time2,sur_data2,by="sample")

  library(survival)
  library(survminer) 
  
  fit <- survfit(Surv(time, status) ~ NAAG, data = sur_data3)
  
  surv_diff <- survdiff(Surv(time, status) ~ NAAG, data = sur_data3)
  
  #p.val <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  
  NAAG <- ggsurvplot(
    fit,                     # survfit object with calculated statistics.
    pval = TRUE,             # show p-value of log-rank test.
    conf.int = FALSE,         # show confidence intervals for 
    pval.method = TRUE,
    log.rank.weights = "survdiff", #"survdiff", "1", "n", "sqrtN", "S1","S2"
    conf.int.style = "step",  # customize style of confidence intervals
    xlab = "Time in days",   # customize X axis label.
    break.time.by = 200,     # break X axis in time intervals by 200.
    #ggtheme = theme_light(), # customize plot and risk table with a theme.
    ggtheme = theme(text = element_text(family="Times New Roman",size = 6),
                    panel.background=element_blank()),
    risk.table = "absolute",  # absolute number and percentage at risk.
    risk.table.y.text.col = T,# colour risk table text annotations.
    risk.table.y.text = FALSE,# show bars instead of names in text annotations
    # in legend of risk table.
    ncensor.plot = FALSE,      # plot the number of censored subjects at time t
    surv.median.line = "hv",  # add the median survival pointer.
    legend.labs = c("High", "Low"),    # change legend labels.
    palette = c("#E7B800", "#2E9FDF") # custom color palettes.
  )
  