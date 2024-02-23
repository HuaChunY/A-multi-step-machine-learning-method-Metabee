
source("Chunk01-R script for searching feature.R")

k1 <- 5 #folds for outer loop
k2 <- 5 # folds for inner loop
##################################################################

###################01 Read the data

##################################################################

######Training data
hDEM<-openxlsx::read.xlsx("DEMannotaion.xlsx") ####differential metabolites from Table S4

DMs_Training <- readRDS("DMs_Training.rds")

Proteins_Training<- readRDS("Proteins_Training.rds")

######Testing data

DMs_Test <- readRDS("返修/归一化结果/DMs_Test.rds")

Proteins_Test<- readRDS("返修/归一化结果/Proteins_Test.rds")


###################merge the metabolites and proteins
nor.eset<-rbind(DMs_Training,Proteins_Training)

lable<-colnames(nor.eset)

lable<-substr(lable,1,3)

dim(nor.eset)

nor.eset1 <-log(nor.eset+1)

A_Z_score <- scale(nor.eset1)

mypal1 <- terrain.colors(ncol(A_Z_score))

par(mfrow=c(2,1))

boxplot(nor.eset,col = mypal1,outline=F,names=F,ylab="Abundance")

boxplot(A_Z_score, col = mypal1,outline=F,names=F,ylab="Abundance")

hDEG_order <- hDEM %>% 
  group_by(Ion) %>% 
  arrange(Rank.score) %>% # or  Rank.score
  slice(1:10)

hDEG_order$Comp <- paste0(hDEG_order$FeatureName,hDEG_order$Ion)

dem<-c(hDEG_order$Comp,"β-HCG")

eset.mat <- A_Z_score[dem,]

eset.mat <- as.data.frame(t(eset.mat))

input21 <- cbind(as.factor(lable), eset.mat)

input21 <- as.data.frame(input21)

#saveRDS(input21,"input21.rds")

##############################

test_dt<-rbind(DMs_Test,Proteins_Test)

test_dt1 <-log(test_dt+1)

test_dt1 <- scale(test_dt1)

test21 <- t(test_dt1[dem,])

test21 <- as.data.frame(test21)

#saveRDS(test21,"test21.rds")

##########################



####################################################################################
#######02 the calssification effect of β-HCG and 20 DMs for germinoma and non-iGCT
####################################################################################

########training data
HCG <- t(Proteins_Training)

X <- HCG[,"β-HCG"]

X<-as.data.frame(X)

colnames(X)<-"β-HCG"

y <- substr(rownames(X),1,3)
########test data

X.test<-as.data.frame(Proteins_Test)["β-HCG",]

X.test<-t(X.test)

y.test1<-as.factor(substr(rownames(X.test),1,3))

#############β-HCG: SVM

folds <- k1

rets <- c("threshold", "specificity", "sensitivity", "accuracy")

############repeating 100 for svm model

all.AUC.tables <- lapply(1:100, function(t) {
  
  set.seed(t)
  ############
  
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
  
  pROC::ci(res.roc)
  
  all.ci <- pROC::ci.coords(res.roc, x="best", input = "threshold",best.policy="random",ret=rets)
  
  c(pROC::ci(res.roc),all.ci$threshold,
    all.ci$specificity,
    all.ci$sensitivity,
    all.ci$accuracy)
  
})

full.AUC.HCG <- do.call(rbind, all.AUC.tables)

colnames(full.AUC.HCG) <- paste0(c(rep("AUC",3),
                                   rep("threshold",3),
                                   rep("specificity",3),
                                   rep("sensitivity",3),
                                   rep("accuracy",3)),
                                 rep(c(".low",".median",".high"),5))


HCG_AUC <-  rbind(full.AUC.HCG,apply(full.AUC.HCG,2,mean))

rownames(HCG_AUC) <- c(paste0(rep("seed",100),1:100),"Mean")

#############################################

########training data

Dms20 <- DMs_Training[hDEG_order$Comp,]

dim(Dms20)

input <- t(Dms20)# training data

y <- substr(rownames(input),1,3)

########test data

X.test20<-as.data.frame(DMs_Test)[hDEG_order$Comp,]

X.test20<-t(X.test20)

y.test1<-as.factor(substr(rownames(X.test20),1,3))

#############20 DMs: SVM

rets <- c("threshold", "specificity", "sensitivity", "accuracy")

for (GM in 1:20) {
  
  deg.cer<-colnames(input)[GM]
 
  X <- input[,deg.cer]
  
  X<-as.data.frame(X)
  
  rownames(X)<-rownames(input)
  
  colnames(X)<-deg.cer

  #####test data

  X.test<-as.data.frame(X.test20[,deg.cer])
  rownames(X.test) <- rownames(X.test20) 
  colnames(X.test)<-deg.cer
  
  ##############
  folds <- k1
  all.AUC.tables <- lapply(1:100, function(t) {
    
    
      set.seed(t)
      ############
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
    
    all.ci <- pROC::ci.coords(res.roc, x="best", input = "threshold",best.policy="random",ret=rets)
    
    c(pROC::ci(res.roc),all.ci$threshold,
      all.ci$specificity,
      all.ci$sensitivity,
      all.ci$accuracy)
    
  })
  
  full.AUC <- do.call(rbind, all.AUC.tables)
  
  colnames(full.AUC) <- paste0(c(rep("AUC",3),
                                 rep("threshold",3),
                                 rep("specificity",3),
                                 rep("sensitivity",3),
                                 rep("accuracy",3)),
                               rep(c(".low",".median",".high"),5))
  
  
  full.AUC <-  rbind(full.AUC,apply(full.AUC,2,mean))
  
  rownames(full.AUC) <- c(paste0(rep("seed",100),1:100),"Mean")
  
  assign(paste0("AUC",deg.cer),full.AUC)
  
}

AUC21 <- list(get(paste0("AUC",dem[1])),
              get(paste0("AUC",dem[2])),
              get(paste0("AUC",dem[3])),get(paste0("AUC",dem[4])),
              get(paste0("AUC",dem[5])),get(paste0("AUC",dem[6])),
              get(paste0("AUC",dem[7])),get(paste0("AUC",dem[8])),
              get(paste0("AUC",dem[9])),get(paste0("AUC",dem[10])),
              get(paste0("AUC",dem[11])),get(paste0("AUC",dem[12])),
              get(paste0("AUC",dem[13])),get(paste0("AUC",dem[14])),
              get(paste0("AUC",dem[15])),get(paste0("AUC",dem[16])),
              get(paste0("AUC",dem[17])),get(paste0("AUC",dem[18])),
              get(paste0("AUC",dem[19])),get(paste0("AUC",dem[20])),
              HCG_AUC)

names(AUC21) <- c(paste0(dem[1], "AUC"),paste0(dem[2], "AUC"),
                  paste0(dem[3], "AUC"),paste0(dem[4], "AUC"),
                  paste0(dem[5], "AUC"),paste0(dem[6], "AUC"),
                  paste0(dem[7], "AUC"),paste0(dem[8], "AUC"),
                  paste0(dem[9], "AUC"),paste0(dem[10], "AUC"),
                  paste0(dem[11],"AUC"),paste0(dem[12], "AUC"),
                  paste0(dem[13],"AUC"),paste0(dem[14], "AUC"),
                  paste0(dem[15],"AUC"),paste0(dem[16], "AUC"),
                  paste0(dem[17],"AUC"),paste0(dem[18], "AUC"),
                  paste0(dem[19],"AUC"),paste0(dem[20], "AUC"),
                  "β-HCG+AUC")
##############

############# 21 biomarkers: SVM 

rets <- c("threshold", "specificity", "sensitivity", "accuracy")

folds <- k1

all.AUC.tables <- lapply(1:100, function(t) {
  
  set.seed(t)
  ############
  
  test.fold <- split(sample(1:length(y)), 1:folds) 
  
  set.seed(10)
  
  all.pred.tables <- lapply(1:folds, function(si) {
    
    X <- input21[,-1]
    
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

full.AUC <- do.call(rbind, all.AUC.tables)

colnames(full.AUC) <- paste0(c(rep("AUC",3),
                               rep("threshold",3),
                               rep("specificity",3),
                               rep("sensitivity",3),
                               rep("accuracy",3)),
                             rep(c(".low",".median",".high"),5))

inte21AUC <-  rbind(full.AUC,apply(full.AUC,2,mean))

rownames(inte21AUC) <- c(paste0(rep("seed",100),1:100),"Mean")

AUC22 <- AUC21
AUC22[[22]] <- inte21AUC
names(AUC22)[22] <- "21 biomarkers"

#########################

library(ggplot2)

data.list <- list()

for (i in names(AUC22)) {
  temp <- cbind(rep(i,100),AUC22[[i]][1:100,"AUC.median"])
  temp1 <- cbind(temp,AUC22[[i]][1:100,"specificity.median"])
  temp2 <- cbind(temp1,AUC22[[i]][1:100,"sensitivity.median"])
  temp3 <- cbind(temp2,AUC22[[i]][1:100,"accuracy.median"])
  colnames(temp3) <- c("name","AUC.median","specificity.median","sensitivity.median","accuracy.median")
  data.list[[i]] <- temp3
}


data.merge <- as.data.frame(matrix(nrow = 0,ncol = 5))
for (i in names(data.list)) {
  data.merge <- rbind(data.merge,data.list[[i]])
}

data.merge$name <- as.character(data.merge$name)
data.merge$name[1:100] <- rep("NAAG-",100)
data.merge$name[101:200] <- rep("2PY-",100)
data.merge$name[1001:1100] <- rep("2PY+",100)
data.merge$name[1201:1300] <- rep("NAAG-",100)



data.merge$AUC.median <- as.numeric(as.character(data.merge$AUC.median))
data.merge$specificity.median <- as.numeric(as.character(data.merge$specificity.median))
data.merge$sensitivity.median <- as.numeric(as.character(data.merge$sensitivity.median))
data.merge$accuracy.median <- as.numeric(as.character(data.merge$accuracy.median))

data.merge$name <- gsub("AUC","",data.merge$name)

p1 <- ggplot(data=data.merge, aes(x= reorder(name, AUC.median,decreasing = T),y= AUC.median))+
  geom_boxplot(color="#336699", size=0.5, width=0.3,outlier.size = 0.1)+
  scale_color_manual("#336699")+
  #facet_wrap(~ Name, scales="free",nrow =5)+
  theme_bw()+
  ylab("AUC")+
  theme(text = element_text(family="serif",size = 10), #"Times New Roman","serif"
        panel.grid = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,angle = 60, hjust = 1),
        axis.text.y = element_text(size=6),
        legend.position ="top",
        legend.text=element_text(size=6),
        legend.title = element_text(size=6)
  )

p2 <- ggplot(data=data.merge, aes(x= reorder(name, specificity.median,decreasing = F),y= specificity.median))+
  geom_boxplot(color="#336699", size=0.5,width=0.3, outlier.size = 0.1)+
  scale_color_manual("#336699")+
  #facet_wrap(~ Name, scales="free",nrow =5)+
  theme_bw()+
  ylab("specificity")+
  theme(text = element_text(family="serif",size = 10), #"Times New Roman","serif"
        panel.grid = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,angle = 60, hjust = 1),
        axis.text.y = element_text(size=6),
        legend.position ="top",
        legend.text=element_text(size=6),
        legend.title = element_text(size=6)
  )


p3 <- ggplot(data=data.merge, aes(x= reorder(name, sensitivity.median,decreasing = F),y= sensitivity.median))+
  geom_boxplot(color="#336699", size=0.5, width=0.3,outlier.size = 0.1)+
  scale_color_manual("#336699")+
  #facet_wrap(~ Name, scales="free",nrow =5)+
  theme_bw()+
  ylab("sensitivity")+
  theme(text = element_text(family="serif",size = 10), #"Times New Roman","serif"
        panel.grid = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,angle = 60, hjust = 1),
        axis.text.y = element_text(size=6),
        legend.position ="top",
        legend.text=element_text(size=6),
        legend.title = element_text(size=6)
  )


p4 <- ggplot(data=data.merge, aes(x= reorder(name, accuracy.median,decreasing = F),y= accuracy.median))+
  geom_boxplot(color="#336699", width=0.3,size=0.5, outlier.size = 0.1)+
  scale_color_manual("#336699")+
  #facet_wrap(~ Name, scales="free",nrow =5)+
  theme_bw()+
  ylab("accuracy")+
  theme(text = element_text(family="serif",size = 10), #"Times New Roman","serif"
        panel.grid = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,angle = 60, hjust = 1),
        axis.text.y = element_text(size=6),
        legend.position ="top",
        legend.text=element_text(size=6),
        legend.title = element_text(size=6)
  )

ggpubr::ggarrange(p1,p2,p3,p4, ncol=1,nrow = 4)



