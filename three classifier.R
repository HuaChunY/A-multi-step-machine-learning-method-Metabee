setwd("/Users/huachunyin/desktop/Xinqiao/tumor/NO/")

library(caret)
library(pROC)
test.mat <- openxlsx::read.xlsx("test.xlsx",rowNames=TRUE)

############### 
#eset.mat is original data

label <- substr(colnames(test.mat),1,3) #set the label 


###############

eset.mat <- openxlsx::read.xlsx("eset.mat.xlsx",rowNames=TRUE)

sam.lab <- substr(colnames(eset.mat),1,3) #group 

sam.lab <- as.factor(sam.lab)

input <- cbind(sam.lab, as.data.frame(t(eset.mat)))

input <- as.data.frame(input)

#View(input)

##################


deg.cer<-c("N-Acetylaspartylglutamate (NAAG)-","Dodecanoic acid-",
           
           "β-HCG")

sam.iris <- input[, c("sam.lab", deg.cer)]


head(sam.iris)

sam.iris<-as.data.frame(sam.iris)
#sam.iris2<-sam.iris[,2]

X <- sam.iris[, names(sam.iris) != "sam.lab"]

X<-as.data.frame(sam.iris[,-1])

rownames(X)<-rownames(sam.iris)

colnames(X)<-c("NAAG","DAA","hCG")

y <- as.character(sam.iris$sam.lab)

#####test data

y.test1<- label
X.test1<- t(test.mat[deg.cer,])
colnames(X.test1)<-c("NAAG","DAA","hCG")

folds <- 10
set.seed(10)

test.fold <- split(sample(1:length(y)), 1:folds) 

set.seed(10)

##################RF

all.pred.tables <- lapply(1:folds, function(si) {
  
  test.id <- test.fold[[si]]
  
  X.train <- as.data.frame(X[-test.id, ])
  
  rownames(X.train)<-rownames(X)[-test.id]
  
  colnames(X.train)<-colnames(X)
  
  y.train <- as.factor(y[-test.id])
  
  X.t<-as.data.frame(X[test.id, ])
  
  colnames(X.t)<-colnames(X)
  
  rownames(X.t)<-rownames(X)[test.id]
  
  X1.test <-rbind(X.test1,X.t)
  
  model <- randomForest::randomForest(y.train ~ ., X.train,
                                      proximity = TRUE,)
  
  predict.test <- predict(model, X1.test, type = "prob")
  
  prob.benign <- data.frame(predict.test)
  
  
  data.frame(y.test =c(as.character(y.test1),y[test.id]), y.pred = prob.benign[,2]) # returning this
  
  
})

full.pred.table <- do.call(rbind, all.pred.tables)


res.roc.RF <- roc(full.pred.table$y.test, 
               full.pred.table$y.pred, 
               plot = TRUE, 
               xlab = "False Positive Percentage", 
               ylab = "TRUE Positive Percentage", 
               legacy.axes = TRUE)
rets <- c("threshold", "specificity", "sensitivity", "accuracy")
pROC::ci.coords(res.roc.RF, x="best", input = "threshold",best.policy="random",ret=rets)

##################

set.seed(10)

##################NB

all.pred.tables <- lapply(1:folds, function(si) {
  
  test.id <- test.fold[[si]]
  
  X.train <- as.data.frame(X[-test.id, ])
  
  rownames(X.train)<-rownames(X)[-test.id]
  
  colnames(X.train)<-colnames(X)
  
  y.train <- as.factor(y[-test.id])
  
  X.t<-as.data.frame(X[test.id, ])
  
  colnames(X.t)<-colnames(X)
  
  rownames(X.t)<-rownames(X)[test.id]
  
  X1.test <-rbind(X.test1,X.t)
  
  
  model = train(X.train,y.train,'nb',trControl=trainControl(method='cv',number=10))
  
  
  predict.test <- predict(model$finalModel, X1.test, type = "prob")
  
  prob.benign <- data.frame(predict.test)
  
  data.frame(y.test =c(as.character(y.test1),y[test.id]), y.pred = prob.benign[,2]) # returning this
  
})

full.pred.table <- do.call(rbind, all.pred.tables)


res.roc.NB <- roc(full.pred.table$y.test, 
               full.pred.table$y.pred, 
               plot = TRUE, 
               xlab = "False Positive Percentage", 
               ylab = "TRUE Positive Percentage", 
               legacy.axes = TRUE)
pROC::ci.coords(res.roc.NB, x="best", input = "threshold",best.policy="random",ret=rets)

##################

set.seed(10)

##################knn

all.pred.tables <- lapply(1:folds, function(si) {
  
  test.id <- test.fold[[si]]
  
  X.train <- as.data.frame(X[-test.id, ])
  
  rownames(X.train)<-rownames(X)[-test.id]
  
  colnames(X.train)<-colnames(X)
  
  y.train <- as.factor(y[-test.id])
  
  X.t<-as.data.frame(X[test.id, ])
  
  colnames(X.t)<-colnames(X)
  
  rownames(X.t)<-rownames(X)[test.id]
  
  X1.test <-rbind(X.test1,X.t)
 
   X.train$y.train <- y.train
  
  model = train(y.train ~ ., 
                X.train, 'knn',
                trControl=trainControl(method='cv',number=10))
 
  predict.test <- predict(model, X1.test,type="prob")
  
  prob.benign <- data.frame(predict.test)
  
  data.frame(y.test =c(as.character(y.test1),y[test.id]), y.pred = prob.benign[,2]) # returning this
  
})

full.pred.table <- do.call(rbind, all.pred.tables)


res.roc.KNN <- roc(full.pred.table$y.test, 
               full.pred.table$y.pred, 
               plot = TRUE, 
               xlab = "False Positive Percentage", 
               ylab = "TRUE Positive Percentage", 
               legacy.axes = TRUE)

res.roc.KNN
pROC::ci.coords(res.roc.KNN, x="best", input = "threshold",best.policy="random",ret=rets)

##################

##################

ggroc(list("Random Forest"=res.roc.RF,
           "Naive Bayes"=res.roc.NB,
           "k-Nearest Neighbors"=res.roc.KNN),
      alpha=0.5,#colour=c("#377EB8","#A254CC","#CCC012","#F06139"),
      pe=2,size=1,legacy.axes = TRUE)+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               color="grey", linetype="dashed")+
  #ggtitle(deg.cer)+
  theme(text = element_text(family="Times New Roman"),
        plot.title = element_text(family = "Times New Roman", 
                                  size = 8, 
                                  hjust = 0.5, 
                                  angle = 0, 
        ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=8),
        axis.text = element_text(size=6))+
  xlab("False Positive Percentage")+
  ylab("TRUE Positive Percentage")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0.005,0))
  #annotate("text",x=0.25,y=0.95,label=paste0("AUC:",round(res.roc$auc,3)),size=2,fonteface="blod", family="Times New Roman")
