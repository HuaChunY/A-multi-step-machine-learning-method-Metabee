

source("Chunk01-R script for searching feature.R")

k1 <- 5 #folds for outer loop
k2 <- 5 # folds for inner loop

#######################################################################

#01 read the data

#######################################################################


NGGCTtrain <- readRDS("NGGCTtraining.rds")
NGGCTtest <- readRDS("NGGCTtest.rds")


######################################################################

#02 linear interpolation

######################################################################

X.aug <- NGGCTtrain[,1:4]

label <- substr(colnames(X.aug),1,3) #set the label 

lab_con <- which(label%in%"Con")

nSample <- length(lab_con) #sum the sample number by groups

nliner <- 10 #linear interpolation number
nfrom <- 0.02# linear interpolation start

augx.M <- c()

index <- lab_con

combn.df <- t(combn(index,2))
  
  liner.index <- seq(nfrom,(1-nfrom),length.out=nliner)
  
  augx.M.m <- purrr::map_dfc(1:nrow(combn.df), .f = function(t){
    
    purrr::map_dfc(1:nliner, .f = function(r){
      
      combn.df2 <- X.aug[,combn.df[t,]]
      
      tmp1 <- combn.df2[,1]
      tmp2 <- combn.df2[,2]
      
      augx <- matrix(tmp1*liner.index[r]+tmp2*liner.index[nliner+1-r], # xtmp1+(1-x)tmp2
                     ncol = 1)
      
      colnames(augx) <- paste0(colnames(combn.df2)[1],"+",colnames(combn.df2)[2],"+",r)
      augx
    })
  })
  augx.M <- cbind(augx.M, as.matrix(augx.M.m))



rownames(augx.M)<- rownames(X.aug)


colnames(augx.M) <- paste0("Con",21:(ncol(augx.M)+20))

X.aug <- cbind(augx.M,NGGCTtrain)

y.aug <- substr(colnames(X.aug),1,3)


#################t-SNE

#Perplexity parameter<(ncol(eset.mat) - 1 )/ 3

set.seed(2)
initial_value<-50
theta_value <- 0.01
plex_value <- 10

tsne_out<-Rtsne::Rtsne(t(X.aug),
                       dims = 2,initial_dims = initial_value,
                       pca = FALSE,
                       perplexity = plex_value,
                       theta = theta_value, 
                       max_iter = 1000
)

Group <- as.factor(c(paste0("fake_",substr(colnames(augx.M),1,3)),
                     substr(colnames(NGGCTtrain),1,3)))

tsne <- data.frame(tSNE1=tsne_out$Y[,1],
                   tSNE2=tsne_out$Y[,2],
                   Group=Group)

#View(TSNE.data)
ggplot(tsne,aes(tSNE1,tSNE2))+
  geom_point(aes(color=Group))+
  scale_color_lancet()+
  scale_color_manual(values=c("#99BFAB","#EBDEB3","#00468B","#F89B9B"))+
  theme_bw()+
  xlim(-50,50)+
  ylim(-100,100)+
  theme(plot.margin = unit(rep(1.5,4),"lines"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.9,0.2), # legend???????????????
        legend.background = element_rect(size = 1, colour = "white"))

############################
###########################
X.aug1 <- X.aug[c("Inosine-","N-Acetylaspartylglutamate (NAAG)+","β-HCG"),]

View(X.aug1)
X.aug1 <- t(X.aug1)

NGGCTtest1 <- NGGCTtest[c("Inosine-","N-Acetylaspartylglutamate (NAAG)+","β-HCG"),]

NGGCTtest1 <- t(NGGCTtest1)

y.test <- str_replace_all(substr(rownames(NGGCTtest1),1,3),
                          "Ger","tet")


###########################################################################

#03 SVM

###########################################################################

folds <- k1

rets <- c("threshold", "specificity", "sensitivity", "accuracy")

all.AUC.tables <- lapply(1:100, function(t) {
  
  set.seed(t)
  ############
  
  test.fold <- split(sample(1:length(y.aug)), 1:folds) 
  
  set.seed(10)
  
  all.pred.tables <- lapply(1:folds, function(si) {
    
    test.id <- test.fold[[si]]
    
    X.train <- as.data.frame(X.aug1[-test.id, ])
    
    rownames(X.train)<-rownames(X.aug1)[-test.id]
    
    colnames(X.train)<-colnames(X.aug1)
    
    y.train <- as.factor(y.aug[-test.id])
    
    X.t<-as.data.frame(X.aug1[test.id, ])
    
    colnames(X.t)<-colnames(X.aug1)
    
    rownames(X.t)<-rownames(X.aug1)[test.id]
    
    X1.test <-rbind(NGGCTtest1,X.t)
    set.seed(10)
    model <- svm(X.train, y.train, 
                 kernel = "radial", 
                 gamma = 0.85,cachesize=500,
                 type="C-classification", 
                 prob = TRUE, cross = k2) 
    #cost=10, cachesize=500
    
    predict.test <- predict(model, X1.test, prob = TRUE)
    
    prob.benign <- attr(predict.test, "probabilities")[, 2]
    
    data.frame(y.test =c(as.character(y.test),as.character(y.aug[test.id])), y.pred = prob.benign) # returning this
    
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

NGGCT.AUC <- do.call(rbind, all.AUC.tables)

colnames(NGGCT.AUC) <- paste0(c(rep("AUC",3),
                                rep("threshold",3),
                                rep("specificity",3),
                                rep("sensitivity",3),
                                rep("accuracy",3)),
                              rep(c(".low",".median",".high"),5))


NGGCTAUC <-  rbind(NGGCT.AUC,apply(NGGCT.AUC,2,mean))

rownames(NGGCTAUC) <- c(paste0(rep("seed",100),1:100),"Mean")

#########################

library(ggplot2)
data.merge <- NGGCTAUC[-101,c("AUC.median","specificity.median", "sensitivity.median","accuracy.median" )]

data.merge2<-reshape2::melt(data.merge,
                            variable.name = "name", 
                            value.name = "value")
ggplot(data=data.merge2, aes(x= Var2,y= value))+
  geom_boxplot(aes(color=Var2),size=0.5, outlier.size = 0.1)+
  scale_color_manual(values = c("#336699","#CBBB7F","#D15C62","#669999"))+
  #facet_wrap(~ Name, scales="free",nrow =5)+
  theme_bw()+
  #ylab("AUC")+
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

#################

