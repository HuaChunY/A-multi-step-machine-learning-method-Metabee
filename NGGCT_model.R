
NGGCTtrain <- readRDS("NGGCTtrian.rds")

y.train <- substr(rownames(NGGCTtrain),1,3)

NGGCTtest <- readRDS("NGGCTtest.rds")

y.test <-  substr(rownames(NGGCTtest),1,3)

#################

###############linear interpolation

X.aug <- t(NGGCTtrain[1:4,])

label <- substr(colnames(X.aug),1,3) #set the label 
#eset.mat is original data

lab_con <- which(label%in%"Con")
lab_tet <- which(label%in%"tet")

nSample <- length(lab_con) #sum the sample number by groups

nliner <- 10 #linear interpolation number
nfrom <- 0.02# linear interpolation start

augx.M <- c()

for(l in 1:2){
  a <- list(lab_con,lab_tet)
  nSample <- length(a[[l]])
  if(l == 1){
    if(nSample < 2){
      next
    }else{
      index <- 1:(length(lab_con) - 1)}
  }else{
    if(nSample < 2){
      next
    }else{
      index <- (length(lab_con)+1):(length(c(lab_con,lab_tet)) - 1)}
  }
  
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
}


rownames(augx.M)<- rownames(X.aug)


colnames(augx.M) <- paste0("Con",21:(ncol(augx.M)+20))

X.aug <- rbind(t(augx.M),NGGCTtrain)


y.aug <- c(substr(colnames(augx.M),1,3),y.train)

folds <- 10
set.seed(10)

test.fold <- split(sample(1:length(y.aug)), 1:folds) 


#0.7882

gama <- seq(0.1,20,0.05)

C.num <- seq(0.1,20,0.05)

com.auc <- c()
AUC.dt <- c()
All.AUC <- lapply(1:length(C.num),function(G){
  
  all.AUC <- lapply(1:length(gama),function(g){
    set.seed(10)
    all.pred.tables <- lapply(1:folds, function(si) {
      
      test.id <- test.fold[[si]]
      
      X.train <- as.data.frame(X.aug[-test.id, ])
      
      rownames(X.train)<-rownames(X.aug)[-test.id]
      
      colnames(X.train)<-colnames(X.aug)
      
      y.train <- as.factor(y.aug[-test.id])
      
      X.t<-as.data.frame(X.aug[test.id, ])
      
      colnames(X.t)<-colnames(X.aug)
      
      rownames(X.t)<-rownames(X.aug)[test.id]
      
      X1.test <-rbind(NGGCTtest,X.t)
      
      model <- svm(X.train, y.train, 
                   kernel = "radial", 
                   gamma = gama[g], C=C.num[G],cachesize=500,
                   type="C-classification", 
                   prob = TRUE, cross = 10) 
      #cost=10, cachesize=500
      
      predict.test <- predict(model, X1.test, prob = TRUE)
      
      prob.benign <- attr(predict.test, "probabilities")[, 2]
      
      data.frame(y.test =c(as.character(y.test),y.aug[test.id]), y.pred = prob.benign) # returning this
      
    })
    
    full.pred.table <- do.call(rbind, all.pred.tables)
    
    res.roc <- roc(full.pred.table$y.test, 
                   full.pred.table$y.pred)
    
    com.auc.m <- data.frame(gama =gama[g],C=C.num[G], AUC = res.roc$auc) # returning this
    
    rbind(com.auc,as.matrix(com.auc.m))
    
  })
  
  
  AUC.table <- do.call(rbind, all.AUC)
  
  rbind(AUC.dt,as.matrix(AUC.table))
})

All.AUC.table <- do.call(rbind, All.AUC)

################SVM
set.seed(10)
all.pred.tables <- lapply(1:folds, function(si) {
  
  test.id <- test.fold[[si]]
  
  X.train <- as.data.frame(X.aug[-test.id, ])
  
  rownames(X.train)<-rownames(X.aug)[-test.id]
  
  colnames(X.train)<-colnames(X.aug)
  
  y.train <- as.factor(y.aug[-test.id])
  
  X.t<-as.data.frame(X.aug[test.id, ])
  
  colnames(X.t)<-colnames(X.aug)
  
  rownames(X.t)<-rownames(X.aug)[test.id]
  
  X1.test <-rbind(NGGCTtest,X.t)
  
  model <- svm(X.train, y.train, 
               kernel = "radial", 
               gamma = 0.85,cachesize=500,
               type="C-classification", 
               prob = TRUE, cross = 10) 
  #cost=10, cachesize=500
  
  predict.test <- predict(model, X1.test, prob = TRUE)
  
  prob.benign <- attr(predict.test, "probabilities")[, 2]
  
  data.frame(y.test =c(as.character(y.test),y.aug[test.id]), y.pred = prob.benign) # returning this
  
})

full.pred.table <- do.call(rbind, all.pred.tables)

res.roc <- roc(full.pred.table$y.test, 
               full.pred.table$y.pred)



ggroc(res.roc,
      alpha=0.5,colour="#377EB8",
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
  scale_x_continuous(expand = c(0.005,0))+
 annotate("text",x=0.25,y=0.95,label=paste0("AUC:",round(res.roc$auc,3)),size=2,fonteface="blod", family="Times New Roman")

rets <- c("threshold", "specificity", "sensitivity", "accuracy")
pROC::ci.coords(res.roc, x="best", input = "threshold",best.policy="random",ret=rets)
pROC::ci(res.roc)

#################

#################t-SNE

#Perplexity parameter<(ncol(eset.mat) - 1 )/ 3

set.seed(2)
initial_value<-40
theta_value <- 0.01
plex_value <- 10

tsne_out<-Rtsne(X.aug,
                dims = 2,initial_dims = initial_value,
                pca = FALSE,
                perplexity = plex_value,
                theta = theta_value, 
                max_iter = 1000
)

Group <- as.factor(c(paste0("fake_",substr(colnames(augx.M),1,3)),
                     substr(rownames(NGGCTtrain),1,3)))

tsne <- data.frame(tSNE1=tsne_out$Y[,1],
                   tSNE2=tsne_out$Y[,2],
                   Group=Group
)

#View(TSNE.data)
ggplot(tsne,aes(tSNE1,tSNE2))+
  geom_point(aes(color=Group))+
  scale_color_lancet()+
  scale_color_manual(values=c("#99BFAB","#EBDEB3","#00468B","#F89B9B"))+
  theme_bw()+
  xlim(-150,150)+
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


