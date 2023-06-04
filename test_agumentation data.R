setwd("/Users/huachunyin/desktop/Xinqiao/tumor/NO/")

test.mat <- openxlsx::read.xlsx("test.xlsx",rowNames=TRUE)

###############linear interpolation

#eset.mat is original data

label <- substr(colnames(test.mat),1,3) #set the label 

lab_con <- which(label%in%"Con")
lab_tet <- which(label%in%"tet")

nSample <- length(lab_tet) #sum the sample number by groups

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
      
      combn.df2 <- test.mat[,combn.df[t,]]
      
      tmp1 <- combn.df2[,1]
      tmp2 <- combn.df2[,2]
      
      augx <- matrix(tmp1*liner.index[r]+tmp2*liner.index[nliner+1-r], # xtmp1+(1-x)tmp2
             ncol = 1)
      
      colnames(augx) <- paste0(colnames(combn.df2)[1],"+",colnames(combn.df2)[2],"+",r)
      augx
    })
  }
  )
  augx.M <- cbind(augx.M, as.matrix(augx.M.m))
}



rownames(augx.M)<- rownames(test.mat)
#View(augx.M)

#################t-SNE

#Perplexity parameter<(ncol(eset.mat) - 1 )/ 3

set.seed(2)
initial_value<-40
theta_value <- 0.01
plex_value <- 10

X.tsn<- cbind(augx.M,test.mat)

tsne_out<-Rtsne(t(X.tsn),
                dims = 2,initial_dims = initial_value,
                pca = FALSE,
                perplexity = plex_value,
                theta = theta_value, 
                max_iter = 1000
)

Group <- as.factor(c(paste0("fake_",substr(colnames(augx.M),1,3)),
                     substr(colnames(test.mat),1,3)))

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
  xlim(-350,350)+
  ylim(-300,300)+
  theme(plot.margin = unit(rep(1.5,4),"lines"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.9,0.2), # legend???????????????
        legend.background = element_rect(size = 1, colour = "white"))



###############

eset.mat <- openxlsx::read.xlsx("eset.mat.xlsx",rowNames=TRUE)

label <- substr(colnames(eset.mat),1,3) #group 

sam.lab <- as.factor(label)

input <- cbind(sam.lab, as.data.frame(t(eset.mat)))

input <- as.data.frame(input)

#View(input)

##################

#"N-Acetylaspartylglutamate (NAAG)-","Cytidine+","L-Pipecolic acid+",
#"Acamprosate-" "3-Methylhistidine+","S-Methyl-5'-thioadenosine+", "Dodecanoic acid-",
#"Linoleic acid-","L-Glutamate+","S-Methyl-5'-thioadenosine+",
 


deg.cer<-c("N-Acetylaspartylglutamate (NAAG)-","Dodecanoic acid-",

            "β-HCG")

sam.iris <- input[, c("sam.lab", deg.cer)]

head(sam.iris)

sam.iris<-as.data.frame(sam.iris)
#sam.iris2<-sam.iris[,2]

X <- sam.iris[, names(sam.iris) != "sam.lab"]

X<-as.data.frame(sam.iris[,-1])

rownames(X)<-rownames(sam.iris)

colnames(X)<-deg.cer

y <- as.character(sam.iris$sam.lab)

#####test data

X.test.aug<- cbind(augx.M,test.mat$Con)
colnames(X.test.aug)[ncol(X.test.aug)] <- "Con"
y.test1<- substr(colnames(X.test.aug),1,3)
X.test1<- t(X.test.aug[deg.cer,])


folds <- 10
set.seed(10)

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
  
  X1.test <-rbind(X.test1,X.t)
  
  model <- svm(X.train, y.train,  gamma = 0.01,kernel = "radial", prob = TRUE, cross = 10) # some tuning may be needed
  
  predict.test <- predict(model, X1.test, prob = TRUE)
  
  prob.benign <- attr(predict.test, "probabilities")[, 2]
  
  data.frame(y.test =c(as.character(y.test1),y[test.id]), y.pred = prob.benign) # returning this
  
})

full.pred.table <- do.call(rbind, all.pred.tables)


res.roc <- roc(full.pred.table$y.test, 
               full.pred.table$y.pred, 
               plot = TRUE, 
               xlab = "False Positive Percentage", 
               ylab = "TRUE Positive Percentage", 
               legacy.axes = TRUE)
res.roc

P.AUC <-  ggroc(res.roc,alpha=0.5,colour="#377EB8",
                pe=2,size=1,legacy.axes = TRUE)+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               color="grey", linetype="dashed")+
  #ggtitle(deg.cer)+
  theme(text = element_text(family="Times New Roman"),
        plot.title = element_text(family = "serif", #??????
                                  size = 8, #????????????
                                  hjust = 0.5, #?????????????????????
                                  angle = 0, #?????????????????????
        ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=8),
        axis.text = element_text(size=6))+
  xlab("False Positive Percentage")+
  ylab("TRUE Positive Percentage")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0.005,0))+
  annotate("text",x=0.25,y=0.95,label=paste0("AUC:",round(res.roc$auc,3)),
           size=2,fonteface="blod", family="Times New Roman")
P.AUC

pROC::ci(res.roc)

#rets=c("threshold", "specificity", "sensitivity","accuracy")
pROC::ci.coords(res.roc, x="best", input = "threshold",best.policy="random",ret=rets)

##############################

dataAUC<-openxlsx::read.xlsx("agumentation_marker.xlsx")

#######TOPSIS全称Technique for Order Preference by Similarity to an Ideal Solution
library(tibble)
library(dplyr)
library(readr)

# 标准化变量值
z_value <- function(x){
  x / sqrt(sum(x^2))
}

# 计算最优距离
dist <-function(x, std){
  res <- c()
  for ( i in 1 : nrow(x)) {
    res[i] = sqrt(sum((unlist(x[i,-1])-std)^2))
  }
  
  return(res)
}

# load sample data

# 按列对数据进行标准化

dat_z <- dataAUC %>% dplyr::mutate(across(c(3:7), z_value))
dat_z <- dat_z[,-1]
## unlist 转换tibble为vector
z_max <- dat_z %>% summarise(across(c(2:6), max)) %>% unlist
z_min <- dat_z %>% summarise(across(c(2:6), min)) %>% unlist

# dat_z %>% select(2:4) %>% rowwise() %>% mutate(du = dist(., z_max), dn= dist(., z_min)) 
du <- dist(dat_z, z_max)
dn <- dist(dat_z, z_min)

# 计算CI并按照降序排序
dat_z <- dat_z %>% add_column(du = du, dn = dn) %>% 
  mutate(ci= dn/(du+dn)) %>%
  arrange(-ci)



