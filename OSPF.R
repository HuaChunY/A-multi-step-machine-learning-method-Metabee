################################################################################
#    &&&....&&&    % Project: a multi-step machine learning method-OSPF        #
#  &&&&&&..&&&&&&  % Author: Bo Li, Huachun Yin                                #
#  &&&&&&&&&&&&&&  % Date: Mar. 1st, 2022                                      #
#   &&&&&&&&&&&&   %                                                           #
#     &&&&&&&&     % Environment: R version 3.5.3;                             #
#       &&&&       % Platform: x86_64-pc-linux-gnu (64-bit)                    #
#        &         %                                                           #
################################################################################
 
source("Chunk01-R script for RNA-seq data simulation based on improved compcodeR.R")
# the source download from https://github.com/libcell/MSPJ
library(preprocessCore)
library(ggplot2)
library(pROC)
library(FSinR)
library(FSA)
library(caret)
library(purrr)
library(ropls)
library(dplyr)
library(ggdendro)
library(ggpubr)
library(Rtsne)
library(ggsci)
library("FactoMineR")
library("factoextra")
library(RankAggreg)
library(ggsignif)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(pheatmap)
library(psych)
library(pheatmap)



# save pdf
save_pdf <- function(x, filename, width=6, height=4) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# Calculating with many distance metrics.
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


#################t-SNE

iris_unique1<-openxlsx::read.xlsx("POS1.xlsx",rowNames=TRUE)# the positive ion modes
iris_unique2<-openxlsx::read.xlsx("NEG1.xlsx",rowNames=TRUE)# the negative ion modes


m<-t(iris_unique2)

sam.lab<-substr(rownames(m),1,3)

sam.lab1<-c(rep("Con",22),
            str_replace_all(substr(rownames(m)[23:46],1,3),
                            "Con","NGGCT"))

c.p <- which(sam.lab1%in%"Con")
e.p <- which(sam.lab1%in%"tet")

m<-m[c(e.p,c.p),]

sam.lab<-substr(rownames(m),1,3)


set.seed(2)
initial_value<-40
theta_value <- 0.01
plex_value <- 12
tsne_out<-Rtsne(log2(m+1),
                dims = 2,initial_dims = initial_value,
                pca = FALSE,
                perplexity = plex_value,
                theta = theta_value, 
                max_iter = 1000
                )

tsne <- data.frame(tSNE1=tsne_out$Y[,1],
                   tSNE2=tsne_out$Y[,2],
                   Group=sam.lab)

  ggplot(tsne,aes(tSNE1,tSNE2))+
  geom_point(aes(color=Group))+
  scale_color_lancet()+
  scale_color_manual(values=c("#00468B","#F89B9B"))+
  theme_bw()+
  xlim(-150,150)+
  ylim(-50,50)+
  theme(plot.margin = unit(rep(1.5,4),"lines"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.9,0.2), # legend???????????????
        legend.background = element_rect(size = 1, colour = "white"))



tsne_out$Y[1,]
#[1]  34.21848 -39.11180


iris.pca <- PCA(log2(m+1), graph = FALSE)

ind.p <-fviz_pca_ind(iris.pca,
                     geom.ind = "point",
                     col.ind = sam.lab,
                     palette = "simpsons",
                     addEllipses = TRUE,
                     ellipse.level = 0.80,
                     legend.title = "Groups"
                     
)
 

ggpubr::ggpar(ind.p,
              legend.title = "Group",palette = c("#00468B","#F89B9B"))+
              theme(plot.margin = unit(rep(1.5,4),"lines"),
                    panel.grid = element_blank(),
                    axis.title = element_blank(),
                    axis.ticks = element_blank(),
                    axis.ticks.x= element_blank(),
                    axis.text = element_blank(),
                    panel.border = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    legend.position = c(0.9,0.2), # legend???????????????
                    legend.background = element_rect(size = 1, colour = "white"))
              

#####################

rawPOS <-openxlsx::read.xlsx("POS.xlsx") #POS value of train
rawNEG <-openxlsx::read.xlsx("NEG.xlsx")#NEG value of train 

POStest<-openxlsx::read.xlsx("POStest.xlsx")#POS value of test
NEGtest<-openxlsx::read.xlsx("NEGtest.xlsx")#POS value of test

NEGananotaion<-openxlsx::read.xlsx("NEGananotaion.xlsx") #ananotion information
POSannotation<-openxlsx::read.xlsx("POSannotation.xlsx")

Manao<-rbind(POSannotation,NEGananotaion)

Man.pos<-which(duplicated(Manao$description)==FALSE) 

Manao<-Manao[Man.pos,]

Clinicdt <-openxlsx::read.xlsx("Clinic.xlsx")

colnames(rawPOS)

sam.lab <- substr(colnames(rawNEG),1,3)

sam.lab <- sam.lab[-(1:2)]

M<-c("rawPOS","rawNEG")

for (G in 1:2) {
  
  rawmatrix<-get(M[G])
  
  c.d <- which(sam.lab == "Con")
  e.d <- which(sam.lab == "tet")
  
  M_matrix<-rawmatrix[,2:ncol(rawmatrix)]

  POS.NA<-which(is.na(M_matrix$description)==FALSE)
    
  M_matrix<-M_matrix[POS.NA,]
    
  rownames(M_matrix)<-M_matrix$description
    
  M_matrix<-M_matrix[,-1]
  
  M_matrix<-M_matrix[,c(c.d,e.d)] 
  
  sam.lab <- substr(colnames(M_matrix),1,3) 

  sam.lab<-as.factor(sam.lab)
    
   #########FC
   
   c.p <- which(sam.lab == "Con")
   e.p <- which(sam.lab == "tet")
   
   FC <-apply(M_matrix, 1, function(x) {mean(x[e.p])/mean(x[c.p])})
   
   FC <-as.data.frame(FC)
   
   FC$"log2(FC)"<-log2(FC$FC)
   
   FC$FeatureName<-rownames(FC)
   
   ##########

   M_matrix<-log2(M_matrix+0.1)
   
   M_matrix <-as.data.frame(t(M_matrix))
   
   input<- cbind(sam.lab, M_matrix)
   
   input <- as.data.frame(input)

   
   ###############OPLSDA
   
   set.seed(1)
   
   oplsda = opls(input[,-1], sam.lab, predI = 1, orthoI = NA, crossvalI=10)
   
   sample.score <- oplsda@scoreMN %>%as.data.frame() %>% mutate(label = sam.lab,
                                                                o1 = oplsda@orthoScoreMN[,1])
   
   vip <- getVipVn(oplsda)
   
   vip <- as.data.frame(vip)
   
   vip$FeatureName<-rownames(vip)
   
   ##########SVM-RFE#####
   
   nfold <- 10
   
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
   
   assign(paste0("deg.svm.",M[G]),deg.svm)
   
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
     
     top.features1<-top.features[1:100,]
     
     top.features1<-merge(FC,top.features1)
     
     top.features1<-merge(top.features1,deg.count1)
     
     sumdata<-merge(top.features1,vip)
     
     sumdata<-sumdata[sumdata$vip>1 & sumdata$Pvalue < 0.05,]
     #
     sumdata<-sumdata[sumdata$FC>1.3 | sumdata$FC<0.6,]
     
     assign(paste0("DEM.",M[G]),sumdata)
     
}

DEM1<-rbind(DEM.rawNEG,DEM.rawPOS)

DEM1 <- DEM1 %>% 
  mutate(
    Ion=c(rep("-",nrow(DEM.rawNEG)),rep("+",nrow(DEM.rawPOS)))
    )


hDEM<-merge(DEM1,Manao,by.x = "FeatureName",by.y = "description")


openxlsx::write.xlsx(hDEM,"DEMannotaion.xlsx")

########

hNEG<-hDEM[which(hDEM$Ion=="-"),]
hPOS<-hDEM[which(hDEM$Ion=="+"),]
hNEG<-merge(hNEG,rawNEG,by.x = "FeatureName",by.y = "description")
hPOS<-merge(hPOS,rawPOS,by.x = "FeatureName",by.y = "description")

hDEG<-rbind(hNEG,hPOS)


hDEG <- hDEG %>% 
  mutate(
    Comp=case_when(
      Ion=="+" ~ paste0(hDEG$FeatureName,"+"),
      
      Ion=="-" ~ paste0(hDEG$FeatureName,"-")
    )
  )

hDEG1<-hDEG[,17:62]

hDEG1<-as.data.frame(lapply(hDEG1, as.numeric))

hDEG1<-log2(hDEG1)

rownames(hDEG1)<-hDEG$Comp

hDEG1<-hDEG1[,c(c.d,e.d)]#####delet the NGGCT

annorow<-hDEG[,c(12:13)]

rownames(annorow)<-rownames(hDEG1)

annocol<-as.data.frame(sam.lab)###sam.lab1

rownames(annocol)<-colnames(hDEG1)

colnames(annocol)<-"Group"
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

col_colors<-get_anno_for_heatmap2(annorow,annocol)


library(d3heatmap)
library(philentropy)

meth<-getDistMethods()


n<-9
m<-meth[n]

a<-hclust(dist_n(t(hDEG1),p=1,mtd=m))

df2<-dendro_data(a,type="rectangle")

Group<- substr(df2$labels$label,1,3)
ggplot(segment(df2))+
  geom_text(data=df2$labels,aes(x=x,y=y,label=label,color=Group),
            angle=90,hjust=1,vjust=0.3,size=3)+
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend))+
 theme(panel.background = element_rect(fill = "white"))


############heatmap

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
  
  
  
# if(only.color){
#   anno_res<-annocolor
# }else{anno_res<-list(annocol=annocol,
#                      annorow=annorow,
#                      annocolor=annocolor)}
# return(anno_res)
  
  
  
}

col_colors<-get_anno_for_heatmap2(annorow,annocol)

################
tree_name <- df2$labels$label


pheatmap::pheatmap(hDEG1[,as.character(tree_name)],
            show_rownames = T,
            show_colnames = T,
            cluster_cols = F,
            cluster_rows= T,
            fontsize_row=6, #??????????????????
            border_color = "NA",
            scale = "row",
            angle_col=45, #??????????????????
            annotation_row = annorow,
            annotation_col = annocol,
            annotation_colors = col_colors,
            color =colorRampPalette(c("#0225A2", "white","#FF2400"))(100), 
            clustering_distance_rows = 'maximum', 
            clustering_distance_cols='manhattan', 
            clustering_method = 'ward.D2'
)

#"ward", "single", "complete", "average", "mcquitty", "median" or "centroid".
#########################pie

annorow <- na.omit(annorow)

df <- group_by(annorow, Class)%>%
  dplyr::summarize(percent = n() / nrow(annorow)) %>%
  arrange(desc(percent))%>%mutate(col=col_colors$Class[Class])


pie(df$percent,labels = df$Class)

pie(df$percent, col = df$col,
    labels = with(df, paste0(round(percent, 2) * 100, "%")))


df1 <- group_by(annorow,SuperClass )%>%
  dplyr::summarize(percent = n() / nrow(annorow)) %>%
  arrange(desc(percent))%>%mutate(col=col_colors$SuperClass[SuperClass])

pie(df1$percent,labels = df1$SuperClass)

pie(df1$percent, col = df1$col,
    labels = with(df1, paste0(round(percent, 2) * 100, "%")))

######################


#############normalizated the metabolites and elisa resultes

elsia1<-openxlsx::read.xlsx("threeelsia.xlsx",rowNames=TRUE)

elsia<-elsia1[,c(c.d,e.d)]###delet the NGGCT

colnames(elsia)<-colnames(hDEG1)

nor.eset<-rbind(hDEG1,elsia)

lable<-colnames(nor.eset)

lable<-substr(lable,1,3)

nor.eset1 <-log(nor.eset+1)#[,-(39)]

tmp1 <- normalize.quantiles(as.matrix(nor.eset1))

colnames(tmp1)<-colnames(nor.eset1)
rownames(tmp1)<-rownames(nor.eset1)

mypal1 <- terrain.colors(ncol(nor.eset))

boxplot(tmp1, col = mypal1,outline=F)


################SVM########


##############preparing training data
View(hDEG)
hDEG_order <- hDEG %>% 
  group_by(Ion) %>% 
  arrange(AvgRank) %>% # or  desc(AvgRank)
  slice(1:10)# gain the value by index 

dem<-c(hDEG_order$Comp,rownames(elsia))

eset.mat <- tmp1[dem,]

eset.mat <- as.data.frame(t(eset.mat))

input <- cbind(sam.lab, eset.mat)

input <- as.data.frame(input)

##############################

##############preparing testing data

elsiatest<-openxlsx::read.xlsx("elisatest.xlsx",rowNames=TRUE)

rownames(POStest) <- paste0(POStest$Name,"+")
POStest <- POStest[,-1]

rownames(NEGtest) <- paste0(NEGtest$Name,"-")

NEGtest <- NEGtest[,-1]

test_dt<-rbind(POStest,NEGtest,elsiatest)

test_dt<-test_dt[dem,]

test_dt <- na.omit(test_dt)


###########################

##########################

X.test1<- test_dt

X.test2 <-log(X.test1+1)

tmp2 <- normalize.quantiles(as.matrix(X.test2))

colnames(tmp2)<-colnames(X.test1)
rownames(tmp2)<-rownames(X.test1)

#boxplot(tmp2, col = mypal1,outline=F)


for (GM in 1:length(dem)) {
 
  deg.cer<-c("Dodecanoic acid-",
             "N-Acetylaspartylglutamate (NAAG)-",
             "β-HCG" ) 
  
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
  X.test<-tmp2[,c(9:13)]#test group: 4 germinoma
  
  # colnames(X.test2)
  
  X.test<-t(X.test)
  
  y.test1<-rownames(X.test)
  y.test1<-substr(y.test1,1,3)
  
  y.test1<-str_replace_all(y.test1,c("con"="Con","ger"="tet"))
  
  X.test<-as.data.frame(X.test[,deg.cer])
 
  colnames(X.test)<-deg.cer
 
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
    
    X1.test <-rbind(X.test,X.t)
    
    model <- svm(X.train, y.train, kernel = "radial", prob = TRUE, cross = 10)  
    
    predict.test <- predict(model, X1.test, prob = TRUE)
    
    prob.benign <- attr(predict.test, "probabilities")[, 2]
    
    data.frame(y.test =c(as.character(y.test1),y[test.id]), y.pred = prob.benign)  
    
  })
  
  full.pred.table <- do.call(rbind, all.pred.tables)
  
  
  
  res.roc <- roc(full.pred.table$y.test, 
                 full.pred.table$y.pred, 
                 plot = TRUE, 
                 xlab = "False Positive Percentage", 
                 ylab = "TRUE Positive Percentage", 
                 legacy.axes = TRUE)
  
  P.AUC <-  ggroc(res.roc,alpha=0.5,colour="#377EB8",
         pe=2,size=1,legacy.axes = TRUE)+
         geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                  color="grey", linetype="dashed")+
        ggtitle(deg.cer)+
        theme(text = element_text(family="serif"),
              plot.title = element_text(family = "serif", 
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
          annotate("text",x=0.25,y=0.95,label=paste0("AUC:",round(res.roc$auc,3)),
                   size=2,fonteface="blod", family="serif")
   P.AUC
  
   assign(paste0("AUC",deg.cer),P.AUC)
 
   feedback<-try(pROC::ci.coords(res.roc, x="best", input = "threshold",best.policy="random",ret=rets),silent=TRUE)
  
   error_or<-"Error in enforce.best.policy(res, best.policy) : \n  More than one \"best\" threshold was found, aborting. Change 'best.policy' to alter this behavior.\n"
  
  if(error_or%in%feedback==TRUE){
    
    feedback1<-try(pROC::ci.coords(res.roc, x="best", input = "threshold", best.policy="random",ret=rets),silent=TRUE)
    
     if(error_or%in%feedback1==TRUE){
      
      thresholdx<-0.9
      ciorder<-pROC::ci.coords(res.roc, x=thresholdx, input = "threshold", best.policy="random",ret=rets)
      
    
      }else{
       
       ciorder<-feedback1
     }
     
  
    }else{
    
    ciorder<-feedback
    
  }
  
  order_ci<-c(as.character(deg.cer),pROC::ci(res.roc),ciorder)
  
  order_ci<-t(as.matrix(unlist(order_ci)))
  
  write.table(order_ci,"ciorderbset.txt",row.names = FALSE, col.names = FALSE, append = TRUE)
  
  
}



  ggarrange(get(paste0("AUC",dem[1])),get(paste0("AUC",dem[2])),
            get(paste0("AUC",dem[3])),get(paste0("AUC",dem[4])),
            get(paste0("AUC",dem[5])),get(paste0("AUC",dem[6])),
            get(paste0("AUC",dem[7])),get(paste0("AUC",dem[8])),
            get(paste0("AUC",dem[9])),get(paste0("AUC",dem[10])),
            get(paste0("AUC",dem[11])),get(paste0("AUC",dem[12])),
            get(paste0("AUC",dem[13])),get(paste0("AUC",dem[14])),
            get(paste0("AUC",dem[15])),get(paste0("AUC",dem[16])),
            get(paste0("AUC",dem[17])),get(paste0("AUC",dem[18])),
            get(paste0("AUC",dem[19])),get(paste0("AUC",dem[20])),
            get(paste0("AUC",dem[21])),get(paste0("AUC",dem[22])),
            get(paste0("AUC",dem[23])),P.AUC,ncol=4,nrow=6)


########boxplot
elsia<-openxlsx::read.xlsx("threeelsia.xlsx",rowNames=TRUE)
 
colnames(elsia)<-colnames(hDEG1)

dim(elsia)

samplot<-hDEG_order[,17:(ncol(hDEG_order)-1)]

samplot$tet24<-as.numeric(samplot$tet24)

samplot <- samplot#[,c(c.d,e.d)]

rownames(samplot)<-hDEG_order$Comp

dim(samplot)

row.pos <- which(rownames(samplot)%in%c("N1-Methyl-2-pyridone-5-carboxamide-",
                                        "N-Acetylaspartylglutamate (NAAG)-",
                                        "N1-Methyl-2-pyridone-5-carboxamide+",
                                        "N-Acetylaspartylglutamate (NAAG)+",
                                        "Dodecanoic acid-"))

rownames(samplot)[row.pos] <- c("PY2-","NAAG-","DAA-","PY2+","NAAG+")


samplot<-as.data.frame(t(samplot))

samplot$Group<-sam.lab1#sam.lab

samplot<-reshape2::melt(samplot,variable.name = "Name")

colnames(samplot)<-c("Group","Name","Expression")

samplot$Expression<-log2(samplot$Expression)

elsia1<-as.data.frame(t(elsia))

elsia1$Group<-sam.lab1#sam.lab

elsia1<-reshape2::melt(elsia1,variable.name = "Name")

colnames(elsia1)<-c("Group","Name","Expression")

elsia1$Expression<-log2(elsia1$Expression+1)

samplot1<-rbind(samplot,elsia1)

sig<-compare_means(Expression~Group, data=samplot1,group.by = "Name",method = "t.test")
View(sig)
samplot1 <- samplot1 %>% 
  dplyr::mutate(
    namelength=str_length(as.character(samplot1$Name))
    )
  
samplot1 <- samplot1 [order(samplot1$namelength),]

p<-ggplot(data=samplot1, aes(x=reorder(Name, namelength),y=Expression))+
  geom_boxplot(aes(color=Group),size=0.8, outlier.size = 0.1)+
  scale_color_manual(values = c("#336699","#D15C62","#CBBB7F"))+
  stat_compare_means(aes(group = Group),
                     method = "wilcox.test",
                     label = "p.signif",
                     label.x.npc = "right")+
 # facet_wrap(~ Name, scales="free")+
  theme_bw()+
  scale_fill_discrete(labels=c("Control", "NGGCT","Germinoma"))+
  theme(#text = element_text(family="serif",size = 8),
        panel.grid = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=9,angle = 45),
        axis.text.y = element_text(size=6),
        legend.position ="top",
        legend.text=element_text(size=11),
        legend.title = element_text(size=11)
  )

p

sam_loc<-hDEG_order[,17:(ncol(hDEG_order)-1)]

sam_loc <- sam_loc[,c(c.d,e.d)]

rownames(sam_loc)<-hDEG_order$Comp

sam_loc<-as.data.frame(t(sam_loc))

row.pos <- which(colnames(sam_loc)%in%c("N1-Methyl-2-pyridone-5-carboxamide-",
                                        "N-Acetylaspartylglutamate (NAAG)-",
                                        "N1-Methyl-2-pyridone-5-carboxamide+",
                                        "N-Acetylaspartylglutamate (NAAG)+",
                                        "Dodecanoic acid-"))

colnames(sam_loc)[row.pos] <- c("PY2-","NAAG-","DAA-","PY2+","NAAG+")


sam_loc$Group<-sam.lab

sam_loc<-sam_loc[which(sam_loc$Group%in%"tet"),]

loc.label<-c("Suprasellar",
             "Multiple regions",
             "Multiple regions",
             "Multiple regions",
             "Multiple regions",
             "Suprasellar",
             "Suprasellar",
             "Pineal",
             "Multiple regions",
             "Multiple regions",
             "Multiple regions",
             "Multiple regions",
             "Pineal",
             "Multiple regions",
             "Suprasellar",
             "T",
             "Pineal",
             "Suprasellar",
             "Multiple regions",
             "Pineal"
)


sam_loc$Group<-loc.label

sam_loc<-sam_loc[-which(sam_loc$Group%in%"T"),]

sam_loc <- tidyr::pivot_longer(sam_loc, cols = 1:(ncol(sam_loc)-1), names_to = "Name", values_to = 'Expression')

sam_loc$Expression<-as.numeric(sam_loc$Expression)

sam_loc$Expression<-log2(sam_loc$Expression)


sig<-compare_means(Expression~Group, data=sam_loc,group.by = "Name",method = "t.test")

sam_loc <- sam_loc %>% 
  dplyr::mutate(
    namelength=str_length(as.character(sam_loc$Name))
  )

sam_loc <- sam_loc [order(sam_loc$namelength),]

p1<-ggplot(data=sam_loc, aes(x=reorder(Name, namelength),y=Expression))+
  geom_boxplot(aes(color=Group),size=0.8, outlier.size = 0.1)+
  scale_color_manual(values = c("#336699","#CBBB7F","#D15C62"))+
  stat_compare_means(aes(group = Group),
                     method = "wilcox.test",
                     label = "p.signif",
                     label.x.npc = "right")+
  #facet_wrap(~ Name, scales="free")+
  theme_bw()+
  theme(text = element_text(family="serif",size = 8),
        panel.grid = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=9, angle = 45),
        axis.text.y = element_text(size=6),
        legend.position ="top",
        legend.text=element_text(size=11),
        legend.title = element_text(size=11)
  )

p1


ggarrange(p,p1,nrow = 2,ncol = 1)

###############

dataAUC<-openxlsx::read.xlsx("marker.xlsx")# load sample data

#######TOPSIS:Technique for Order Preference by Similarity to an Ideal Solution
library(tibble)
library(dplyr)
library(readr)

 
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

###############sample bias#########

sigMarker<-dataAUC$Signature[1]

sigMarker <- hDEG$Comp
sig_df<-hDEG[which(hDEG$Comp%in%sigMarker),17:62]

sig_df<-sig_df[,-which(sam.lab1%in%"NGGCT")]

rownames(sig_df) <- sigMarker


#sig_df<-rbind(sig_df,elsia[2,])

sig_df<-as.data.frame(t(sig_df))

sig_df$Sample<-rownames(sig_df)

rownames(Clinicdt)<-Clinicdt$Sample

sig_df <- merge(sig_df,Clinicdt,by="Sample")

sig_df<-as.data.frame(sig_df)

rownames(sig_df)<-sig_df[,1]
tumorsize <- sig_df
sig_df <- sig_df[,1:11]


sig_df1<-reshape2::melt(sig_df,
                        id = c("Sample","Gender","Year","BMI","Hydrocephalus","Group1","Group2"))

colnames(sig_df1)[9:10]<-c("Name","Expression")

sig_df2 <- sig_df1[which(sig_df1$Group2%in%c("Con","tet")),]

sig_df2$Expression<-log2(as.numeric(sig_df2$Expression)+1)

sig_df2 <- na.omit(sig_df2)


library(ggplot2)

library(RColorBrewer)

rhg_cols1 <- c("#CD2626","#53868B","#E9967A",
               "#8B864E","#FFE7BA","#9370DB","#8B636C","#458B74")


Clinicdt <- openxlsx::read.xlsx("Clinic2.xlsx")

sig_df2<-Clinicdt[,c("Gender","Group1")]

sig_df2<-table(sig_df2)

sig_df2<-reshape2::melt(sig_df2,value.name = "value")


sig_df2 <- sig_df2 %>% 
  mutate(Proportion=case_when(Group1== "con" ~ value/22,
                              Group1== "tet" ~ value/20)
  )



pgender<-ggplot(sig_df2,aes(x=Group1,y=Proportion,fill=Gender))+
         geom_bar(stat="identity",position="stack", color="white", width=0.4,size=0.25)+
         coord_flip()+
         scale_fill_manual(values = c("#E58686","#7795B5"))+
         scale_x_discrete(breaks=c("con", "tet"),
                         labels=c("Control","Germinoma"))+
         labs(x=NULL)+
         geom_text(aes(label=value),position = position_stack(vjust = .5),size=3,family="serif")+
         theme_classic()+
         theme(text = element_text(family="serif",size = 10))
         

#######

sig_df3<-Clinicdt[,c("Hydrocephalus","Group1")]

sig_df3<-table(sig_df3)

sig_df3<-reshape2::melt(sig_df3,value.name = "value")

sig_df3 <- sig_df3 %>% 
  mutate(Proportion=case_when(Group1== "con" ~ value/22,
                              Group1== "tet" ~ value/19)
  )

pHydro <-ggplot(sig_df3,aes(x=Group1,y=Proportion,fill=Hydrocephalus))+
         geom_bar(stat="identity",position="stack", color="white", width=0.4,size=0.25)+
         coord_flip()+
         scale_fill_manual(values = c("#E58686","#7795B5"))+
         scale_x_discrete(breaks=c("con", "tet"),
                          labels=c("Control","Germinoma"))+
         labs(x=NULL)+
         geom_text(aes(label=value),position = position_stack(vjust = .5),size=3,family="serif")+
         theme_classic()+
         theme(text = element_text(family="serif",size = 10))

ggarrange(page,pgender,pHydro,nrow = 1,ncol = 3)


#########
  sig_df <- sig_df[,c("Sample","tumor.size")]
  sig_df <- sig_df %>% 
    mutate(tumortype=case_when(as.numeric(tumor.size) < 4 ~ "<4",
                               as.numeric(tumor.size) > 4  ~ ">4",
                               TRUE ~ 'NA')
    )
  
  
  sig_cor1<-sig_df[,c("Sample","tumor.size","tumortype")]
  
  cor_1<-hDEG[,39:63]

  cor_1<-cor_1[, -which(sam.lab1[23:46]%in%"NGGCT")]

  rownames(cor_1)<-cor_1[,"Comp"]
  
  cor_1 <- cor_1[,-21]
  
  cor_1<-as.data.frame(lapply(cor_1, as.numeric))
  
  cor_1<-t(cor_1)
  
  colnames(cor_1) <- hDEG$Comp
  
  sig_cor1 <- sig_df
  sig_cor1 <- sig_cor1 %>% 
    transmute(
      "<4" = case_when(
        tumortype == ">4" ~ 0,
        tumortype == "<4" ~ 1
      ),
      ">4" = case_when(
        tumortype == ">4" ~ 1,
        tumortype == "<4" ~ 0
      )
    )
  
    
  cor_matrix <- corr.test(cor_1, sig_cor1, use = "pairwise", method = "spearman",adjust="holm")
  cor_p<-cor_matrix$p.adj
  
  cor_p<-round(cor_p,3)
  
  cor_p<-as.data.frame(cor_p)
  
  cor_p<-cor_p %>% 
    mutate(
      ">4"=paste0("p.adj: ",cor_p[,1]),
      "<4"=paste0("p.adj: ",cor_p[,2])
    )
  
  b=cor_matrix$r

  p<-pheatmap(b,
              show_rownames = T,
              show_colnames = T,
              cluster_cols = F,
              cluster_rows=T,
              fontsize_row=5, #??????????????????
              border_color = "#FFFFFF",
              scale = "none",
              #annotation_row = annorow,
              angle_col=0, #??????????????????
              color =colorRampPalette(c("#0225A2", "white","#FF2400"))(20),
              clustering_distance_rows = 'euclidean', 
              clustering_method = 'single',
              display_numbers = cor_p,
              fontsize_number = 3.5,
              number_format = "%.2f", number_color = "black"
            
  )  

#########

#############survival 
  
  sur_data<-sam.iris[sam.iris$sam.lab=="tet",]
  
  sur_time <- openxlsx::read.xlsx("survival.xlsx")
  
  View(sur_data)
  
  #summary(sur_data)
  
  
  sur_data1<-sur_data %>% 
    mutate(NAAG=case_when(`N-Acetylaspartylglutamate (NAAG)-` < median(sur_data$`N-Acetylaspartylglutamate (NAAG)-`) ~"L",
                 TRUE ~"H"),
       
        DDA=case_when(
        `Dodecanoic acid-` < median(sur_data$`Dodecanoic acid-`) ~"L",TRUE ~"H"),
        cy=case_when(
        `Cytidine+` < median(sur_data$`Cytidine+`) ~"L",TRUE ~"H"),
        HCGG=case_when(
        `β-HCG` < median(sur_data$`β-HCG`) ~"L",TRUE ~"H"),
      LA=case_when(
        `Linoleic acid-` < median(sur_data$`Linoleic acid-`) ~"L",TRUE ~"H")
      
    )

  sur_data2 <-cbind(sur_data1,sur_time[which(substr(sur_time$sample,1,3)%in%"tet"),])
  
  colnames(sur_data2)  
  
library(survival)
library(survminer) 
  
colnames(sur_data2)

for (sur1 in 1:length(deg.cer)) {
  

  fit <- survfit(Surv(time, status) ~ NAAG, data = sur_data2)
  
  surv_diff <- survdiff(Surv(time, status) ~ NAAG, data = sur_data2)
  p.val <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  
  #p1=ggsurvplot(fit,
  #              pval = TRUE, 
  #              conf.int = FALSE,
  #              risk.table = TRUE, # Add risk table
  #              risk.table.col = "strata", # Change risk table color by groups
  #              linetype = "strata", # Change line type by groups
  #              surv.median.line = "hv", # Specify median survival
  #              #ggtheme = theme_bw(), # Change ggplot2 theme
  #              palette = c("#E7B800", "#2E9FDF"))
  
  p=ggsurvplot(
    fit,                     # survfit object with calculated statistics.
    pval = TRUE,             # show p-value of log-rank test.
    conf.int = FALSE,         # show confidence intervals for 
    # point estimaes of survival curves.
    conf.int.style = "step",  # customize style of confidence intervals
    xlab = "Time in days",   # customize X axis label.
    break.time.by = 200,     # break X axis in time intervals by 200.
    #ggtheme = theme_light(), # customize plot and risk table with a theme.
    ggtheme = theme(text = element_text(family="serif",size = 10),
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
  
  surname <- colnames(sur_data2)[length(deg.cer)+1+sur1]
  
  assign(paste0("survial",surname),p)
}



res <-arrange_ggsurvplots(list(get(paste0("survial",colnames(sur_data2)[length(deg.cer)+2])),
                               get(paste0("survial",colnames(sur_data2)[length(deg.cer)+3])),
                               get(paste0("survial",colnames(sur_data2)[length(deg.cer)+4])),
                               get(paste0("survial",colnames(sur_data2)[length(deg.cer)+5]))))
ggsave("ssurvival.pdf", res,width =16,height=6)


   