################################################################################
#The raw data covered in this chuank are available on request from the authors.
################################################################################

source("Chunk01-R script for RNA-seq data simulation based on improved compcodeR.R")

save_pdf <- function(x, filename, width=6, height=4) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}



##################################################################

#01 t-SNE

#################################################################


M_matrix <-openxlsx::read.xlsx("POS1.xlsx",rowNames=TRUE) # raw data for positive ion

M_matrix2 <-openxlsx::read.xlsx("NEG1.xlsx",rowNames=TRUE) # raw data for negative ion

M_matrix$tet24<-as.numeric(M_matrix$tet24)

M_matrix2$tet24<-as.numeric(M_matrix2$tet24)

iris_unique1<-unique(M_matrix)
iris_unique2<-unique(M_matrix2)


m<-t(iris_unique2)

sam.lab<-substr(rownames(m),1,3)

sam.lab1<-c(rep("Con",22),
            str_replace_all(substr(rownames(m)[23:46],1,3),
                            "Con","NGGCT"))

c.p <- which(sam.lab1%in%"Con")
e.p <- which(sam.lab1%in%"tet")

m<-m[c(e.p,c.p),]

sam.lab<-substr(rownames(m),1,3)

set.seed(5)
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

##################################################################

#01 PCA

#################################################################

iris.pca <- PCA(log2(m+1), graph = FALSE)

ind.p <-fviz_pca_ind(iris.pca,
                     # show points only (nbut not "text") ????????????????????????????????????????????????
                     geom.ind = "point",
                     # ??????????????????
                     col.ind = sam.lab,
                     # ????????????
                     palette = "simpsons",
                     # ???????????? Concentration ellipses
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

##################################################################

#03 Metabee

#################################################################

rawPOS1 <-openxlsx::read.xlsx("removalPOS.xlsx") #POS value of train
rawNEG1 <-openxlsx::read.xlsx("removalNEG.xlsx")#NEG value of train 

rawPOS1$tet24<-as.numeric(rawPOS1$tet24)

rawNEG1$tet24<-as.numeric(rawNEG1$tet24)

POStest<-openxlsx::read.xlsx("POStest.xlsx")#POS value of test
NEGtest<-openxlsx::read.xlsx("NEGtest.xlsx")#POS value of test

POStest <- na.omit(POStest)
NEGtest <- na.omit(NEGtest)

comm_POS<- intersect(POStest$Name,rawPOS1$description)
comm_NEG<- intersect(NEGtest$Name,rawNEG1$description)


rawPOS<-rawPOS1[which(rawPOS1$description%in%comm_POS),]
rawNEG<-rawNEG1[which(rawNEG1$description%in%comm_NEG),]


NEGananotaion<-openxlsx::read.xlsx("NEGananotaion.xlsx") #ananotion information
POSannotation<-openxlsx::read.xlsx("POSannotation.xlsx")

Manao<-rbind(POSannotation,NEGananotaion)

Man.pos<-which(duplicated(Manao$description)==FALSE) 

Manao<-Manao[Man.pos,]

Clinicdt <-openxlsx::read.xlsx("Clinic.xlsx")

sam.lab1<-c(rep("Con",22),
            str_replace_all(substr(colnames(rawPOS)[25:48],1,3),
                            "Con","NGGCT"))

colnames(rawPOS)
#all samples 

M<-c("rawPOS","rawNEG")#POS value of train
M2<-c("rawPOS1","rawNEG1")


for (G in 1:2) {
  
  
  ################
  rawmatrix<-get(M[G])#POS value of train
  
  M_matrix<-rawmatrix[,2:ncol(rawmatrix)]
  
  POS.NA<-which(is.na(M_matrix$description)==FALSE)
  
  M_matrix<-M_matrix[POS.NA,]
  
  rownames(M_matrix)<-M_matrix$description
  
  M_matrix<-M_matrix[,-1]
  
  c.d <- which(sam.lab1 == "Con")#delet the NGGCT
  e.d <- which(sam.lab1 == "tet")
  M_matrix<-M_matrix[,c(e.d,c.d)]
  
  dim(M_matrix)
  colnames(M_matrix)
  ######
  rawmatrix<-get(M2[G])
  
  rownames(rawmatrix)<-rawmatrix$description
  
  rawmatrix<-rawmatrix[,-c(1:2)]
  M_matrix1<-rawmatrix[,c(e.d,c.d)]
  
  sam.lab <- substr(colnames(M_matrix1),1,3)#delet the NGGCT
  
  sam.lab<-as.factor(sam.lab)
  
  M_matrix1<-log2(M_matrix1+0.1)
  
  M_matrix1 <-as.data.frame(t(M_matrix1))
  dim(M_matrix1)
  input1<- cbind(sam.lab, M_matrix1)
  
  input1 <- as.data.frame(input1)
  
  ###############OPLSDA
  
  set.seed(1)
  
  oplsda = opls(input1[,-1], sam.lab, predI = 1, orthoI = NA, crossvalI=5)###fist step for batch1
  
  sample.score <- oplsda@scoreMN %>%as.data.frame() %>% mutate(label = sam.lab,
                                                               o1 = oplsda@orthoScoreMN[,1])
  
  vip <- getVipVn(oplsda)
  
  vip <- as.data.frame(vip)
  
  vip$FeatureName<-rownames(vip)
  
  
  #########FC###second step for common metablites on batch1 and batch2 
  
  e.p <- which(substr(colnames(M_matrix),1,3)%in%"tet") #get a new label
  c.p <- which(substr(colnames(M_matrix),1,3)%in%"Con")
  
  FC <-apply(M_matrix, 1, function(x) {mean(x[e.p])/mean(x[c.p])})
  
  FC <-as.data.frame(FC)
  
  FC$"log2(FC)"<-log2(FC$FC)
  
  FC$FeatureName<-rownames(FC)
  
  
  ##########
  
  sam.lab <- substr(colnames(M_matrix),1,3)#delet the NGGCT
  
  sam.lab<-as.factor(sam.lab)
  
  M_matrix<-log2(M_matrix+0.1)
  
  M_matrix <-as.data.frame(t(M_matrix))
  
  input<- cbind(sam.lab, M_matrix)
  
  input <- as.data.frame(input)
  
  
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
  
  fdrp=p.adjust(deg.count1$Pvalue, "BH")
  fdrp <-data.frame(FeatureName=deg.count1$FeatureName,
                    Pvalue=deg.count1$Pvalue,
                    adjP=fdrp)
  
  vip1 <- vip[fdrp$FeatureName,]
  
  P1 <- ggplot(data=vip1)+
    geom_histogram(aes(x=vip),breaks=seq(0,3,0.05),binwidth = 0.1,
                   fill="steelblue")+
    geom_vline(aes(xintercept=0.9),colour="#BB0000",linetype="dashed")+
    theme_bw()+
    # geom_jitter(shape=16, position = position_jitter(0.2))+
    #scale_y_continuous(limits=c(0, 60), breaks = seq(5,55,5))+
    theme(text = element_text(family="Times New Roman",size = 10),#"Times New Roman","serif"
          panel.grid = element_blank(),
          legend.position = "none",
          # axis.ticks.x.bottom = element_blank(),
          axis.title.x = element_text(size=10),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10)
    )
  P2 <- ggplot(data=FC)+
    geom_histogram(aes(x=FC),breaks=seq(0,3,0.05),binwidth = 0.1,
                   fill="steelblue")+
    geom_vline(aes(xintercept=0.6),colour="#BB0000",linetype="dashed")+
    geom_vline(aes(xintercept=1.3),colour="#BB0000",linetype="dashed")+
    theme_bw()+
    # geom_jitter(shape=16, position = position_jitter(0.2))+
    #scale_y_continuous(limits=c(0, 60), breaks = seq(5,55,5))+
    theme(text = element_text(family="Times New Roman",size = 10),#"Times New Roman","serif"
          panel.grid = element_blank(),
          legend.position = "none",
          # axis.ticks.x.bottom = element_blank(),
          axis.title.x = element_text(size=10),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10)
    )
  
  P3 <- ggplot(data=fdrp)+
    geom_histogram(aes(x=Pvalue),breaks=seq(0,1,0.001),binwidth = 0.2,
                   fill="steelblue")+
    geom_vline(aes(xintercept=0.05),colour="#BB0000",linetype="dashed")+
    theme_bw()+
    # geom_jitter(shape=16, position = position_jitter(0.2))+
    #scale_y_continuous(limits=c(0, 60), breaks = seq(5,55,5))+
    theme(text = element_text(family="Times New Roman",size = 10),#"Times New Roman","serif"
          panel.grid = element_blank(),
          legend.position = "none",
          # axis.ticks.x.bottom = element_blank(),
          axis.title.x = element_text(size=10),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10)
    )
  P4 <- ggplot(data=fdrp)+
    geom_histogram(aes(x=adjP),breaks=seq(0,1,0.001),binwidth = 0.2,
                   fill="steelblue")+
    geom_vline(aes(xintercept=0.15),colour="#BB0000",linetype="dashed")+
    theme_bw()+
    # geom_jitter(shape=16, position = position_jitter(0.2))+
    #scale_y_continuous(limits=c(0, 60), breaks = seq(5,55,5))+
    theme(text = element_text(family="Times New Roman",size = 10),#"Times New Roman","serif"
          panel.grid = element_blank(),
          legend.position = "none",
          # axis.ticks.x.bottom = element_blank(),
          axis.title.x = element_text(size=10),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10)
    )
  
  P <-  ggarrange(P1,P2,P3,P4, ncol=2,nrow = 2)
  # c( "BH","fdr")
  #################
  
  top.features1<-top.features[1:100,]
  
  top.features1<-merge(FC,top.features1)
  
  top.features1<-merge(top.features1,fdrp)
  
  sumdata<-merge(top.features1,vip)
  
  sumdata<-sumdata[sumdata$vip>0.9 & sumdata$Pvalue < 0.05,]
  #
  sumdata<-sumdata[sumdata$FC>1.3 | sumdata$FC<0.6,]
  
  assign(paste0("DEM.",M[G]),sumdata)
  assign(paste0("P.",M[G]),P)
}

DEM1<-rbind(DEM.rawNEG,DEM.rawPOS)

DEM1 <- DEM1 %>% 
  mutate(
    Ion=c(rep("-",nrow(DEM.rawNEG)),rep("+",nrow(DEM.rawPOS)))
  )


hDEM<-merge(DEM1,Manao,by.x = "FeatureName",by.y = "description")

hDEM <- hDEM[hDEM$adjP < 0.15,]

table(hDEM$Ion)#########DMs were identified

#openxlsx::write.xlsx(hDEM,"DEMannotaion.xlsx")

########51 DMs information

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

rownames(hDEG) <- hDEG$Comp

hDEG1<-hDEG[,(which(colnames(hDEG)=="name")+1):(which(colnames(hDEG)=="name")+46)]

#extracting abandance of DMs

dim(hDEG1)

colnames(hDEG1) #check the colnames
rownames(hDEG1)

hDEG1<-hDEG1[,c(which(sam.lab1 == "Con"),which(sam.lab1 == "tet"))]#####delet the NGGCT

hDEG1<-log2(hDEG1)

table(substr(colnames(hDEG1),1,3))

annorow<-hDEG[,c("SuperClass","Class")]

rownames(annorow)<-rownames(hDEG1)

annocol<-as.data.frame(substr(colnames(hDEG1),1,3))###sam.lab1

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
  
  
  
  # if(only.color){
  #   anno_res<-annocolor
  # }else{anno_res<-list(annocol=annocol,
  #                      annorow=annorow,
  #                      annocolor=annocolor)}
  # return(anno_res)
  
  
  
}

col_colors<-get_anno_for_heatmap2(annorow,annocol)



source( "dist_n.R")

###  ????????????
library(d3heatmap)
library(philentropy)

meth<-getDistMethods()
2
n<-2

for(i in 1:46){
  
  m<-meth[i]
  
  a<-hclust(dist_n(t(hDEG1),mtd=m))#wavehedges 
  
  df2<-dendro_data(a,type="rectangle") #"rectangle", "triangle"
  
  Group<- substr(df2$labels$label,1,3)
  
  table(Group)
  
  p<- ggplot(segment(df2))+
    geom_text(data=df2$labels,aes(x=x,y=y,label=label,color=Group),
              angle=90,hjust=1,vjust=0.3,size=3)+
    geom_segment(aes(x=x,y=y,xend=xend,yend=yend))+
    theme(panel.background = element_rect(fill = "white"))
  
  ggsave(paste(meth[i],".pdf"),p,width = 10, height = 6 )
  
}
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
                   clustering_method = 'ward'
)

#"ward", "single", "complete", "average", "mcquitty", "median" or "centroid".
#########################pie
annorow <- na.omit(annorow)
table(annorow$Class)
df <- group_by(annorow, Class)%>%
  dplyr::summarize(percent = n() / nrow(annorow)) %>%
  arrange(desc(percent))%>%mutate(col=col_colors$Class[Class])
###########################
pie(df$percent,labels = df$Class)

pie(df$percent, col = df$col,
    labels = with(df, paste0(round(percent, 2) * 100, "%")))


df1 <- group_by(annorow,SuperClass )%>%
  dplyr::summarize(percent = n() / nrow(annorow)) %>%
  arrange(desc(percent))%>%mutate(col=col_colors$SuperClass[SuperClass])


pie(df1$percent,labels = df1$SuperClass)

pie(df1$percent, col = df1$col,
    labels = with(df1, paste0(round(percent, 2) * 100, "%")))

######################Fig 3A and B
input21 <- readRDS("input21.rds")
HCG <- readRDS("Proteins_Training.rds")
DM<- readRDS("DMs_Training.rds")
NGGCT_DM_Training<- readRDS("NGGCT_DM_Training.rds")
NGGCT_HCG_Training<- readRDS("NGGCT_Proteins_Training.rds")

NGGCT_DM <- NGGCT_DM_Training[,which(colnames(NGGCT_DM_Training)%in%c("Con9","Con15","Con19","Con20"))]

NGGCT_HCG <- NGGCT_HCG_Training[,which(colnames(NGGCT_HCG_Training)%in%c("Con9","Con15","Con19","Con20"))]

colnames(NGGCT_DM) <- paste0(rep("NGGCT",4),1:4)
colnames(NGGCT_HCG) <- paste0(rep("NGGCT",4),1:4)

DM <- DM[rownames(NGGCT_DM),]
DM <- cbind(DM,NGGCT_DM)
HCG <- cbind(HCG,NGGCT_HCG)
colnames(DM)
colnames(HCG)
df <- rbind(DM,HCG)
colnames(df)
dem21 <- colnames(input21[,-1])
dem20 <- dem21[-which(dem21%in%"β-HCG")]
df1 <- df[c(dem20,rownames(HCG)),]
df1 <- log(df1+1)
df1$name <- rownames(df1)
df2<- reshape2::melt(df1,
                     id.vars="name",
                     variable.name = "sample", 
                     value.name = "value")

df2$sample <- substr(df2$sample,1,3)

row.pos <- which(df2$name%in%"N1-Methyl-2-pyridone-5-carboxamide-")

df2$name[row.pos] <- "PY2-"

row.pos <- which(df2$name%in%c("N1-Methyl-2-pyridone-5-carboxamide+"))

df2$name[row.pos] <- "PY2+"

row.pos <- which(df2$name%in%c("N-Acetylaspartylglutamate (NAAG)-"))

df2$name[row.pos] <- "NAAG-"
row.pos <- which(df2$name%in%c("N-Acetylaspartylglutamate (NAAG)+"))

df2$name[row.pos] <- "NAAG+"

df3 <- df2 %>% 
  dplyr::mutate(
    namelength=str_length(as.character(df2$name))
  )

df3 <- df3 [order(df3$namelength),]

library(ggpubr)

sig<-compare_means(value~sample, data=df3,group.by = "name",method = "t.test")

ggplot(data=df3, aes(x= reorder(name, namelength),y= value))+
  geom_boxplot(aes(color=sample), size=0.5, width=0.3,outlier.size = 0.1)+
  scale_color_manual(values = c("#336699","#FCD042","#D45151"))+
  #facet_wrap(~ batch, scales="free",ncol =2)+
  theme_bw()+
  ylab("Abundance (log)")+
  scale_fill_discrete(breaks=c("Con","NGG","tet"),
                      labels=c("Non-iGCT","NGGCT","Germinoma"))+
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
############################

df3 <- df1[which(df1$name%in%dem20),]
col.pos <- which(substr(colnames(df3),1,3)%in%"tet")
df3 <- df3[,col.pos]
df3 <- as.data.frame(t(df3))
loc.label<-c("Suprasellar",
             "Multiple regions",
             "Multiple regions",
             "Multiple regions",
             "Multiple regions",
             "Suprasellar",
             "Suprasellar",
             "Pineal",
             "Pineal",
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

df3$Group<-loc.label

sam_loc<-df3[-which(df3$Group%in%"T"),]


sam_loc<-reshape2::melt(sam_loc,
                        id.vars="Group",
                        variable.name = "Name", 
                        value.name = "Abundance")


sam_loc$Abundance<-as.numeric(sam_loc$Abundance)

sig<-compare_means(Abundance~Group, data=sam_loc,group.by = "Name",method = "t.test")#
View(sig)
sam_loc <- sam_loc %>% 
  dplyr::mutate(
    namelength=str_length(as.character(sam_loc$Name))
  )

sam_loc <- sam_loc [order(sam_loc$namelength),]

ggplot(data=sam_loc, aes(x=reorder(Name, namelength),y=Abundance))+
  geom_boxplot(aes(color=Group),size=0.5, outlier.size = 0.1)+
  scale_color_manual(values = c("#336699","#D15C62","#CBBB7F"))+
  stat_compare_means(aes(group = Group),
                     method = "t.test",
                     label = "p.signif",
                     hide.ns = TRUE,
                     label.x.npc = "right")+
  #facet_wrap(~ Name, scales="free",nrow =5)+
  theme_bw()+
  ylab("Abundance (log)")+
  theme(text = element_text(family="Times New Roman",size = 6), #"Times New Roman","serif"
        panel.grid = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=6,angle = 45),
        axis.text.y = element_text(size=6),
        legend.position ="top",
        legend.text=element_text(size=6),
        legend.title = element_text(size=6)
  )







##################################################################

#04 Methylation data analysis

#################################################################


################
require(GEOquery)   
require(Biobase)
GSE70787 <- getGEO("GSE70787",destdir = './')
#GSE70787

GSE70783 <- exprs(GSE70787[[1]])
GSE70782 <- exprs(GSE70787[[2]])

##########GSE70783
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(impute)
library(wateRmelon)
library(gplots)
library(cluster)
library(stringr)


group1<-GSE70787$`GSE70787-GPL13534_series_matrix.txt.gz`$title

label1<-GSE70787$`GSE70787-GPL13534_series_matrix.txt.gz`$characteristics_ch1

GPL13534 <- read.table("GPL13534-11288.txt",sep = "\t",header = TRUE)

anno34 <- GPL13534[,c("ID","UCSC_RefGene_Name")]

anno34$UCSC_RefGene_Name<-sapply(str_split(anno34$UCSC_RefGene_Name,";"), "[",1)


GSE70783 <- as.data.frame(GSE70783)

GSE70783$ID<-rownames(GSE70783)

GSE70783_1 <- merge(GSE70783,anno34,by="ID")

GSE70783_1<-GSE70783_1[-which(GSE70783_1$UCSC_RefGene_Name%in%""),]

GSE70783_1 <- GSE70783_1[,-1]

GSE70783_1 <- aggregate(.~GSE70783_1$UCSC_RefGene_Name,data = GSE70783_1[,1:c(ncol(GSE70783_1)-1)], mean)

colnames(GSE70783_1)[1]<-"ID"

gp1<-as.data.frame(group1)

label1<-gsub("tissue: Germinoma","tet",label1)

label1[which(label1!="tet")]<-"Con"
library(dplyr)
gp1 <- gp1 %>% 
  mutate(
    group=stringr::str_extract(group1,"GCT"),
    label=label1
  )

rownames(gp1)<-colnames(GSE70783_1)[2:ncol(GSE70783_1)]

View(GSE70783_1)

GSE70783_1 <- GSE70783_1[,c(1,which(gp1$group=="GCT")+1)]

gp1<-na.omit(gp1)

label1 <- gp1$label

###########GSE70782


group2<-GSE70787$`GSE70787-GPL17052_series_matrix.txt.gz`$title

label2<-GSE70787$`GSE70787-GPL17052_series_matrix.txt.gz`$characteristics_ch1

GPL17052 <- read.table("GPL17052-12905.txt",sep = "\t",header = TRUE,fill = T)

anno52 <- GPL17052[,c("ID","GENE_SYMBOL")]

anno52$GENE_SYMBOL<-sapply(str_split(anno52$GENE_SYMBOL,";"), "[",1)

GSE70782 <- as.data.frame(GSE70782)

GSE70782$ID<-rownames(GSE70782)

GSE70782_1<-merge(GSE70782,anno52,by="ID")

GSE70782_1<-GSE70782_1[-which(GSE70782_1$GENE_SYMBOL%in%""),]

GSE70782_1 <- GSE70782_1[,-1]

GSE70782_1 <- aggregate(.~GSE70782_1$GENE_SYMBOL, data = GSE70782_1[,1:c(ncol(GSE70782_1)-1)], mean)

colnames(GSE70782_1)[1]<-"ID"

group2<-group2[1:78]
label2<-label2[1:78]

label2 <- gsub("tissue: Germinoma","tet",label2)

label2[which(label2!="tet")]<-"Con"

GSE70782_1 <- GSE70782_1[,1:79]

##########################


beta.m<-merge(GSE70783_1,GSE70782_1,by="ID")

dim(beta.m)

rownames(beta.m)<-beta.m$ID

beta.m<-beta.m[,-1]

beta.m <- as.matrix(beta.m)

beta.m1 <- impute.knn(beta.m)#NA

beta.m1 <- beta.m1$data +0.000001

beta.m1 <- beta.m1[rowMeans(beta.m1)>0.005,]

beta.m2 <- betaqn(beta.m1)

boxplot(beta.m2,xaxt="n",outlier=F)

group_all<-c(label1,label2)

pdf(file = "densityBeanPlot.pdf")

par(oma=c(2,10,2,2))

densityBeanPlot(beta.m2,sampGroups = group_all, sampNames = colnames(beta.m2))

grest <- makeGenomicRatioSetFromMatrix(beta.m2,what="Beta")


#########

library(FEM)

deg.limma <- read.csv("GSE19348/GSE19348_DEG.csv",row.names = 1)
label1

GSE<- GSE70783[,which(gp1$group=="GCT")]

anntrans <- ann[,11:12]

rownames(ann.mean.eset) <- ann.mean.eset$ENTREZ_GENE_ID

anntrans <- anntrans[which(anntrans$Gene.Symbol%in%rownames(ann.mean.eset)),]

anntrans <- unique(anntrans)

anntrans$ENTREZ_GENE_ID <- as.character(anntrans$ENTREZ_GENE_ID)

Na.DEG.pos <- which(anntrans$ENTREZ_GENE_ID %in%"")

anntrans <- anntrans[-Na.DEG.pos,]

colnames(ann.mean.eset)[1] <- "Gene.Symbol"

ann.mean.eset <- merge(ann.mean.eset,anntrans,by="Gene.Symbol")

rownames(ann.mean.eset)<- ann.mean.eset$Gene.Symbol

DEl.col <- c("Gene.Symbol","ENTREZ_GENE_ID")

DEl.pos <- which(colnames(ann.mean.eset)%in%DEl.col)

DEGeset <- ann.mean.eset[,-DEl.pos]# RNA-seq

label1 <-  str_replace_all(label1, "Con","control")

label1 <-  str_replace_all(label1, "tet","case")

GenStatM.YHC <- function (dnaM.m, pheno.v, chiptype = "450k") {
  if (chiptype == "450k") {
    data("probe450kfemanno")
    probefemanno <- probe450kfemanno
  }
  else if (chiptype == "EPIC") {
    data("probeEPICfemanno")
    probefemanno <- probeEPICfemanno
  }
  else {
    print("ERROR: Please indicate correct data type!")
    break
  }
  extractFn <- function(tmp.v, ext.idx) {
    return(tmp.v[ext.idx])
  }
  map.idx <- match(rownames(dnaM.m), probefemanno$probeID)
  probeInfo.lv <- lapply(probefemanno, extractFn, map.idx)
  beta.lm <- list()
  for (g in 1:6) {
    group.idx <- which(probeInfo.lv[[3]] == g)
    tmp.m <- dnaM.m[group.idx, ]
    tmp.m$ID <- probeInfo.lv$eid[group.idx]
    tmp.m<- aggregate(.~tmp.m$ID, data = tmp.m[,1:c(ncol(tmp.m)-1)], mean)
    rownames(tmp.m) <- tmp.m$ID
    tmp.m <- tmp.m[,-1]
    sel.idx <- which(is.na(rownames(tmp.m)) == FALSE)
    tmp.m <- tmp.m[sel.idx, ]
    nL <- length(factor(rownames(tmp.m)))
    nspg.v <- summary(factor(rownames(tmp.m)), maxsum = nL)
    beta.lm[[g]] <- rowsum(tmp.m, group = rownames(tmp.m))/nspg.v
    print(paste("Done for regional gene group ", g, sep = ""))
  }
  unqEID.v <- unique(c(rownames(beta.lm[[2]]), rownames(beta.lm[[4]]), 
                       rownames(beta.lm[[1]])))
  avbeta.m <- matrix(nrow = length(unqEID.v), ncol = ncol(dnaM.m))
  colnames(avbeta.m) <- colnames(dnaM.m)
  rownames(avbeta.m) <- unqEID.v
  avbeta.m <- as.data.frame(avbeta.m)
  for (gr in c(1, 4, 2)) {
    avbeta.m[match(rownames(beta.lm[[gr]]), rownames(avbeta.m)), ] <- beta.lm[[gr]]
  }
  avbeta.m <- as.matrix(avbeta.m)
  data.m <- avbeta.m
  sampletype.f <- as.factor(pheno.v)
  design.sample <- model.matrix(~0 + sampletype.f)
  colnames(design.sample) <- levels(sampletype.f)
  sampletypes.v <- levels(sampletype.f)
  lmf.o <- lmFit(data.m, design.sample)
  ntypes <- length(levels(sampletype.f))
  ncomp <- 0.5 * ntypes * (ntypes - 1)
  cont.m <- matrix(0, nrow = ncol(design.sample), ncol = ncomp)
  tmp.v <- vector()
  c <- 1
  for (i1 in 1:(ntypes - 1)) {
    for (i2 in (i1 + 1):ntypes) {
      cont.m[i1, c] <- -1
      cont.m[i2, c] <- 1
      tmp.v[c] <- paste(sampletypes.v[i2], "--", sampletypes.v[i1], 
                        sep = "")
      c <- c + 1
    }
  }
  rownames(cont.m) <- sampletypes.v
  colnames(cont.m) <- tmp.v
  lmf2.o <- contrasts.fit(lmf.o, cont.m)
  bay.o <- eBayes(lmf2.o)
  top.lm <- list()
  for (c in 1:ncol(cont.m)) {
    top.lm[[c]] <- topTable(bay.o, coef = c, adjust.method = "fdr", 
                            number = nrow(data.m))
  }
  return(list(top = top.lm, cont = cont.m, avbeta = avbeta.m))
}

statM.o <- GenStatM.YHC(GSE,label1,chiptype = "450k")

View(statM.o)
statM.o$top[[1]]# difference methylation level


statR.o=GenStatR(DEGeset,group)

load("hprdAsigH-13Jun12.Rd")

re=DoIntFEM450k(statM.o,statR.o,hprdAsigH.m,1,1,"avbeta")#the 1 present the 1st of statM.o and statR.o 

Mlistre <- data.frame(re$statM)

fdrp=p.adjust(Mlistre$P.Value, "BH")

fdrp <-data.frame(Mlistre,
                  adjP=fdrp)

DoFEMbi.YHC <- function (intFEM.o, nseeds = 100, gamma = 0.5, nMC = 1000, sizeR.v = c(1, 
                                                                                      100), minsizeOUT = 10, writeOUT = TRUE, nameSTUDY = "X", 
                         ew.v = NULL,pvalue_pvN=0.1) {
  PasteVector <- function(v) {
    vt <- v[1]
    if (length(v) > 1) {
      for (g in 2:length(v)) {
        vt <- paste(vt, v[g], sep = " ")
      }
    }
    vt <- paste(vt, " EnD", sep = "")
    out.v <- sub(" EnD", "", vt)
    out.v <- sub("NA , ", "", out.v)
    out.v <- sub(" , NA", "", out.v)
    out.v <- sub(" , NA , ", " , ", out.v)
    return(out.v)
  }
  Heaviside <- function(v) {
    out.v <- v
    out.v[which(v >= 0)] <- 1
    out.v[which(v < 0)] <- 0
    return(out.v)
  }
  WriteOutPval <- function(pv.v, round.min = 3, round.max = 50) {
    round.v <- round.min:round.max
    th.v <- 10^(-round.v)
    outpv.v <- vector(length = length(pv.v))
    done.idx <- which(pv.v >= th.v[1])
    outpv.v[done.idx] <- round(pv.v[done.idx], round.min)
    todo.idx <- setdiff(1:length(pv.v), done.idx)
    for (i in todo.idx) {
      if (length(which(th.v <= pv.v[i])) > 0) {
        outpv.v[i] <- round(pv.v[i], round.v[min(which(th.v <= 
                                                         pv.v[i]))])
      }
      else {
        outpv.v[i] <- 0
      }
    }
    return(outpv.v)
  }
  adj.m <- intFEM.o$adj
  statM.m <- intFEM.o$statM
  statR.m <- intFEM.o$statR
  if (!identical(rownames(statM.m), rownames(adj.m))) {
    print("Please ensure that rownames of statM.m and adj.m are identical")
  }
  if (!identical(rownames(statR.m), rownames(adj.m))) {
    print("Please ensure that rownames of statR.m and adj.m are identical")
  }
  nameSTUDY <- paste("Both-", nameSTUDY, sep = "")
  x <- org.Hs.egSYMBOL
  mapped_genes <- mappedkeys(x)
  xx <- as.list(x[mapped_genes])
  mapEIDtoSYM.v <- unlist(xx)
  map.idx <- match(rownames(adj.m), names(mapEIDtoSYM.v))
  anno.m <- cbind(rownames(adj.m), mapEIDtoSYM.v[map.idx])
  colnames(anno.m) <- c("EntrezID", "Symbol")
  statM.v <- statM.m[, 1]
  statR.v <- statR.m[, 1]
  sdD <- sqrt(var(statM.v))
  sdR <- sqrt(var(statR.v))
  statR.v <- statR.v * sdD/sdR
  statI.v <- (Heaviside(statM.v) * Heaviside(-statR.v) + Heaviside(-statM.v) * 
                Heaviside(statR.v)) * abs(statM.v - statR.v)
  statI.v[which(sign(statM.v) == sign(statR.v))] <- 0
  names(statI.v) <- rownames(statM.m)
  ntop <- nseeds
  rankN.s <- sort(statI.v, decreasing = TRUE, index.return = TRUE)
  seedsN.v <- names(statI.v)[rankN.s$ix[1:ntop]]
  print("Constructing weighted network")
  tmpA.m <- adj.m
  gr.o <- graph.adjacency(tmpA.m, mode = "undirected")
  tmpE.m <- get.edgelist(gr.o)
  if (is.null(ew.v)) {
    tmpW.v <- vector(length = nrow(tmpE.m))
    for (e in 1:nrow(tmpE.m)) {
      map.idx <- match(tmpE.m[e, ], rownames(tmpA.m))
      tmpW.v[e] <- mean(statI.v[map.idx])
      print(e)
    }
  }else {
    tmpW.v <- ew.v
  }
  minval <- min(setdiff(tmpW.v, 0))
  tmpW.v[tmpW.v == 0] <- minval
  grW.o <- set.edge.attribute(gr.o, "weight", value = tmpW.v)
  adjW.m <- get.adjacency(grW.o, attr = "weight")
  print("Running Spin-Glass algorithm")
  sizeN.v <- vector()
  sgcN.lo <- list()
  for (v in 1:ntop) {
    sgcN.o <- spinglass.community(graph=gr.o, weights = tmpW.v, 
                                  spins = 25, start.temp = 1, stop.temp = 0.1, cool.fact = 0.99, 
                                  update.rule = c("config"), gamma = gamma, vertex = rankN.s$ix[v])
    sizeN.v[v] <- length(sgcN.o$comm)
    sgcN.lo[[v]] <- sgcN.o
    print(paste("Done for seed ", v, sep = ""))
  }
  names(sizeN.v) <- seedsN.v
  print("Module Sizes=")
  print(sizeN.v)
  modN.v <- vector()
  for (v in 1:ntop) {
    subgr.o <- induced.subgraph(grW.o, sgcN.lo[[v]]$comm)
    modN.v[v] <- mean(get.edge.attribute(subgr.o, name = "weight"))
  }
  names(modN.v) <- seedsN.v
  print("Modularity values=")
  print(modN.v)
  print("Starting Monte Carlo Runs")
  modNmc.m <- matrix(nrow = ntop, ncol = nMC)
  for (m in 1:ntop) {
    subgr.o <- induced.subgraph(gr.o, sgcN.lo[[m]]$comm)
    nN <- sizeN.v[m]
    if ((nN > sizeR.v[1]) && (nN < sizeR.v[2])) {
      tmpEL.m <- get.edgelist(subgr.o)
      for (run in 1:nMC) {
        permN.idx <- sample(1:nrow(tmpA.m), nrow(tmpA.m), 
                            replace = FALSE)
        tmpEW.v <- vector()
        for (e in 1:nrow(tmpEL.m)) {
          map.idx <- match(tmpEL.m[e, ], rownames(tmpA.m)[permN.idx])
          tmpEW.v[e] <- mean(statI.v[map.idx])
        }
        subgrW.o <- set.edge.attribute(subgr.o, "weight", 
                                       value = tmpEW.v)
        modNmc.m[m, run] <- mean(get.edge.attribute(subgrW.o, 
                                                    name = "weight"))
      }
    }
    print(paste("Done for seed/module ", m, sep = ""))
  }
  modNpv.v <- rep(1, ntop)
  for (v in 1:ntop) {
    if ((sizeN.v[v] > sizeR.v[1]) && (sizeN.v[v] < sizeR.v[2])) {
      modNpv.v[v] <- length(which(modNmc.m[v, ] > modN.v[v]))/nMC
    }
  }
  names(modNpv.v) <- seedsN.v
  print(modNpv.v)
  print("Summarising and generating output")
  selpvN.idx <- which(modNpv.v < pvalue_pvN)
  selSize.idx <- which(sizeN.v >= minsizeOUT)
  selMod.idx <- intersect(selpvN.idx, selSize.idx)
  print(selMod.idx)
  print(seedsN.v)
  topmodN.m <- matrix(nrow = length(selMod.idx), ncol = 6)
  map.idx <- match(seedsN.v[selMod.idx], anno.m[, 1])
  seedsSYM.v <- anno.m[map.idx, 2]
  topmodN.m[, 1] <- seedsN.v[selMod.idx]
  topmodN.m[, 2] <- seedsSYM.v
  topmodN.m[, 3:5] <- cbind(sizeN.v[selMod.idx], modN.v[selMod.idx], 
                            modNpv.v[selMod.idx])
  mi <- 1
  for (m in selMod.idx) {
    tmpEID.v <- rownames(tmpA.m)[sgcN.lo[[m]]$comm]
    genes.v <- anno.m[match(tmpEID.v, anno.m[, 1]), 2]
    topmodN.m[mi, 6] <- PasteVector(genes.v)
    mi <- mi + 1
  }
  colnames(topmodN.m) <- c("EntrezID(Seed)", "Symbol(Seed)", 
                           "Size", "Mod", "P", "Genes")
  if (writeOUT) {
    write.table(topmodN.m, file = paste("topFEM-", nameSTUDY, 
                                        ".txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
  }
  seltopmodN.lm <- list()
  for (m in 1:length(selMod.idx)) {
    tmpEID.v <- rownames(tmpA.m)[sgcN.lo[[selMod.idx[m]]]$comm]
    map.idx <- match(tmpEID.v, anno.m[, 1])
    map1.idx <- match(tmpEID.v, rownames(tmpA.m))
    seltopmodN.m <- cbind(anno.m[map.idx, 1:2], statM.m[map1.idx, 
    ], statR.m[map1.idx, ], statI.v[map1.idx])
    seltopmodN.lm[[m]] <- seltopmodN.m
    colnames(seltopmodN.lm[[m]]) <- c("EntrezID", "Symbol", 
                                      "stat(DNAm)", "P(DNAm)", "stat(mRNA)", "P(mRNA)", 
                                      "stat(Int)")
  }
  names(seltopmodN.lm) <- seedsSYM.v
  if (writeOUT) {
    for (m in 1:length(selMod.idx)) {
      out.m <- seltopmodN.lm[[m]]
      out.m[, 3] <- round(as.numeric(out.m[, 3]), 2)
      out.m[, 4] <- WriteOutPval(as.numeric(out.m[, 4]), 
                                 round.max = 100)
      out.m[, 5] <- round(as.numeric(out.m[, 5]), 2)
      out.m[, 6] <- WriteOutPval(as.numeric(out.m[, 6]), 
                                 round.max = 100)
      out.m[, 7] <- round(as.numeric(out.m[, 7]), 2)
      write(paste(seedsSYM.v[m], " (", nrow(seltopmodN.lm[[m]]), 
                  " genes)", sep = ""), file = paste("topFEMLists-", 
                                                     nameSTUDY, ".txt", sep = ""), ncolumns = 1, append = TRUE)
      write.table(out.m, file = paste("topFEMLists-", nameSTUDY, 
                                      ".txt", sep = ""), quote = FALSE, sep = "\t", 
                  row.names = FALSE, append = TRUE)
    }
  }
  return(list(size = sizeN.v, mod = modN.v, pv = modNpv.v, 
              selmod = selMod.idx, fem = topmodN.m, topmod = seltopmodN.lm, 
              sgc = sgcN.lo, ew = tmpW.v, adj = intFEM.o$adj))
}

DoFEMbi=DoFEMbi.YHC(re, nseeds =100, gamma = 0.5, 
                    nMC = 1000, sizeR.v = c(1,120), 
                    minsizeOUT = 3, writeOUT = TRUE, 
                    nameSTUDY = "TEST", ew.v = NULL,pvalue_pvN=0.05) #if error, check the available nseeds


for (n1 in 1:length(names(DoFEMbi$topmod))) {
  
  net_tpm <- DoFEMbi$topmod
  
  net_tpm1 <- net_tpm[[n1]]
  
  FemModShow(net_tpm1,name=names(DoFEMbi$topmod)[n1], DoFEMbi)
  
}

View(DoFEMbi$topmod$ROBO1)
#Visualize the module information



##################################################################

#05 Transcripte level analysis

#################################################################

#DNA microarray

############################

library(affy)
library(annotate)
getGEOSuppFiles("GSE19348") 

setwd("GSE19348/")

untar("GSE19348_RAW.tar")  
files <- dir(pattern="gz$")  

data <- ReadAffy(filenames=files)
affydb<-annPkgName(data@annotation,type="db")
require(affydb, character.only=TRUE)
eset<-rma(data,verbose=FALSE)
eset <- exprs(eset)  

eset<-eset[which(apply(eset, 1, min) >0.000001),]

tmp1 <- normalize.quantiles(eset)

rownames(tmp1) <- rownames(eset)

colnames(tmp1) <- colnames(eset)

eset <- tmp1

mypal1 <- terrain.colors(ncol(eset))
boxplot(eset, col = mypal1, main = "DNA microarray data")


ann <- read.csv("GPL570.csv") #annotation of GPL570
#download from GEO

probe.list <- row.names(eset)
anndataframe <- NULL
for (prob in probe.list) {
  
  Gene.Symbol <- as.character(ann$Gene.Symbol[which(ann$ID == prob )])
  rowone <- data.frame("Symbol"= Gene.Symbol,"Probe"=prob)
  anndataframe <- rbind(anndataframe,rowone )
  
}

 
removeops.prob <- grep("///",anndataframe$Symbol)

anndataframe <- anndataframe[-removeops.prob, ]

 
ann.eset <- eset[anndataframe$Probe, ]

ann.eset <-  cbind(rownames(ann.eset),ann.eset)

ann.eset <- cbind(anndataframe,ann.eset)

a <- apply(ann.eset,1,as.character)
a <- apply(a,1,as.numeric)
colnames(a) <-colnames(ann.eset)
a <-as.data.frame(a)
class(a[10,10])
a$Symbol<-ann.eset$Symbol           

ann.mean.eset <- aggregate(.~a$Symbol,data = a[  ,-c(1:3)], mean)

rownames(ann.mean.eset) <- ann.mean.eset$`a$Symbol`

ann.mean.eset <- ann.mean.eset[,-1]

colnames(ann.mean.eset)
ann.mean.eset <- ann.mean.eset[,c(1:5,8:13)]

group<-c(rep("Germinoma",5),rep("Control",6))

design <- model.matrix( ~ group)

fit <- lmFit(ann.mean.eset, design)

fit2 <- eBayes(fit, trend = FALSE)  

limmaDEGs <- topTable(fit2, coef = 2, number = Inf)

limmaDEGs$logFC <- round(limmaDEGs$logFC,2)

deg.limma <- limmaDEGs[abs(limmaDEGs$logFC) >= 2 & limmaDEGs$adj.P.Val < 0.05, ]

#########################

#RNA-seq

#########################

NGGCT<-openxlsx::read.xlsx("NGGCTCount.xlsx")#raw data

eset <- aggregate(.~NGGCT$symbol,data = NGGCT[,2:ncol(NGGCT)], mean)

rownames(eset) <- eset$`NGGCT$symbol`

eset <- eset[,-1]

eset2<-eset[which(apply(eset, 1, min) >0.000001),]

mypal1 <- terrain.colors(ncol(eset2))

boxplot(log(eset2), col = mypal1)

colnames(eset2)

eset3 <- eset2[,c("Count_.T1", "Count_.T4","Count_.T8", "Count_.T6",
                  "Count_.T2","Count_.T9")]

group <- c("Control",
           "Control",
           "Germinoma",
           "Control",
           "Germinoma",
           "Germinoma")


group <- as.factor(group)

design <- model.matrix( ~ group)

y<-voom(eset3, design) #RNA-seq 

fit <- lmFit(y, design) #RNA-seq 

fit2 <- eBayes(fit)  

limmaDEGs1 <- topTable(fit2, coef = 2, number = Inf)

deg.limma1 <- limmaDEGs[abs(limmaDEGs1$logFC) > 1.5&limmaDEGs1$P.Value < 0.05, ]

deg.limma1 <- deg.limma[order(deg.limma1$adj.P.Val), ]


#######metabolistDEGs

GM_gene <- openxlsx::read.xlsx("metabolistDEGs.xlsx")# 


colnames(limmaDEGs)[4] <- "PValue"


limmaDEGs <- limmaDEGs %>% mutate(
  state=case_when(
    logFC > 2 & PValue <0.05 & !rownames(limmaDEGs)%in%GM_gene$Gene.Symbol~ 'up',
    
    logFC < -2 & PValue <0.05 & !rownames(limmaDEGs)%in%GM_gene$Gene.Symbol ~ 'down',
    
    rownames(limmaDEGs)%in%GM_gene$Gene.Symbol & logFC > 2 ~ 'up1',
    
    rownames(limmaDEGs)%in%GM_gene$Gene.Symbol & logFC < -2 ~ 'down1',
    
    TRUE ~ 'none'
  )
)


library(ggplot2)
library(tidyverse)
library(ggrepel)

ggplot(limmaDEGs,aes(x=logFC,y=-log10(PValue),color=state))+
  geom_point()+
  scale_color_manual(values = c("ins"="grey",
                                "up"="#EDD3D6",
                                "down"="#CFDEF4",
                                "up1"="#990099",
                                "down1"="#990099"))+
  #geom_text_repel(data = GM_gene,
  #                aes(label = Gene.Symbol),
  #                size = 2,segment.color = "black",
  #                max.overlaps = 100)+
  
  theme_classic()+
  ylab('-log10 (P value)')+
  xlab('log(FoldChange)')+
  #ylim(0,5)+
  xlim(-8,8)+
  geom_vline(xintercept=c(-1.5,1.5),
             lty=3,
             col="black",
             lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),
             lty=3,
             col="black",
             lwd=0.5)


###################

#06 GSVA analysis

####################

library("KEGGREST")

keggGetclass<-function(x){
  
  gs<-keggGet(x)
  
  class<-gs[[1]]$CLASS
  
  if(length(class)>0) {
    
    a <- try(unlist(strsplit(class,";")))
    
    b<-a[2]
    
  } else {
    b<-NA}
  return(b)
  
}

keggGetmetagene<-function(x){
  
  gs<-keggGet(x)
  
  COM<-names(gs[[1]]$COMPOUND)
  
  if(length(COM)>0){
    
    Gene<-gs[[1]]$GENE
    
    a<-sapply(str_split(Gene,";"), "[",1)
    
    b<-grep("[A-Z]", a, value = TRUE) #grep the upper 
    
    a <- setdiff(a,b)
    
    a<-na.omit(a)
    
    a<-as.character(a)
    
    return(a)
  }
  
}

keggGetgene<-function(x){
  
  gs<-keggGet(x)
  
  Gene<-gs[[1]]$GENE
  
  a<-sapply(str_split(Gene,";"), "[",1)
  
  b<-grep("[A-Z]", a, value = TRUE) #grep the upper 
  
  a <- setdiff(a,b)
  
  a<-na.omit(a)
  
  a<-as.character(a)
  
  return(a)
  
} #entrz ID

#keggList("organism")

KEGGlist<-keggList("pathway","hsa")

keggnames<-names(KEGGlist)

keggnames<-substring(keggnames,"6")

keggnames_list <- list()

keggnames_list <- c(keggnames_list, keggnames)

KEGGclass <- lapply(keggnames_list, keggGetclass)#obtain the class of metabolism pathway

KEGGgene <- lapply(keggnames_list, keggGetmetagene)

#KEGGCom <- lapply(keggnames_list, keggGetcompound)#obtain the metabolits in metabolism pathway

KEGGlist1<-as.data.frame(KEGGlist)


names(KEGGgene)<-KEGGlist1$KEGGlist
names(KEGGclass)<-KEGGlist1$KEGGlist
KEGGc1=KEGGgene[!sapply(KEGGgene,is.null)]


class<-unlist(KEGGclass[names(KEGGc1)])

class<-as.data.frame(class)

KEGGc2 <- KEGGc1[which(names(KEGGc1)%in%rownames(class))]

gsva_es <- GSVA::gsva(as.matrix(log2(DEGeset+1)), KEGGc2, method= "zscore", kcdf="Gaussian")
#method="gsva","ssgsea", "zscore", "plage".  kcdf=c("Gaussian", "Poisson", "none")

class2<-unlist(KEGGclass[names(KEGGc2)])

class2<-as.data.frame(class2)

rownames(class2) <- sapply(str_split(rownames(class2)," - H"), "[",1)

rownames(gsva_es) <- sapply(str_split(rownames(gsva_es)," - H"), "[",1)

M.index<-intersect(rownames(gsva_es),rownames(class2))

anno<-class2[M.index,]

anno<-as.data.frame(anno)

rownames(anno)<-M.index

annocol<-as.data.frame(group)###sam.lab1

rownames(annocol)<-colnames(DEGeset)

colnames(annocol)<-"Group"
get_anno_for_heatmap2<-function(annocol,annorow=NULL,color=NULL,only.color=F){
  require(plyr)
  require(stringr)
  if(is.null(color)){
    require(RColorBrewer)
    color=c(brewer.pal(12,"Set3"),brewer.pal(12,"Paired"),brewer.pal(8,"Set2"),brewer.pal(9,"Set1"),brewer.pal(8,"Dark2"))
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

col_colors<-get_anno_for_heatmap2(anno,annocol)



pheatmap(gsva_es,
         show_rownames = T,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows=T,
         fontsize_row=10, #??????????????????
         border_color = "NA",#FFFFFF",
         scale = "column",
         angle_col=90, #??????????????????
         annotation_row = anno,
         annotation_col = annocol,
         annotation_colors = col_colors,
         color =colorRampPalette(c("#0225A2", "white","#FF2400"))(100)
         #color =colorRampPalette(c("#5BB93B","#030303", "#C22D1D"))(100)
)


##############

Mlistre 

Mlistre$"ENTREZ_GENE_ID" <- rownames(Mlistre)

GM_heatmap <- merge(GM_gene,fdrp,by="ENTREZ_GENE_ID")


GM_heatmap <- GM_heatmap[abs(GM_heatmap$t.y)>2 & GM_heatmap$adjP<0.05, ]


GM_heatmap <- GM_heatmap[,c(2,3,9)]


rownames(GM_heatmap) <- GM_heatmap$Gene.Symbol

GM_heatmap <- GM_heatmap[,-1]

#load("intergene.RData")

GM_heatmap <- GM_heatmap[intergene,]

GM_heatmap <- na.omit(GM_heatmap)

dim(GM_heatmap)

pheatmap(GM_heatmap,
         show_rownames = T,
         show_colnames = T,
         cluster_cols = F,
         cluster_rows=T,
         fontsize_row=6, #??????????????????
         border_color = "NA",#FFFFFF",
         scale = "none",
         angle_col=90, #??????????????????
         #annotation_row = anno,
         #annotation_col = annocol,
         # annotation_colors = col_colors,
         color =colorRampPalette(c("#0225A2", "white","#FF2400"))(100)
         #color =colorRampPalette(c("#5BB93B","#030303", "#C22D1D"))(100)
)



###################

#06 Cibersort

####################

source("Cibersort.R")
LM22 <- read.table("LM22.txt",sep = "\t",row.names = 1,header = TRUE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggpubr)
eng <- rownames(DEGeset)

SC_GSE52529_A <- cbind(eng,DEGeset)

gene.df <- bitr(eng,
                fromType = "ENTREZID", 
                toType =  "SYMBOL",
                OrgDb= org.Hs.eg.db ) 




SC_GSE52529_A <- SC_GSE52529_A[gene.df$ENTREZID,]

SC_GSE52529_A <- cbind(gene.df,SC_GSE52529_A)

SC_GSE52529_A <- aggregate(.~SC_GSE52529_A$SYMBOL,data=SC_GSE52529_A[,-c(1,2,3)],mean)

rownames(SC_GSE52529_A ) <- SC_GSE52529_A$`SC_GSE52529_A$SYMBOL`

SC_GSE52529_A <- SC_GSE52529_A[,-1]

DEGeset <- SC_GSE52529_A


CBres <- Cibersort(bulkdata = DEGeset,signature = LM22)


library(reshape2)

gd1_long1 <- melt(CBres[c(1:5,8:13),1:22],
                  value.name='Proportion')


colnames(gd1_long1) <- c("Sample_ID","Cell_Type","Proportion")
 
library(ggplot2)

library(RColorBrewer)
library(ggsci)

rhg_cols1 <- pal_igv("default", alpha = 0.66)(22)
 
lable <- unlist(lapply(strsplit(colnames(DEGeset),"[.]"),"[",1))
 

p2 <- ggplot(gd1_long1,aes(x=Sample_ID,y=Proportion,fill=Cell_Type))+
  geom_bar(stat="identity")+
  theme_classic()+
  scale_fill_manual(values = rhg_cols1)+
  scale_x_discrete(labels=lable)+
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))



Experimental <- apply(CBres[1:5,1:22], 2, mean)

Control <- apply(CBres[8:13,1:22], 2, mean)

imm <- data.frame("Germinoma" = Experimental,"Control" = Control )



library(ggpubr)

gd1_long2 <- melt(imm,
                  value.name='Proportion')

colnames(gd1_long2) <- c("Group","Proportion")
gd1_long2$celltype <- rep("Total immune cells",nrow(gd1_long2))
gd1 <- gd1_long2

for(cc1 in 1:22){
  
  cell_t <- colnames(CBres)[cc1]
  
  gd1_long2 <- data.frame("Group"= c(rep("Germinoma",5),rep("Control",6)),
                          "Proportion"=as.numeric(c(CBres[1:5,cc1],CBres[8:13,cc1])),
                          "celltype"=cell_t)
  gd1 <- rbind(gd1,gd1_long2)
}


gd1 <- gd1 %>% 
  dplyr::mutate(
    namelength=str_length(as.character(gd1$celltype))
  )

gd1 <- gd1 [order(gd1$namelength),]


sig<-compare_means(Proportion~Group, data=gd1,group.by = "celltype",method = "wilcox.test")

p<-ggplot(data=gd1, aes(x=reorder(celltype, namelength),y=Proportion))+
  geom_boxplot(aes(color=Group),size=0.8, outlier.size = 0.1)+
  scale_color_manual(values = c("#CBBB7F","#336699"))+
  stat_compare_means(aes(group = Group),
                     method = "wilcox.test",
                     label = "p.signif",
                     label.x.npc = "right")+
  # facet_wrap(~ Name, scales="free")+
  theme_bw()+
  theme(text = element_text(family="serif",size = 8),
        panel.grid = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=6,angle = 90),
        axis.text.y = element_text(size=6)
  )

p

