rm(list=ls())
library(paleomorph)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggpubr)

setwd("C://Users/pc/Desktop/Temp from Archive/")

#Set up the threshold for the Tricontacts (ttc)
My_Threshold <- 400

#Import AfterAPC for WT and dA
NPC_WT_AfterAPC <- read.csv2("R_outpout/NPC_WT/NPC_WT_AfterAPC.csv", stringsAsFactors = F)
NPC_dA_AfterAPC <- read.csv2("R_outpout/NPC_dA//NPC_dA_AfterAPC.csv", stringsAsFactors = F)
NPC_dA_AfterAPC[,2] <- NPC_dA_AfterAPC[,2]*-1
NPC_dA_AfterAPC[,3] <- NPC_dA_AfterAPC[,3]*-1


#Calculate the distances between the points after PCA : AfterAPC

#For NPC_WT
rm(tAfterAPC)
NPC_WT <- rbind(read.csv2("Merge_summary/All_replicates/NPC_wt_rep1_merge.csv", header = T, stringsAsFactors = F),
                read.csv2("Merge_summary/All_replicates/NPC_WT_rep3_merge.csv", header = T, stringsAsFactors = F))
NPC_WT[,-1]<- apply(NPC_WT[,-1],2, as.numeric)

#Calculate distances then remove dist>2.2µm for NPC WT
#dist FC1FC2
NPC_WT$F1F2 <- sqrt(((NPC_WT$xFC2-NPC_WT$xFC1)*40)**2+((NPC_WT$yFC2-NPC_WT$yFC1)*40)**2+((NPC_WT$zFC2-NPC_WT$zFC1)*140)**2)
#dist FC2FC3
NPC_WT$F2F3 <- sqrt(((NPC_WT$xFC3-NPC_WT$xFC2)*40)**2+((NPC_WT$yFC3-NPC_WT$yFC2)*40)**2+((NPC_WT$zFC3-NPC_WT$zFC2)*140)**2)
#dist FC1FC3
NPC_WT$F1F3 <- sqrt(((NPC_WT$xFC3-NPC_WT$xFC1)*40)**2+((NPC_WT$yFC3-NPC_WT$yFC1)*40)**2+((NPC_WT$zFC3-NPC_WT$zFC1)*140)**2)
cut_off_mydist<- function(a){
  a <- a[a$F1F2 %in% subset(a$F1F2, a$F1F2 < 2200),]
  a <- a[a$F2F3 %in% subset(a$F2F3, a$F2F3 < 2200),]
  a <- a[a$F1F3 %in% subset(a$F1F3, a$F1F3 < 2200),]
}
NPC_WT <- cut_off_mydist(NPC_WT)
tAfterAPC <- data.frame(ID=rep(NA,NROW(NPC_WT)),
                        Xprom=rep(NA,NROW(NPC_WT)),
                        Yprom=rep(NA,NROW(NPC_WT)),
                        XEnhA=rep(NA,NROW(NPC_WT)),
                        YEnhA=rep(NA,NROW(NPC_WT)),
                        XEnhb=rep(NA,NROW(NPC_WT)),
                        YEnhb=rep(NA,NROW(NPC_WT)))
u <- 1
for (i in seq(1,NROW(NPC_WT_AfterAPC),by=3)) {
  tAfterAPC[u,1] <- NPC_WT_AfterAPC[i,1]
  tAfterAPC[u,2:3] <- NPC_WT_AfterAPC[i,2:3]
  tAfterAPC[u,4:5] <- NPC_WT_AfterAPC[i+1,2:3]
  tAfterAPC[u,6:7] <- NPC_WT_AfterAPC[i+2,2:3]
  u <- u+1
  
}
#dist Prom EnhA
tAfterAPC$d_Prom_EnhA <- sqrt((tAfterAPC$XEnhA-tAfterAPC$Xprom)**2+(tAfterAPC$YEnhA-tAfterAPC$Yprom)**2)
#dist EnhA EnhB
tAfterAPC$d_EnhA_EnhB <- sqrt((tAfterAPC$XEnhb-tAfterAPC$XEnhA)**2+(tAfterAPC$YEnhb-tAfterAPC$YEnhA)**2)
#dist Prom EnhB
tAfterAPC$d_Prom_Enhb <- sqrt((tAfterAPC$XEnhb-tAfterAPC$Xprom)**2+(tAfterAPC$YEnhb-tAfterAPC$Yprom)**2)


##Collect coordinates x and y for distances shorter then threshold
tAfterAPC_NPC_WT <- tAfterAPC
#First threshold 


NPC_WT_close_distance_Prom_ENhA <- tAfterAPC_NPC_WT[tAfterAPC_NPC_WT$d_Prom_EnhA < My_Threshold,c(2,3,4,5,8)] 
NPC_WT_close_distance_Prom_ENhB <- tAfterAPC_NPC_WT[tAfterAPC_NPC_WT$d_Prom_Enhb < My_Threshold,c(2,3,6,7,10)]
NPC_WT_close_distance_EnhA_Enhb <- tAfterAPC_NPC_WT[tAfterAPC_NPC_WT$d_EnhA_EnhB < My_Threshold,c(4,5,6,7,9)]
NPC_WT_Tricontact <- tAfterAPC_NPC_WT[tAfterAPC_NPC_WT$d_Prom_EnhA< My_Threshold &
                          tAfterAPC_NPC_WT$d_EnhA_EnhB<My_Threshold &
                          tAfterAPC_NPC_WT$d_Prom_Enhb<My_Threshold, 1]
NPC_WT_Tricontact <- NPC_WT_AfterAPC[NPC_WT_AfterAPC$ID %in% NPC_WT_Tricontact,]


#Calculate the distances between the points after PCA : AfterAPC
#For NPC_dA
rm(tAfterAPC)
NPC_dA <- rbind(read.csv2("Merge_summary/All_replicates/NPC_dA_rep1_merge.csv", header = T, stringsAsFactors = F),
                read.csv2("Merge_summary/All_replicates/NPC_dA_rep2_merge.csv", header = T, stringsAsFactors = F),
                read.csv2("Merge_summary/All_replicates/NPC_dA_rep3_merge.csv", header = T, stringsAsFactors = F))
NPC_dA[,-1]<- apply(NPC_dA[,-1],2, as.numeric)
#str(NPC_dA)
#Calculate distances then remove dist>2.2µm for NPC dA
#dist FC1FC2
NPC_dA$F1F2 <- sqrt(((NPC_dA$xFC2-NPC_dA$xFC1)*40)**2+((NPC_dA$yFC2-NPC_dA$yFC1)*40)**2+((NPC_dA$zFC2-NPC_dA$zFC1)*140)**2)
#dist FC2FC3
NPC_dA$F2F3 <- sqrt(((NPC_dA$xFC3-NPC_dA$xFC2)*40)**2+((NPC_dA$yFC3-NPC_dA$yFC2)*40)**2+((NPC_dA$zFC3-NPC_dA$zFC2)*140)**2)
#dist FC1FC3
NPC_dA$F1F3 <- sqrt(((NPC_dA$xFC3-NPC_dA$xFC1)*40)**2+((NPC_dA$yFC3-NPC_dA$yFC1)*40)**2+((NPC_dA$zFC3-NPC_dA$zFC1)*140)**2)
cut_off_mydist<- function(a){
  a <- a[a$F1F2 %in% subset(a$F1F2, a$F1F2 < 2200),]
  a <- a[a$F2F3 %in% subset(a$F2F3, a$F2F3 < 2200),]
  a <- a[a$F1F3 %in% subset(a$F1F3, a$F1F3 < 2200),]
}
NPC_dA <- cut_off_mydist(NPC_dA)
#rm(tAfterAPC)
tAfterAPC <- data.frame(ID=rep(NA,NROW(NPC_dA)),
                        Xprom=rep(NA,NROW(NPC_dA)),
                        Yprom=rep(NA,NROW(NPC_dA)),
                        XEnhA=rep(NA,NROW(NPC_dA)),
                        YEnhA=rep(NA,NROW(NPC_dA)),
                        XEnhb=rep(NA,NROW(NPC_dA)),
                        YEnhb=rep(NA,NROW(NPC_dA)))
u <- 1
for (i in seq(1,NROW(NPC_dA_AfterAPC),by=3)) {
  tAfterAPC[u,1] <- NPC_dA_AfterAPC[i,1]
  tAfterAPC[u,2:3] <- NPC_dA_AfterAPC[i,2:3]
  tAfterAPC[u,4:5] <- NPC_dA_AfterAPC[i+1,2:3]
  tAfterAPC[u,6:7] <- NPC_dA_AfterAPC[i+2,2:3]
  u <- u+1
  
}
#dist Prom EnhA
tAfterAPC$d_Prom_EnhA <- sqrt((tAfterAPC$XEnhA-tAfterAPC$Xprom)**2+(tAfterAPC$YEnhA-tAfterAPC$Yprom)**2)
#dist EnhA EnhB
tAfterAPC$d_EnhA_EnhB <- sqrt((tAfterAPC$XEnhb-tAfterAPC$XEnhA)**2+(tAfterAPC$YEnhb-tAfterAPC$YEnhA)**2)
#dist Prom EnhB
tAfterAPC$d_Prom_Enhb <- sqrt((tAfterAPC$XEnhb-tAfterAPC$Xprom)**2+(tAfterAPC$YEnhb-tAfterAPC$Yprom)**2)
#summary(tAfterAPC[,8:10])

##Collect coordinates x and y for distances shorter then threshold
#First threshold 

NPC_dA_close_distance_Prom_ENhA <- tAfterAPC[tAfterAPC$d_Prom_EnhA < My_Threshold,c(2,3,4,5,8)] 
NPC_dA_close_distance_Prom_ENhB <- tAfterAPC[tAfterAPC$d_Prom_Enhb < My_Threshold,c(2,3,6,7,10)]
NPC_dA_close_distance_EnhA_Enhb <- tAfterAPC[tAfterAPC$d_EnhA_EnhB < My_Threshold,c(4,5,6,7,9)]
NPC_dA_Tricontact <- tAfterAPC[tAfterAPC$d_Prom_EnhA< My_Threshold &
                          tAfterAPC$d_EnhA_EnhB<My_Threshold &
                          tAfterAPC$d_Prom_Enhb<My_Threshold, 1]
NPC_dA_Tricontact <- NPC_dA_AfterAPC[NPC_dA_AfterAPC$ID %in% NPC_dA_Tricontact,]



#Generate the plot with all the aligned points 
#Keep the same number points for both conditions, same as WT
a <- ggplot(NPC_WT_AfterAPC, aes(x=Dim1, y=Dim2))+
  geom_polygon(aes(group=ID), fill=NA, colour="grey51", size=0.01)+
  geom_point(aes(colour=molecular_ID), size=1)+
  #stat_density2d(
  #  data=AfterAPC, 
  #  aes(x=Dim1,y=Dim2, fill=molecular_ID), 
  #  alpha=0.2, geom="polygon", bins=15, show.legend=F)+
  # xlim(-1500,1000)+ylim(-1500,1000)+
  xlab("Distance (nm)")+ylab("Distance (nm)")+
  scale_colour_discrete(type=c("#a1291b","#7bac25","#2570a3"),name="Molecular\nIdentity", labels= c("Enhancer A","Enhancer B","Promoter"))+
  theme_minimal()+
  theme(legend.position = c(0.9,0.2))
a





a_and_dist <- a +  geom_segment(data=NPC_WT_close_distance_Prom_ENhA,
                 aes(x=NPC_WT_close_distance_Prom_ENhA$Xprom, 
                     xend=NPC_WT_close_distance_Prom_ENhA$XEnhA,
                     y= NPC_WT_close_distance_Prom_ENhA$Yprom, 
                     yend=NPC_WT_close_distance_Prom_ENhA$YEnhA),
                 col="#c53bc7")+
  
  geom_segment(data=NPC_WT_close_distance_Prom_ENhB,
               aes(x=NPC_WT_close_distance_Prom_ENhB$Xprom, 
                   xend=NPC_WT_close_distance_Prom_ENhB$XEnhb,
                   y= NPC_WT_close_distance_Prom_ENhB$Yprom, 
                   yend=NPC_WT_close_distance_Prom_ENhB$YEnhb),
               col="#23f7d3")+
  
  geom_segment(data=NPC_WT_close_distance_EnhA_Enhb,
               aes(x=NPC_WT_close_distance_EnhA_Enhb$XEnhA, 
                   xend=NPC_WT_close_distance_EnhA_Enhb$XEnhb,
                   y= NPC_WT_close_distance_EnhA_Enhb$YEnhA, 
                   yend=NPC_WT_close_distance_EnhA_Enhb$YEnhb),
               col="#f6880c")+
  
  geom_polygon(data=NPC_WT_Tricontact, aes(x=Dim1, y=Dim2, group=ID),  fill=NA, colour="red")+
  xlim(-500,500)+ylim(-500,500)+ggtitle("NPC WT")
a_and_dist

#DO the same plot with the same number of n for WT and dA
#I keep WT like that as it s the smallest one
#I took only the first 484 of dA
NPC_dA_AfterAPC_1452 <- NPC_dA_AfterAPC[1:1452,]
tAfterAPC_484 <- tAfterAPC[1:484,]


NPC_dA_close_dist_Prom_ENhA_1452 <- tAfterAPC_484[tAfterAPC_484$d_Prom_EnhA < My_Threshold,c(2,3,4,5,8)] 
NPC_dA_close_dist_Prom_ENhB_1452 <- tAfterAPC_484[tAfterAPC_484$d_Prom_Enhb < My_Threshold,c(2,3,6,7,10)]
NPC_dA_close_dist_EnhA_Enhb_1452 <- tAfterAPC_484[tAfterAPC_484$d_EnhA_EnhB < My_Threshold,c(4,5,6,7,9)]
NPC_dA_Tricontact_1452 <- tAfterAPC_484[tAfterAPC_484$d_Prom_EnhA< My_Threshold &
                                 tAfterAPC_484$d_EnhA_EnhB<My_Threshold &
                                 tAfterAPC_484$d_Prom_Enhb<My_Threshold, 1]
NPC_dA_Tricontact_1452 <- NPC_dA_AfterAPC_1452[NPC_dA_AfterAPC_1452$ID %in% NPC_dA_Tricontact_1452,]



bbis <- ggplot(NPC_dA_AfterAPC_1452, aes(x=Dim1, y=Dim2, colour=molecular_ID))+
  geom_polygon(aes(group=ID), fill=NA, colour="grey51", size=0.01, alpha=0.1)+
  geom_point(size=1)+
  #stat_density2d(
  #  data=AfterAPC, 
  #  aes(x=Dim1,y=Dim2, fill=molecular_ID), 
  #  alpha=0.2, geom="polygon", bins=15, show.legend=F)+
  # xlim(-1500,1000)+ylim(-1500,1000)+
  xlab("Distance (nm)")+ylab("Distance (nm)")+
  scale_colour_discrete(type = c("#a1291b","#7bac25","#2570a3"), name="Molecular\nIdentity", labels= c("Enhancer A","Enhancer B","Promoter"))+
  theme_minimal()+
  theme(legend.position = c(0.9,0.2))
bbis

bbis_and_dist <-  bbis + 
  geom_segment(data=NPC_dA_close_dist_Prom_ENhA_1452,
                 aes(x=NPC_dA_close_dist_Prom_ENhA_1452$Xprom, 
                     xend=NPC_dA_close_dist_Prom_ENhA_1452$XEnhA,
                     y= NPC_dA_close_dist_Prom_ENhA_1452$Yprom, 
                     yend=NPC_dA_close_dist_Prom_ENhA_1452$YEnhA),
                 col="#c53bc7")+
  
  geom_segment(data=NPC_dA_close_dist_Prom_ENhB_1452,
               aes(x=NPC_dA_close_dist_Prom_ENhB_1452$Xprom, 
                   xend=NPC_dA_close_dist_Prom_ENhB_1452$XEnhb,
                   y= NPC_dA_close_dist_Prom_ENhB_1452$Yprom, 
                   yend=NPC_dA_close_dist_Prom_ENhB_1452$YEnhb),
               col="#23f7d3")+
  
  geom_segment(data=NPC_dA_close_dist_EnhA_Enhb_1452,
               aes(x=NPC_dA_close_dist_EnhA_Enhb_1452$XEnhA, 
                   xend=NPC_dA_close_dist_EnhA_Enhb_1452$XEnhb,
                   y= NPC_dA_close_dist_EnhA_Enhb_1452$YEnhA, 
                   yend=NPC_dA_close_dist_EnhA_Enhb_1452$YEnhb),
               col="#f6880c")+
  
  geom_polygon(data=NPC_dA_Tricontact_1452, aes(x=Dim1, y=Dim2, group=ID),  fill=NA, colour="red")+
  xlim(-500,500)+ylim(-500,500)+ggtitle("NPC dA")
bbis_and_dist


df <-
  data.frame(
    distance = rep(c("Prom EnhA", "EnhA EnhB", "EnhB Prom", "Tricontact"), 2),
    condition = c(rep("_WT", 4), rep("dA", 4)),
    count = c(
      NROW(NPC_WT_close_distance_Prom_ENhA),
      NROW(NPC_WT_close_distance_EnhA_Enhb),
      NROW(NPC_WT_close_distance_Prom_ENhB),
      NROW(NPC_WT_Tricontact)/3,
      NROW(NPC_dA_close_dist_Prom_ENhA_1452),
      NROW(NPC_dA_close_dist_EnhA_Enhb_1452),
      NROW(NPC_dA_close_dist_Prom_ENhB_1452),
      NROW(NPC_dA_Tricontact_1452)/3
    )
  )
df$distance <- factor(df$distance, levels = c("Prom EnhA", "EnhA EnhB", "EnhB Prom", "Tricontact"))
#df2 <- rbind(NPC_WT_close_distance_EnhA_Enhb[,c(5,6,7)],
#             NPC_WT_close_distance_Prom_ENhA[,c(5,6,7)],
#             NPC_WT_close_distance_Prom_ENhB[,c(5,6,7)],
#             NPC_WT_Tricontact[,c(5,6,7)],
#             NPC_dA_close_dist_EnhA_Enhb_1452[,c(5,6,7)],
#             NPC_dA_close_dist_Prom_ENhA_1452[,c(5,6,7)],
#             NPC_dA_close_dist_Prom_ENhB_1452[,c(5,6,7)],
#             NPC_dA_Tricontact_1452[,c(5,6,7)])



dfplot <- ggplot(df, aes(x=distance,y=count,fill=condition))+
  geom_bar(stat = 'identity', position=position_dodge())+
 # scale_fill_manual(values = c('purple','yellow','cyan','red'))
  theme_minimal()
  

#pdf(paste0("R_outpout/superimposition_T",My_Threshold,"_zoom.pdf"), height = 10, width = 20)
grid.arrange(a_and_dist, bbis_and_dist, dfplot, ncol=3)
#dev.off()

#Calculate area for WT and dA
#Perimeter
NPC_dA_484 <- NPC_dA[1:484,]
NPC_WT$perimeter <- (NPC_WT$F1F2+ NPC_WT$F1F3+ NPC_WT$F2F3)
NPC_dA_484$perimeter <- (NPC_dA_484$F1F2+ NPC_dA_484$F1F3+ NPC_dA_484$F2F3)
#Semi-Perimeter
NPC_WT$semiperimeter <- (NPC_WT$perimeter/2)
NPC_dA_484$semiperimeter <- (NPC_dA_484$perimeter/2)
#Area
NPC_WT$area <- sqrt(NPC_WT$semiperimeter*(NPC_WT$semiperimeter - NPC_WT$F1F2)*(NPC_WT$semiperimeter - NPC_WT$F2F3)*(NPC_WT$semiperimeter - NPC_WT$F1F3))
NPC_dA_484$area <- sqrt(NPC_dA_484$semiperimeter*(NPC_dA_484$semiperimeter - NPC_dA_484$F1F2)*(NPC_dA_484$semiperimeter - NPC_dA_484$F2F3)*(NPC_dA_484$semiperimeter - NPC_dA_484$F1F3))

NPC_dA_484$condition <- "dA"
NPC_WT$condition <- "_WT"

AllCondition<- rbind(NPC_WT, NPC_dA_484)
AllCondition$area <- AllCondition$area/1000

pdf("R_outpout/area_all_condition.pdf",width = 5, height = 5)
ggviolin(AllCondition, x="condition", y="area", fill="condition", palette = "Dark2",
         add=c("boxplot"), add.params=list(fill="white", width=0.1),
         width = 0.5, order=c("_WT", "dA"), title = "All points")+
  stat_compare_means( method = "wilcox.test" )
dev.off()

area_tricontact_NPC <- tAfterAPC_NPC_WT[tAfterAPC_NPC_WT$ID %in% unique(NPC_WT_Tricontact$ID),]
area_tricontact_dA <-  tAfterAPC_484[tAfterAPC_484$ID %in% unique(NPC_dA_Tricontact_1452$ID),]
#Perimeter
area_tricontact_NPC$perimeter <- (area_tricontact_NPC$d_Prom_EnhA+ area_tricontact_NPC$d_EnhA_EnhB+ area_tricontact_NPC$d_Prom_Enhb)
area_tricontact_dA$perimeter <- (area_tricontact_dA$d_Prom_EnhA+ area_tricontact_dA$d_EnhA_EnhB+ area_tricontact_dA$d_Prom_Enhb)
#semi perimeter
area_tricontact_NPC$semiperimeter <- (area_tricontact_NPC$perimeter/2)
area_tricontact_dA$semiperimeter <- (area_tricontact_dA$perimeter/2)
#area
area_tricontact_NPC$area <- sqrt(area_tricontact_NPC$semiperimeter*(area_tricontact_NPC$semiperimeter - area_tricontact_NPC$d_Prom_EnhA)*(area_tricontact_NPC$semiperimeter - area_tricontact_NPC$d_EnhA_EnhB)*(area_tricontact_NPC$semiperimeter - area_tricontact_NPC$d_Prom_Enhb))
area_tricontact_dA$area <- sqrt(area_tricontact_dA$semiperimeter*(area_tricontact_dA$semiperimeter - area_tricontact_dA$d_Prom_EnhA)*(area_tricontact_dA$semiperimeter - area_tricontact_dA$d_EnhA_EnhB)*(area_tricontact_dA$semiperimeter - area_tricontact_dA$d_Prom_Enhb))

area_tricontact_NPC$condition <- "_WT"
area_tricontact_dA$condition <- "dA"
All_area_tricontact <- rbind(area_tricontact_NPC, area_tricontact_dA)
All_area_tricontact$area <- All_area_tricontact$area/1000 
pdf(paste0("R_outpout/area_tricontact_T",My_Threshold,".pdf"), width = 5, height = 5)
ggviolin(All_area_tricontact, x="condition", y="area", fill="condition", palette = "Dark2",
         add=c("boxplot"), add.params=list(fill="white", width=0.1),
         width = 0.5, order=c("_WT", "dA"), title = "Tricontact t=250nm")+
  stat_compare_means( method = "wilcox.test", label.y = 75)+
  ylab("Area (µm2)")
dev.off()
