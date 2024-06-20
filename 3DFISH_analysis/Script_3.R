rm(list=ls())
library(paleomorph)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggpubr)
library(tidyverse)
library(splancs)
library(MASS)
library(fields)
setwd("C://Users/pc/Desktop/Temp from Archive/")


#######Prepare the df#########
#Import AfterAPC for WT and dA
NPC_WT_AfterAPC <- read.csv2("R_outpout/NPC_WT/NPC_WT_AfterAPC.csv", stringsAsFactors = F)
NPC_dA_AfterAPC <- read.csv2("R_outpout/NPC_dA//NPC_dA_AfterAPC.csv", stringsAsFactors = F)
NPC_dA_AfterAPC[,2] <- NPC_dA_AfterAPC[,2]*-1
NPC_dA_AfterAPC[,3] <- NPC_dA_AfterAPC[,3]*-1

#Generate a transformed data frame from the after APC df
#This new df should have 3 times less row than the original df after APC

#START OPTION
#Define the condition : WT or Mutant
condition <- NPC_WT_AfterAPC
#END OPTION
#FUtur amélioration : essayer de faire entrer ka génération cette df dans une loop
#changer le nom des positions génomiques pour que ce soit choisi par utilisateur (ici c'est EnhA et EnhB)
tAfterAPC <- data.frame(ID=rep(NA,NROW(paste0(condition))/3),
                        Xprom=rep(NA,NROW(condition)/3),
                        Yprom=rep(NA,NROW(condition)/3),
                        XEnhA=rep(NA,NROW(condition)/3),
                        YEnhA=rep(NA,NROW(condition)/3),
                        XEnhb=rep(NA,NROW(condition)/3),
                        YEnhb=rep(NA,NROW(condition)/3))
u <- 1
for (i in seq(1,NROW(condition),by=3)) {
  tAfterAPC[u,1] <- condition[i,1]
  tAfterAPC[u,2:3] <- condition[i,2:3]
  tAfterAPC[u,4:5] <- condition[i+1,2:3]
  tAfterAPC[u,6:7] <- condition[i+2,2:3]
  u <- u+1
}
#Calculation of the distances
#dist Prom EnhA
tAfterAPC$d_Prom_EnhA <- sqrt((tAfterAPC$XEnhA-tAfterAPC$Xprom)**2+(tAfterAPC$YEnhA-tAfterAPC$Yprom)**2)
#dist EnhA EnhB
tAfterAPC$d_EnhA_EnhB <- sqrt((tAfterAPC$XEnhb-tAfterAPC$XEnhA)**2+(tAfterAPC$YEnhb-tAfterAPC$YEnhA)**2)
#dist Prom EnhB
tAfterAPC$d_Prom_Enhb <- sqrt((tAfterAPC$XEnhb-tAfterAPC$Xprom)**2+(tAfterAPC$YEnhb-tAfterAPC$Yprom)**2)

#Save the new df for latter and delete the temp
tAfterAPC_NPC_WT <- tAfterAPC
rm(condition, tAfterAPC)


#START OPTION
#Define the condition : WT or dA
condition <- NPC_dA_AfterAPC
#END OPTION
tAfterAPC <- data.frame(ID=rep(NA,NROW(paste0(condition))/3),
                        Xprom=rep(NA,NROW(condition)/3),
                        Yprom=rep(NA,NROW(condition)/3),
                        XEnhA=rep(NA,NROW(condition)/3),
                        YEnhA=rep(NA,NROW(condition)/3),
                        XEnhb=rep(NA,NROW(condition)/3),
                        YEnhb=rep(NA,NROW(condition)/3))
u <- 1
for (i in seq(1,NROW(condition),by=3)) {
  tAfterAPC[u,1] <- condition[i,1]
  tAfterAPC[u,2:3] <- condition[i,2:3]
  tAfterAPC[u,4:5] <- condition[i+1,2:3]
  tAfterAPC[u,6:7] <- condition[i+2,2:3]
  u <- u+1
  
}
#Calcule of the distances
#dist Prom EnhA
tAfterAPC$d_Prom_EnhA <- sqrt((tAfterAPC$XEnhA-tAfterAPC$Xprom)**2+(tAfterAPC$YEnhA-tAfterAPC$Yprom)**2)
#dist EnhA EnhB
tAfterAPC$d_EnhA_EnhB <- sqrt((tAfterAPC$XEnhb-tAfterAPC$XEnhA)**2+(tAfterAPC$YEnhb-tAfterAPC$YEnhA)**2)
#dist Prom EnhB
tAfterAPC$d_Prom_Enhb <- sqrt((tAfterAPC$XEnhb-tAfterAPC$Xprom)**2+(tAfterAPC$YEnhb-tAfterAPC$Yprom)**2)

#Save the new df for latter and delete the temp
tAfterAPC_NPC_dA <- tAfterAPC
rm(condition, tAfterAPC)
NPC_dA_AfterAPC <- NPC_dA_AfterAPC[1:1452,]
tAfterAPC_NPC_dA <- tAfterAPC_NPC_dA[1:484,]
tAfterAPC_NPC_dA$condition <- "dA"
tAfterAPC_NPC_WT$condition <-  "WT"


#################################
#######Start the analysis########
#################################



#####################################
####    With 2d Kernel Density   ####
#####################################

Area_kd <- data.frame(step=rep(NA,60),
                      Area.WT=rep(NA,60),
                      Area.dA=rep(NA,60),
                      loci=rep(NA,60))
df_wt <- data.frame(
  x = rep(NA, 60),
  y = rep(NA, 60),
  z = rep(NA, 60),
  loci = rep(NA, 60))
df_dA <- data.frame(
  x = rep(NA, 60),
  y = rep(NA, 60),
  z = rep(NA, 60),
  loci = rep(NA, 60))
j=1
for (loci in unique(NPC_WT_AfterAPC$molecular_ID)) {
    #for WT
    x= NPC_WT_AfterAPC[NPC_WT_AfterAPC$molecular_ID==loci,2]
    y= NPC_WT_AfterAPC[NPC_WT_AfterAPC$molecular_ID==loci,3]
    f1 <- kde2d(x=x, y=y,n=20, lims = c(range(x), range(y)))
    z=interp.surface(f1,loc = matrix(c(x,y), ncol=2))
    df <- data.frame(x=x, y=y, z=z)
    df$loci <- loci
    df_wt <- rbind(df_wt, df)
    u=j
    for (i in seq(from=0.05, to=1, by=0.05)) {
      order <- df[order(df$z, decreasing=T),]
      order <- order[1:(484*i),] 
      hull <- order %>% slice(chull(x, y))
      Area_kd[u,1] <- i
      Area_kd[u,2] <- as.matrix(hull[,1:2]) %>% areapl
      Area_kd[u,4] <- loci
      u=u+1
    }
    
    #for dA
    x= NPC_dA_AfterAPC[NPC_dA_AfterAPC$molecular_ID==loci,2]
    y= NPC_dA_AfterAPC[NPC_dA_AfterAPC$molecular_ID==loci,3]
    
    f1 <- kde2d(x=x, y=y,n=20, lims = c(range(x), range(y)))
    z=interp.surface(f1,loc = matrix(c(x,y), ncol=2))
    df <- data.frame(x=x, y=y,z=z)
    df$loci <- loci
    df_dA <- rbind(df_dA, df)
    u=j
    for (i in seq(from=0.05, to=1, by=0.05)) {
      order <- df[order(df$z, decreasing=T),]
      order <- order[1:(484*i),] 
      hull <- order %>% slice(chull(x, y))
      Area_kd[u,1] <- i
      Area_kd[u,3] <- as.matrix(hull[,1:2]) %>% areapl
      u=u+1
    }
   j=u
  }

tArea_kd <- melt(Area_kd, id.vars = c("step", "loci"))

ggpaired(
  tArea_kd[tArea_kd$loci=="EnhA",], x = "variable", y = "value", color = "variable", palette = "jco", 
  line.color = "gray", line.size = 0.4)+
  yscale("log10")+#, ylim = c(0, 40))
  stat_compare_means(paired = T )
wilcox.test(x=tArea_kd[tArea_kd$loci=="EnhA" & tArea_kd$variable=="Area.WT",4], 
            y=tArea_kd[tArea_kd$loci=="EnhA" & tArea_kd$variable=="Area.dA",4], 
             paired = T, p.adjust.method = "fdr",)


#make the barplots to compare Area in WT and dA for the 3 loci
for (loci in unique(NPC_WT_AfterAPC$molecular_ID)) {
  a <- ggpaired(
    tArea_kd[tArea_kd$loci==loci,], x = "variable", y = "value", color = "variable", palette = "jco", 
    line.color = "gray", line.size = 0.4)+
    yscale("log10")+#, ylim = c(0, 40))
    stat_compare_means(paired = T)+
    ylab("Area nm2")+
    xlab("")+
    ggtitle(loci)+
    theme(legend.position='none')
  print(a)
  pdf(paste0(loci,"_boxplot_area.pdf"), height = 4, width = 3 )
  print(a)
  dev.off()
}

a1 <-
  ggplot(df_wt, aes(x=x, y=y, colour=z))+
  geom_point(size=0.2 )+
  scale_color_gradientn(colours = rainbow(20), limits=c(min(df_dA$z, na.rm = T), max(df_wt$z, na.rm = T )))+
  ylim(-1250,1000)+xlim(-1500, 1250)+
  theme(legend.position='none')+
  ggtitle("wt")

a2 <- 
  ggplot(df_dA, aes(x=x, y=y, colour=z))+
  geom_point(size=0.2)+
  scale_color_gradientn(colours = rainbow(20), limits=c(min(df_dA$z, na.rm = T), max(df_wt$z, na.rm = T )))+
  ylim(-1250,1000)+xlim(-1500, 1250)+
  ggtitle("dA")
#pdf("kernel_density_WT_dA.pdf", height = 5, width = 8)
grid.arrange(a1, a2, ncol=2, widths= c(1.5,2))
#dev.off()

