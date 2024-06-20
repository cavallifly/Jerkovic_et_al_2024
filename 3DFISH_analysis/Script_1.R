
library(paleomorph)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggpubr)

# Import and select the coordinates (x,y,z) of genomic loci
# Prepare input data frame for the Superimpostion
# df is the data frame with the coordinates of the loci : F1(x,y,z), F2(x,y,z), F3(x,y,z)
############################################################################
#####                    DATA FRAME STRUCTURE : df                     #####
############################################################################
    #ID     # xFC1 # yFC1 # zFC1 # xFC2 # yFC2 # zFC2 # xFC3 # yFC3 # zFC3 #
#nucleus 1  # ...  # ...  # ...  # ...  # ...  # ...  # ...  # ...  # ...  #
#nucleus 2  # ...  # ...  # ...  # ...  # ...  # ...  # ...  # ...  # ...  #
#nucleus 3  # ...  # ...  # ...  # ...  # ...  # ...  # ...  # ...  # ...  #
############################################################################

rm(list=ls())

#You must run this script for each conditions you want to study.
#It will create a dataframe with the aligned coordinates 
#then we will use the output of this script for down stream analysis

name_Of_The_Condition <- "NPC_WT" #Here it is an example from the paper where we study 2 conditions NPC_WT and NPC_dA

df[,-1]<- apply(df[,-1],2, as.numeric)

#Calculate distances then remove distances > CutOff
CutOff <-  1500 #Add your cut off in nanometer

#Add the size in nm of the voxel in x,y,z
#In our study, we use confocal microscopie
x=40
y=40
z=140
#dist FC1FC2
df$F1F2 <- sqrt(((df$xFC2-df$xFC1)*x)**2+((df$yFC2-df$yFC1)*y)**2+((df$zFC2-df$zFC1)*z)**2)
#dist FC2FC3
df$F2F3 <- sqrt(((df$xFC3-df$xFC2)*x)**2+((df$yFC3-df$yFC2)*y)**2+((df$zFC3-df$zFC2)*z)**2)
#dist FC1FC3
df$F1F3 <- sqrt(((df$xFC3-df$xFC1)*x)**2+((df$yFC3-df$yFC1)*y)**2+((df$zFC3-df$zFC1)*z)**2)
cut_off_mydist<- function(a){
  a <- a[a$F1F2 %in% subset(a$F1F2, a$F1F2 < CutOff),]
  a <- a[a$F2F3 %in% subset(a$F2F3, a$F2F3 < CutOff),]
  a <- a[a$F1F3 %in% subset(a$F1F3, a$F1F3 < CutOff),]
}
df <- cut_off_mydist(df)
#Select the data corresponding to the coordinates of the 3 different points
B <- df
Bx <- B[, c("xFC1", "xFC2", "xFC3")]
By <- B[, c("yFC1", "yFC2", "yFC3")]
Bz <- B[, c("zFC1", "zFC2", "zFC3")]
B <- cbind(Bx, By, Bz)
B <- B[!B$xFC1=="NaN",]
B[,1:3] <- B[,1:3]*x
B[,4:6] <- B[,4:6]*y
B[,7:9] <- B[,7:9]*z
str(B)
#Generate array with the coordinates for Procrustes Superimposition
num.coords <- dim(B)[2]
coor3D <- 3
ns <- dim(B)[1]
ma <- array(dim=c(num.coords/coor3D, coor3D, ns))
xi <- c(1,2,3); yi <- xi+3; zi <- yi+3
for (i in 1:ns) {
  ma[, 1, i] <- unlist(B[i, xi])
  ma[, 2, i] <- unlist(B[i, yi])
  ma[, 3, i] <- unlist(B[i, zi])
}

# Align the data 
aligned <- procrustes(ma, scale=F)

plotSpecimens(aligned,size= 5)
rm(B, Bx, By, Bz, xi, yi, zi)


##Create a new data frame with the coordinate of the 3 points after 
##procrustres superimposition
molecular_ID = c("Prom", "EnhA", "EnhB")   #you can change the name of the loci
PreAPC <- data.frame(
  ID = rep(NA, dim(aligned)[3]*3),
  X = rep(NA, dim(aligned)[3]*3),
  Y = rep(NA, dim(aligned)[3]*3),
  Z = rep(NA, dim(aligned)[3]*3),
  molecular_ID = rep(molecular_ID,dim(aligned)[3])
)
u <- 1; v <- 2; w <- 3
for (i in 1:dim(aligned)[3]) {
  PreAPC[u:w,1] <- i
  PreAPC[u,2:4] <- aligned[1,1:3,i]
  PreAPC[v,2:4] <- aligned[2,1:3,i]
  PreAPC[w,2:4] <- aligned[3,1:3,i]
  PreAPC[u:w,5] <- molecular_ID
  u <- u+3; v <- v+3; w <- w+3
  
}
rm(u,v,w)
res.pca <- PCA(PreAPC[,2:4],scale.unit = F, ncp=2, graph = FALSE)

           legend.title = "Groups"
#)
#To get the coordinates
ind <- get_pca_ind(res.pca)

AfterAPC <- data.frame(
  ID = PreAPC$ID,
  Dim1 = ind$coord[, 1],
  Dim2 = ind$coord[, 2],
  molecular_ID = PreAPC$molecular_ID
)
rm(PreAPC,res.pca,ind)


write.csv2(AfterAPC, paste0(name_Of_The_Condition,".csv"),row.names = F)