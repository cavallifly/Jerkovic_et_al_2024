#https://stackoverflow.com/questions/74066714/when-having-plotted-a-half-violin-plot-for-2-conditions-in-ggplot-how-do-i-chan

options(warn=-1)

library(see)
#library(scales) # to access break formatting functions
library(tidyverse)
library(ggplot2)
#library(dplyr)
library(ggpubr)
library(ggrepel)

args = commandArgs(trailingOnly=TRUE)

h = as.numeric(args[1])

data <- read.table('_tmp.tab')
colnames(data) <- c("ind1","ind2","rank","tag")

data$loop <- paste(data$ind1, data$ind2, sep="_")

levels = c("Models","CHi-C")
#colors=rainbow(length(levels))
colors=c("red","blue")
print(head(data))

summ <- data %>%
  group_by(across(all_of(c("tag")))) %>%
  summarize(n = n(), score = h)
print(summ)

pdf('_tmp.pdf')

compare_means(rank ~ tag,  data = data, method = "wilcox.test")
write.table(compare_means(rank ~ tag,  data = data, method = "wilcox.test"), file = "clipboard.tab", quote=F)

my_comparisons <- list(c("Mod", "Exp"))

#dotplot2 <- ggplot(data, aes(x = tag, y = rank, fill = tag)) +
#  geom_violinhalf(flip = c(1, 3)) +
#  geom_point() +
#  geom_line(aes(x = tag, y = rank, group = loop)) +
  #  geom_text_repel(aes(label = loop),box.padding = 0.5) + #, max.overlaps = Inf, ) +
#  scale_fill_manual(values = c("red","blue")) + 
#  labs(
#    title = '', 
#    subtitle = ''  
#  ) +
#  xlab("")+                                                                                              
#  theme(text=element_text(family = "Arial", size = 28),                                                                                                       
#        plot.title = element_text(face="bold", hjust=0),
#        plot.subtitle = element_text(face="italic", hjust=0),
#        axis.title.y = element_text(size=28, face="bold"),
#        axis.title.x = element_text(size=28, face="bold"),
#        legend.position="none", # Remove legend
#        aspect.ratio = 1 # Make the image square
#  ) +
#  ylab("Loop ranks") +  
#  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label="p") + #+ guides(color="none")
#  theme_classic(axis.title.y = element_text(size=28, face="bold"), axis.title.x = element_text(size=28, face="bold"))

#print(dotplot2)

dotplot2 <- ggplot(data, mapping=aes(x=tag, y=rank, fill=tag), alpha=0.7) + geom_violin(position=position_dodge(1), trim=F, scale="width") + geom_boxplot(position=position_dodge(1), width=0.1, color="black", fill="white", outlier.shape = NA) + scale_color_manual(values=colors) + scale_fill_manual(values=colors) + theme(panel.background = element_rect(fill = NA), panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(), axis.title = element_text(size=28), axis.text.x = element_text(angle=60, hjust=1), axis.text = element_text(size=24), axis.ticks.x = element_blank(), legend.position = "none", legend.title =  element_blank()) + labs(x="",y="Loop ranks") #+ stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label="p") + #+ guides(color="none")
print(dotplot2)

#geom_text(aes(label=n), color="black", data = summ, position=position_dodge(1))
#geom_text(aes(label=n), color="black", data = summ, position=position_dodge(1)) 

dev.off()
#conflicts()
