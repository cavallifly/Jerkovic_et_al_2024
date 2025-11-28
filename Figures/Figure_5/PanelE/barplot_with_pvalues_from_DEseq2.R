# Script inspired by https://www.datanovia.com/en/lessons/anova-in-r/#two-way-independent-anova

library(tidyverse)
library(ggpubr)
library(rstatix)
library(dunn.test)
library(dplyr)
library(FSA)


args <- commandArgs(trailingOnly = TRUE)

target <- args[[1]]

if(target == "H3K27ac")
{
    yMax = 2500
}
if(target == "Pol2")
{
    yMax = 5000
}
if(target == "Pax6")
{
    yMax = 1600
}
if(target == "ATACseq")
{
    yMax = 4000
}


inFiles = list.files("./",pattern="normalized.*ANOVA.*bed.*")
print(inFiles)
#quit()

#subsetting <- "FALSE"
subsetting <- "TRUE"

pvaluesFromDEseq2 = read.table("results_from_DEseq2_analysis.tsv",header=T)
pvaluesFromDEseq2$pvalue <- format(pvaluesFromDEseq2$pvalue, scientific=TRUE, digits=3)

for(inFile in inFiles)
{
    print(inFile, quote=F)
    allData <- read.table(inFile, header=F)
    ###
    colnames(allData) <- c("locus","values")    
    print(head(allData), quote=F)
    #quit()
    name <- gsub(".bed","",inFile)    

    for(condition1 in c("wt"))
    {
	#print(unique(gsub("C4_enhA_","",allData[grep("enhA",allData$locus),]$locus)))
        for(condition2 in unique(gsub("C4_enhA_","",allData[grep("enhA",allData$locus),]$locus)))
    	{
	    condition2 <- tail(str_split(condition2, pattern="_")[[1]], n=1)
	    print(condition2)
	    if(condition1 == condition2){next;}

	    data <- allData[grep(condition1,allData$locus),]	
	    for(r in 1:nrow(allData))
	    {
		cond <- tail(str_split(allData[r,], pattern="_")[[1]], n=1)
		if(cond == condition2)
		{
		    data  <- rbind(data,allData[r,])
		}		
	    }
	    colnames(data) <- c("locus","values")


    # List of specific comparisons (custom pairs)
    custom_pairs <- list(c("Ci_1_wt", paste0("Ci_1_",condition2)), c("Prom_wt", paste0("Prom_",condition2)), c("C1_dprom_wt", paste0("C1_dprom_",condition2)), c("C2_wt", paste0("C2_",condition2)), c("C3_wt", paste0("C3_",condition2)), c("C4_enhA_wt", paste0("C4_enhA_",condition2)), c("C5_Badj_wt", paste0("C5_Badj_",condition2)), c("C6_enhB_wt", paste0("C6_enhB_",condition2)), c("controlPeak1_wt", paste0("controlPeak1_",condition2)), c("controlPeak2_wt", paste0("controlPeak2_",condition2)), c("controlPeak3_wt", paste0("controlPeak3_",condition2)))
    print(custom_pairs)
    #quit()

    levels <- c("Ci_1_wt", paste0("Ci_1_",condition2), "Prom_wt", paste0("Prom_",condition2), "C1_dprom_wt", paste0("C1_dprom_",condition2), "C2_wt", paste0("C2_",condition2), "C3_wt", paste0("C3_",condition2), "C4_enhA_wt", paste0("C4_enhA_",condition2), "C5_Badj_wt", paste0("C5_Badj_",condition2), "C6_enhB_wt", paste0("C6_enhB_",condition2), "controlPeak1_wt", paste0("controlPeak1_",condition2), "controlPeak2_wt", paste0("controlPeak2_",condition2), "controlPeak3_wt", paste0("controlPeak3_",condition2))
    print(levels)
    data$locus <- factor(data$locus, levels = levels)
    print(data)
    #quit()

    print(paste0("### Get summary statistics ###"), quote=F)
    dataStats <- data %>%
	  group_by(locus) %>%
	  get_summary_stats(values, type = "mean_sd")
    print(as.data.frame(dataStats), quote=F)
    outFileTextMean <- paste0("meanAndStdDev_",name,"_",condition1,"_",condition2,"_with_stats.tsv")
    write.table(as.data.frame(dataStats), file=outFileTextMean,  row.names=F, sep="\t")  
    #quit()

    summ <- data %>%
           group_by(locus) %>%
	   summarize(n = n(), locus=locus, score = ((1+0.15)*min(data$values,-100)))
    summ <- unique(summ)
    print(as.data.frame(summ), quote=F)

    yMax <- yMax + 100
    yMin <- (1+0.05)*unique(summ$score)

    bxp <- ggplot() +
	   geom_point(data=data, aes(x=locus,y=values), position = position_jitter(width=0.2)) +
    	   #ggplot(dataStats, aes(x=locus, y=mean, fill=locus)) +    	   
           geom_bar(data=dataStats, aes(x=locus, y=mean, fill=locus), position=position_dodge(), stat="identity", colour='black', alpha=0.3) +
           geom_errorbar(data=dataStats, position=position_dodge(), aes(x=locus, ymin=mean-sd, ymax=mean+sd), width=0.1) +
	   labs(y=paste0(target," signal"), x = "") +
	   geom_text(data=summ, aes(x=locus, y=score, label = paste0("n = ",n)), size=2., color="black") +    
	   theme_classic() +
   	   theme(axis.text.x = element_text(angle = 60, hjust=1), legend.position = "none") +
	   ylim(yMin,yMax)

	
    outFile <- paste0("barplot_oneWay_Zfp608_",name,"_",condition1,"_",condition2,".pdf")
    pdf(outFile)
    print(bxp)
    dev.off()

    if(summ[summ$locus == paste0("C3_",condition2),]$n == 1)
    {	
        next
    }

    print(paste0("### Check assumptions ###"), quote=F)
    print(paste0("1) Checking the presence of outliers"), quote=F)
    checkOutliers <- data %>%
	  group_by(locus) %>%
	  identify_outliers(values)
    if(nrow(checkOutliers) == 0)
    {
        print(paste0("No outliers found: This assumption is verified"), quote=F)
    } else {
        print(paste0("### WARNING: The following data-points are outiers!"), quote=F)
        print(checkOutliers, quote=F)
        print(paste0("Note that, in the situation where you have extreme outliers, this can be due to:"), quote=F)
        print(paste0("1) data entry errors, measurement errors or unusual values."), quote=F)
        print(paste0("You can include the outlier in the analysis anyway if you do not believe the result will be substantially affected."), quote=F)
        print(paste0("This can be evaluated by comparing the result of the ANOVA test with and without the outlier."), quote=F)
        print(paste0("It’s also possible to keep the outliers in the data and perform robust ANOVA test using the WRS2 package."), quote=F)
        print("", quote=F)     	    
    }
    print("", quote=F)


    print(paste0("2) Checking normality assumption"), quote=F)
    print(paste0("# Building the linear model"), quote=F)
    model  <- lm(values ~ locus,
    	         data = data)
    print(paste0("# Creating a QQ plot of residuals"), quote=F)
    p <- ggqqplot(residuals(model))
    pdf(paste0("qqplot_oneWay_Zfp608_",name,"_",condition1,"_",condition2,".pdf")) 
    print(p)
    dev.off()	
    print(paste0("### Compute Shapiro-Wilk test of normality ###"), quote=F)
    checkNormality <- shapiro_test(residuals(model))
    print(checkNormality, quote=F)
    if(checkNormality$p.value > 0.05)
    {
        print(paste0("In the QQ plot, as all the points fall approximately along the reference line, we can assume normality."), quote=F)
	print(paste0("This conclusion is supported by the Shapiro-Wilk test. The p-value is not significant (p = ",checkNormality$p.value,"),"), quote=F)
	print(paste0("so we can assume normality."), quote=F)
    } else {
        print(paste0("The normality assumption is not supported by the Shapiro-Wilk test. The p-value is significant (p = ",checkNormality$p.value,"),"), quote=F)
        print(paste0("so we cannot assume normality for this dataset."), quote=F)
        #    quit()
    }

    print("", quote=F)
    print(paste0("3) Checking homogeneity of variance assumption"), quote=F)
    print(paste0("This can be checked using the Levene’s test:"), quote=F)
    checkHomogenityOfVariance <- data %>% levene_test(values ~ locus)
    print(checkHomogenityOfVariance, quote=F)
    if(checkHomogenityOfVariance$p > -0.05)
    {
        print(paste0("The Levene’s test is not significant (p > ",checkHomogenityOfVariance$p,"). Therefore, we can assume the homogeneity of variances in the different groups."), quote=F)
	print(paste0("### Computation of the one-way ANOVA test ###"), quote=F)
	print(paste0("In the R code below, the asterisk represents the interaction effect and the main effect of each variable (and all lower-order interactions)."), quote=F)
	res.aov <- data %>% anova_test(values ~ locus)
	print(res.aov, quote=F)

	significantInteraction = res.aov[res.aov$p < 0.05,]
	if(nrow(significantInteraction) > 0)
	{
	    print(paste0("There was a statistically significant interaction between locis"), quote=F)
	    print(significantInteraction, quote=F)
	} else {
	    print(paste0("There was no statistically significant interaction between loci"), quote=F)
	    #next
	}

	print("", quote=F)	
	print(paste0("### Post-hoct tests ###"), quote=F)

	print(paste0("### Procedure for significant one-way interaction ###"), quote=F)
	print(paste0("### Compute pairwise comparisons ###"), quote=F)

	print(paste0("To determine which group means are different. We’ll now perform multiple pairwise comparisons between the different conditions."), quote=F)

	print(paste0("You can run and interpret all possible pairwise comparisons using a Bonferroni adjustment."), quote=F)
	print(paste0("This can be easily done using the function emmeans_test() [rstatix package], a wrapper around the emmeans package,"), quote=F)
	print(paste0("which needs to be installed. Emmeans stands for estimated marginal means (aka least square means or adjusted means)."), quote=F)
	

	print(paste0("Compare the values of the different cellType levels by condition levels:"), quote=F)
	pairwiseTestsCondition <- data %>%
			       pairwise_t_test(
			       values ~ locus, 
			       p.adjust.method = "bonferroni"
			    )		      

	if(subsetting == "TRUE")
	{
	# Filter only desired comparisons
	if(exists("subset_t_test"))
        {
	    rm(subset_t_test)
	}
	i = 0
	ypos <- c()
	for(pair in custom_pairs)
	{
	    delta <- 50
	    yposStart <- yMax - delta
	    #yposStart <- max(data$values)+delta
	    
	    print(paste0(pair[[1]]," ",pair[[2]]))
	    subset <- pairwiseTestsCondition[(pairwiseTestsCondition$group1 == pair[[1]] & pairwiseTestsCondition$group2 == pair[[2]]) | (pairwiseTestsCondition$group1 == pair[[2]] & pairwiseTestsCondition$group2 == pair[[1]]),]
	    print(paste0("P-value from multiple testing ",subset$p.adj))
	    padjust = pvaluesFromDEseq2[pvaluesFromDEseq2$condition2 == pair[[2]],]$pvalue
	    subset$p.adj <- padjust
	    print(paste0("P-value from DEseq2 ",padjust))
	    print(subset)
	    
	    ypos <- c(ypos,yposStart + 0 * delta)
	    if(exists("subset_t_test"))
	    {
	        subset_t_test <- rbind(subset_t_test,subset)
	    } else {
		subset_t_test <- subset	       
	    }	    
	    i = i + 1
	}
	
	# Adjust p-values only for this subset
	#subset_t_test$p.adj <- p.adjust(subset_t_test$p, method = "bonferroni")
	pairwiseTestsCondition <- subset_t_test
	}
	print(as.data.frame(pairwiseTestsCondition), quote=F)
	#quit()
    } else {
       print(paste0("The Levene’s test is significant (p = ",checkHomogenityOfVariance$p," < 0.05). Therefore, we cannot assume the homogeneity of variances in the different groups."), quote=F)

       print(paste0("The Welch one-way test is an alternative to the standard one-way ANOVA"), quote=F)
       print(paste0("in the situation where the homogeneity of variance can’t be assumed (i.e., Levene test is significant)."), quote=F)

       res.aov <- data %>% welch_anova_test(values ~ locus)
       print(res.aov, quote=F)

       print(paste0("In this case, the Games-Howell post hoc test or pairwise t-tests (with no assumption of equal variances) can be used to compare all possible combinations of group differences."), quote=F)

       pairwiseTestsCondition <- data %>% games_howell_test(values ~ locus)
	    
    }
    print("", quote=F)	

    print(paste0("### Visualization: barplots with p-values"), quote=F)
    #pairwiseTestsCondition <- as.data.frame(pairwiseTestsCondition)    		    
    print(as.data.frame(pairwiseTestsCondition), quote=F)
    outFileText <- paste0("barplot_oneWay_Zfp608_",name,"_",condition1,"_",condition2,"_with_stats.tsv")
    print(outFileText)
    write.table(as.data.frame(pairwiseTestsCondition), file=outFileText,  row.names=F, sep="\t")    
    pairwiseTestsCondition <- pairwiseTestsCondition %>% 
    			      			       add_xy_position(x = "locus")
    pairwiseTestsCondition$y.position <- ypos 						      

    print(as.data.frame(pairwiseTestsCondition), quote=F)

    finalBxp <- bxp +
	   stat_pvalue_manual(data=pairwiseTestsCondition, label="p.adj", label.size = 2, bracket.size = 0.1, tip.length = 0, ) +
	   labs(subtitle = get_test_label(res.aov, detailed = TRUE), caption = get_pwc_label(pairwiseTestsCondition)) +
	   theme(legend.position = "none")

    outFile <- paste0("barplot_oneWay_Zfp608_",name,"_",condition1,"_",condition2,"_with_stats.pdf")
    pdf(outFile)
    print(finalBxp)
    dev.off()	
    #quit()
    print("", quote=F)
    }
    }
} # Close cycle over inFile