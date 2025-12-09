# Script inspired by https://www.datanovia.com/en/lessons/anova-in-r/#two-way-independent-anova

library(tidyverse)
library(ggpubr)
library(rstatix)

inFiles = list("ZFP608_ddCt_values_dprom_polyA_for_anova.bed")

for(inFile in inFiles)
{
    print(inFile, quote=F)
    allData <- read.table(inFile, header=F)
    ###
    colnames(allData) <- c("primer","cellType", "condition", "ddCt")
    allData <- allData[grep("NPC",allData$cellType),]
    allData$values <- 2^(-allData$ddCt)
    ###
    print(head(allData), quote=F)

    primers = unique(allData$primer)
    for(primer in primers)
    {
	print(primer, quote=F)
        data = allData[grep(primer,allData$primer),]
        data$condition   <- paste0(data$cellType,"_",data$condition)

	levels <- c("NPC_wt","NPC_dprom","NPC_polyA")	
        data$condition <- factor(data$condition, levels = levels)

	print(paste0("### Get summary statistics ###"), quote=F)
	dataStats <- data %>%
	  group_by(condition) %>%
	  get_summary_stats(values, type = "mean_sd")
	print(dataStats, quote=F)

	bxp <- ggboxplot(
	  data, x = "condition", y = "values",
	  color = "cellType", palette = "jco"
	) +
        geom_text(data=dataStats, aes(y=0.02, label = n), color="black") +
        theme(axis.text.x = element_text(angle = 60, hjust=1))
	
	outFile <- paste0("boxplot_oneWay_Zfp608_",primer,".pdf")
	pdf(outFile)
	print(bxp)
	dev.off()

	print(paste0("### Check assumptions ###"), quote=F)
	print(paste0("1) Checking the presence of outliers"), quote=F)
	checkOutliers <- data %>%
	  group_by(condition) %>%
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
	#next

	print(paste0("2) Checking normality assumption"), quote=F)
	print(paste0("# Building the linear model"), quote=F)
	model  <- lm(values ~ condition,
             data = data)
	print(paste0("# Creating a QQ plot of residuals"), quote=F)
	p <- ggqqplot(residuals(model))
	pdf(paste0("qqplot_oneWay_Zfp608_",primer,".pdf")) 
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
	#next

	print(paste0("3) Checking homogeneity of variance assumption"), quote=F)
	print(paste0("This can be checked using the Levene’s test:"), quote=F)
	checkHomogenityOfVariance <- data %>% levene_test(values ~ condition)
	print(checkHomogenityOfVariance, quote=F)
	#if(checkHomogenityOfVariance$p > 0.00)
	#{
	    print(paste0("The Levene’s test is not significant (p > ",checkHomogenityOfVariance$p,"). Therefore, we can assume the homogeneity of variances in the different groups."), quote=F)

	    print(paste0("### Computation of the one-way ANOVA test ###"), quote=F)
	    print(paste0("In the R code below, the asterisk represents the interaction effect and the main effect of each variable (and all lower-order interactions)."), quote=F)
	    res.aov <- data %>% anova_test(values ~ condition)
	    print(res.aov, quote=F)

	    significantInteraction = res.aov[res.aov$p < 0.052,]
	    if(nrow(significantInteraction) > 0)
	    {
	        print(paste0("There was a statistically significant interaction between cell type and condition for 2^(-ddCt)"), quote=F)
	        print(significantInteraction, quote=F)
	    } else {
	        print(paste0("There was no statistically significant interaction between cell type and condition for 2^(-ddCt)"), quote=F)
		#	        next
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
	         values ~ condition,
		 alternative = "two.sided",
		 p.adjust.method = "bonferroni",
		 pool.sd = T
		 )	
	    #pairwiseTestsCondition <- data %>%
	    #			   tukey_hsd(values ~ condition)
	    #print(paste0(pairwiseTestsCondition), quote=F)	
	   
	#}
	#else {
	#   print(paste0("The Levene’s test is significant (p = ",checkHomogenityOfVariance$p," < 0.05). Therefore, we cannot assume the homogeneity of variances in the different groups."), quote=F)

	#   print(paste0("The Welch one-way test is an alternative to the standard one-way ANOVA"), quote=F)
	#   print(paste0("in the situation where the homogeneity of variance can’t be assumed (i.e., Levene test is significant)."), quote=F)

	#   res.aov <- data %>% welch_anova_test(values ~ condition)
	#   print(res.aov, quote=F)

	#   print(paste0("In this case, the Games-Howell post hoc test or pairwise t-tests (with no assumption of equal variances) can be used to compare all possible combinations of group differences."), quote=F)

	#   pairwiseTestsCondition <- data %>% games_howell_test(values ~ condition)
	#   #print(paste0(pairwiseTestsCondition), quote=F)
	    
	#}
	print("", quote=F)	

	print(paste0("### Visualization: box plots with p-values"), quote=F)
	print(data.frame(pairwiseTestsCondition), quote=F)	
	pairwiseTestsCondition <- pairwiseTestsCondition %>% add_xy_position(x = "condition")
	#	print(pairwiseTestsCondition, quote=F)

	finalBxp <- bxp +
	   stat_pvalue_manual(pairwiseTestsCondition[pairwiseTestsCondition$p.adj < 0.05,]) +
	   labs(
	       subtitle = get_test_label(res.aov, detailed = TRUE),
	       caption = get_pwc_label(pairwiseTestsCondition)
	   )
	outFile <- paste0("boxplot_oneWay_Zfp608_",primer,"_with_stats.pdf")
	pdf(outFile)
	print(finalBxp)
	dev.off()	

	print("", quote=F)	
    } # Close cycle over primer
} # Close cycle over inFile
