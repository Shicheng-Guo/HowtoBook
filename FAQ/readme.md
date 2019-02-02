How to Calculate Genomic Inflation Factor and λgc for GWAS

You have conducted your genome-wide association study (GWAS) and have tested each genetic variant for an association with your trait of interest. Now it is time to investigate if there are any systematic biases that may be present in your association results. A common way to do this is to calculate the genomic inflation factor, also known as lambda gc (λgc). By definition, λgc is defined as the median of the resulting chi-squared test statistics divided by the expected median of the chi-squared distribution. The median of a chi-squared distribution with one degree of freedom is 0.4549364. A λgc value can be calculated from z-scores, chi-square statistics, or p-values, depending on the output you have from the association analysis. Follow these simple steps to calculate lambda GC using R programming language.

(1) Convert your output to chi-squared values
```
	# For z-scores, just square them
	chisq <- data$z^2

# For chi-squared values, keep as is
	chisq <- data$chisq

	# For p-values, calculate chi-squared statistic
	chisq <- qchisq(1-data$pval,1)
```
(2) Calculate lambda gc (λgc)
```
  median(chisq)/qchisq(0.5,1)
```
If analysis results your data follows the normal chi-squared distribution, the expected λgc value is 1. If the λgc value is greater than 1, then this may be evidence for some systematic bias that needs to be corrected in your analysis.
