# gwas-power
**R Functions** to calculate power of GWAS studies for a single associated SNP, under various parameters. Suitable for classical (i.e. single-SNP single-trait) GWAS studies using linear regression models, i.e for quantitative traits.

Uses formulae for power calculation presented in *Appendix A* of *Visscher PM, Wray NR, Zhang Q, et al. 10 Years of GWAS Discovery: Biology, Function, and Translation. Am J Hum Genet 2017;101(1):5-22. doi: 10.1016/j.ajhg.2017.06.005.*
* Assumes that additional covariates, if any, are uncorrelated to the SNP.
* It is common to include genetic PCs as covariates in GWAS to adjust for stratification or admixture; in such cases, if the SNP is severely stratified, these formulae will not be applicable, and a simulation-based approach might be preferrable.
* Uses a chi-square statistic, equivalent to using the z-statistic in a regression (the chi-square statistic is the square of the z-statistic). The results should be very similar to using tests based on t-statistics, as for large sample sizes the difference between a t distribution and normal distribution is negligible.
* Calculates power by calculating the non-centrality parameter (NCP) of the chi-square statistic under the alternative.

## Parameters:
* **n**		= Sample size (has to be >= 0).
* **qsq**	= Fraction of trait variance explained by the SNP (has to be between 0 and 1). Denoted as q^2 or q-squared; other studies use h^2 too, to indicate it's similarity with heritability.
* **beta**	= Effect size of the SNP on the trait, in SD units. Only square of beta is used, so result is symmetric in +ve and -ve beta values. For simplicity, preferable to have all beta to be +ve.
* **maf**	= Minor allele frequency of the SNP (between 0 and 0.5).
* **het**	= Heterozygous genotype frequency. Equal to 2\*maf\*(1-maf) under Hardy-Weinberg Equilibrium. (has to be between 0 and 1, usually between 0 and 0.5).
* **pval**	= P-value threshold for significance. Common values are **5E-8** (for genome-wide significance) or **1E-5** (for suggestive significance). (has to be between 0 and 1, >0 and <1).

## Use as: 
```
source("power_calc_functions.R")
```

#### Tested under R version 3.4.1.

## Functions:

### 1. power_n_hsq
Calculates power based on a range of values of n and qsq. pval has to be fixed.
* Parameters as defined above.
* n and hsq can be vectors.
* pval a single numeric value. Default is the common genome-wide significance threshold 5E-8.
* Returns pow, a matrix containing power values, with n along rows and qsq along columns.
* Can contain NA if NCP is negative, NA or infinite, which may happen e.g. if qsq=1.

**example call:**
```
pow = power_n_hsq(n = (1:5)*1000, qsq = (1:10)/100, pval=5E-8)
```

### 2. power_beta_het
Calculates power based on a range of values of beta and het. n and pval has to be fixed.
* Parameters as defined above.
* beta and het can be vectors.
* n and pval has to be single numeric values.
* Default for pval is the common genome-wide significance threshold 5E-8.
* Returns pow, a matrix containing power values, with beta along rows and het along columns.
* Can contain NA if NCP is negative, NA or infinite, which may happen e.g. if beta is very large.

**example call:**
```
pow = power_beta_het(beta = (5:10)/50, het = (1:5)/10, n = 5000, pval=5E-8)
```

### 3. power_beta_maf
Calculates power based on a range of values of beta and maf. n and pval has to be fixed. Could be used when HWE assumption holds.
* Parameters as defined above.
* beta and maf can be vectors.
* n and pval has to be single numeric values.
* Default for pval is the common genome-wide significance threshold 5E-8.
* Returns pow, a matrix containing power values, with beta along rows and maf along columns.
* Can contain NA if NCP is negative, NA or infinite, which may happen e.g. if beta is very large.

**example call:**
```
pow = power_beta_maf(beta = (5:10)/50, maf = (1:5)/10, n = 5000, pval=5E-8)
```

### 4. power_plot
Plots a power matrix as heatmap
* Need to provide X and Y axis labels in addition to the matrix.
* Uses ggplot2 to plot.
* Also requires packages reshape2 and RColorBrewer.

**example call** (includes saving image as pdf) :
```
pdf("power_heatmap.pdf",width=10,height=9); power_plot(pow, "n", "q-squared"); dev.off();
```

### 5. is.scalar
Accessory function: Check whether the input is numeric scalar.
* Used by other functions for power calculation in checking input consistency.
* The user doesn't need to use this function directly.
