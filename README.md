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

**Example call:**
```
pow = power_n_hsq(n = (1:5)*1000, qsq = (1:10)/100, pval=5E-8)
```
**Example output:**
```
> pow
           0.01      0.02      0.03      0.04      0.05      0.06      0.07      0.08      0.09       0.1
1000 0.01151002 0.1752110 0.5437832 0.8422292 0.9643414 0.9944263 0.9993689 0.9999464 0.9999965 0.9999998
2000 0.16937331 0.8257378 0.9921005 0.9998822 0.9999992 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
3000 0.52133658 0.9911850 0.9999855 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
4000 0.81729560 0.9998307 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
5000 0.95107629 0.9999983 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
```

### 2. power_beta_het
Calculates power based on a range of values of beta and het. n and pval has to be fixed.
* Parameters as defined above.
* beta and het can be vectors.
* n and pval has to be single numeric values.
* Default for pval is the common genome-wide significance threshold 5E-8.
* Returns pow, a matrix containing power values, with beta along rows and het along columns.
* Can contain NA if NCP is negative, NA or infinite, which may happen e.g. if beta is very large.

**Example call:**
```
pow = power_beta_het(beta = (5:10)/50, het = (1:5)/10, n = 5000, pval=5E-8)
```
**Example output:**
```
> pow
              0.1        0.2       0.3       0.4       0.5
0.1  0.0006542167 0.01113106 0.0579168 0.1659726 0.3304165
0.12 0.0028366019 0.04935898 0.2136970 0.4724055 0.7157991
0.14 0.0102316499 0.15495311 0.4947653 0.7979632 0.9433170
0.16 0.0308104754 0.35246887 0.7791710 0.9591605 0.9953441
0.18 0.0778254130 0.60230655 0.9399304 0.9959898 0.9998541
0.2  0.1659725860 0.81559267 0.9903982 0.9998183 0.9999983
```

### 3. power_beta_maf
Calculates power based on a range of values of beta and maf. n and pval has to be fixed. Could be used when HWE assumption holds.
* Parameters as defined above.
* beta and maf can be vectors.
* n and pval has to be single numeric values.
* Default for pval is the common genome-wide significance threshold 5E-8.
* Returns pow, a matrix containing power values, with beta along rows and maf along columns.
* Can contain NA if NCP is negative, NA or infinite, which may happen e.g. if beta is very large.

**Example call:**
```
pow = power_beta_maf(beta = (5:10)/50, maf = (1:5)/10, n = 5000, pval=5E-8)
```
**Example output:**
```
> pow
             0.1        0.2       0.3       0.4       0.5
0.1  0.007170496 0.07424369 0.1951479 0.2944155 0.3304165
0.12 0.032399975 0.26101711 0.5257082 0.6728830 0.7157991
0.14 0.106771502 0.56605303 0.8391246 0.9250314 0.9433170
0.16 0.261017108 0.83523627 0.9726057 0.9925678 0.9953441
0.18 0.485842586 0.96297523 0.9978416 0.9997051 0.9998541
0.2  0.715799127 0.99534406 0.9999252 0.9999955 0.9999983
```

### 4. power_plot
Plots a power matrix as heatmap.
* Need to provide X and Y axis labels in addition to the matrix.
* Uses ggplot2 to plot.
* Also requires packages reshape2 and RColorBrewer.

**Example call** (includes saving image as pdf) :
```
pdf("power_heatmap.pdf",width=10,height=9); power_plot(pow, "n", "q-squared"); dev.off();
```

**Example output** (using pow matrix from step 1): *see attached* **power_heatmap.pdf**


### 5. is.scalar
Accessory function: Check whether the input is numeric scalar.
* Used by other functions for power calculation in checking input consistency.
* The user doesn't need to use this function directly.
