# R Functions to calculate power of GWAS studies for a single associated SNP, under various parameters.
# Suitable for classical (i.e. single-SNP single-trait) GWAS studies using linear regression models, i.e for quantitative traits.
# Uses formulae for power calculation presented in Appendix A of
# Visscher PM, Wray NR, Zhang Q, et al. 10 Years of GWAS Discovery: Biology, Function, and Translation. Am J Hum Genet 2017;101(1):5-22. doi: 10.1016/j.ajhg.2017.06.005.
# Assumes that additional covariates, if any, are uncorrelated to the SNP.
# It is common to include genetic PCs as covariates in GWAS to adjust for stratification or admixture;
# in such cases, if the SNP is severely stratified, these formulae will not be applicable, and a simulation-based approach might be preferrable.
# Uses a chi-square statistic, equivalent to using the z-statistic in a regression (the chi-square statistic is the square of the z-statistic).
# The results should be very similar to using tests based on t-statistics, as for large sample sizes the difference between a t distribution and normal distribution is negligible.
# Calculates power by calculating the non-centrality parameter (NCP) of the chi-square statistic under the alternative.

# Parameters:
# n		= Sample size (has to be >= 0)
# qsq	= Fraction of trait variance explained by the SNP (has to be between 0 and 1). Denoted as q^2 or q-squared; other studies use h^2 too, to indicate it's similarity with heritability.
# beta	= Effect size of the SNP on the trait, in SD units. Only square of beta is used, so result is symmetric in +ve and -ve beta values. For simplicity, preferable to have all beta to be +ve.
# maf	= Minor allele frequency of the SNP (between 0 and 0.5).
# het	= Heterozygous genotype frequency. Equal to 2*maf*(1-maf) under Hardy-Weinberg Equilibrium. (has to be between 0 and 1, usually between 0 and 0.5).
# pval	= P-value threshold for significance. Common values are 5E-8 (for genome-wide significance) or 1E-5 (for suggestive significance). (has to be between 0 and 1, >0 and <1).

# Use as: source("power_calc_functions.R")

# Tested under R version 3.4.1.


is.scalar <- function(v) { is.numeric(v) && length(v)==1; }
# Accessory function: Check whether the input is numeric scalar.
# Used by other functions for power calculation in checking input consistency.
# The user doesn't need to use this function directly.


power_n_hsq <- function(n, qsq, pval = 5E-8) {
# Calculates power based on a range of values of n and qsq. pval has to be fixed.
# Parameters as defined above.
# n and hsq can be vectors.
# pval a single numeric value. Default is the common genome-wide significance threshold 5E-8.
# Returns pow, a matrix containing power values, with n along rows and qsq along columns.
# Can contain NA if NCP is negative, NA or infinite, which may happen e.g. if qsq=1.
# Example call: pow = power_n_hsq(n = (1:5)*1000, qsq = (1:10)/100, pval=5E-8)

if (missing(n)) { stop("Parameter n not found."); }
if (missing(qsq)) { stop("Parameter qsq not found."); }
if ( !( is.vector(n) && is.numeric(n) ) ) { stop("Parameter n not a numeric vector."); }
if ( !( is.vector(qsq) && is.numeric(qsq) ) ) { stop("Parameter qsq not a numeric vector."); }
if ( !( is.scalar(pval) ) ) { stop("Parameter pval not a numeric scalar."); }
if ( any( n<0 | !is.finite(n) ) ) { stop("Parameter n has unacceptable values."); }
if ( any( qsq<0 | qsq>1 | !is.finite(qsq) ) ) { stop("Parameter qsq has unacceptable values."); }
if ( pval<=0 | pval>=1 | !is.finite(pval) ) { stop("Parameter pval has unacceptable value."); }

th=qchisq(pval,df=1,lower.tail=F); # Significance threshold for chi-square, corresponding to P-value threshold

pow=matrix(NA,nrow=length(n),ncol=length(qsq)); # Pre-allocalte matrix to store power values
rownames(pow)=n; colnames(pow)=qsq; # Assign row and column names

for (i in 1:length(n)) {
	for (j in 1:length(qsq)) {
		ncp=n[i]*qsq[j]/(1-qsq[j]); # Calculate NCP parameter
		if (is.finite(ncp) && ncp>=0) { pow[i,j]=pchisq(th,df=1,lower.tail=F,ncp=ncp); } # Calculate power when NCP is finite and >=0
}}

return(pow); # Return calculated power matrix
}




power_beta_het <- function(beta, het, n, pval = 5E-8) {
# Calculates power based on a range of values of beta and het. n and pval has to be fixed.
# Parameters as defined above.
# beta and het can be vectors.
# n and pval has to be single numeric values.
# Default for pval is the common genome-wide significance threshold 5E-8.
# Returns pow, a matrix containing power values, with beta along rows and het along columns.
# Can contain NA if NCP is negative, NA or infinite, which may happen e.g. if beta is very large.
# Example call: pow = power_beta_het(beta = (5:10)/50, het = (1:5)/10, n = 5000, pval=5E-8)

if (missing(beta)) { stop("Parameter beta not found."); }
if (missing(het)) { stop("Parameter het not found."); }
if (missing(n)) { stop("Parameter n not found."); }
if ( !( is.vector(beta) && is.numeric(beta) ) ) { stop("Parameter beta not a numeric vector."); }
if ( !( is.vector(het) && is.numeric(het) ) ) { stop("Parameter het not a numeric vector."); }
if ( !( is.scalar(n) ) ) { stop("Parameter n not a numeric scalar."); }
if ( !( is.scalar(pval) ) ) { stop("Parameter pval not a numeric scalar."); }
if ( any( !is.finite(beta) ) ) { stop("Parameter beta has unacceptable values."); }
if ( any( het<0 | het>1 | !is.finite(het) ) ) { stop("Parameter het has unacceptable values."); }
if ( n<0 | !is.finite(n) ) { stop("Parameter n has unacceptable value."); }
if ( pval<=0 | pval>=1 | !is.finite(pval) ) { stop("Parameter pval has unacceptable value."); }

th=qchisq(pval,df=1,lower.tail=F); # Significance threshold for chi-square, corresponding to P-value threshold

pow=matrix(NA,nrow=length(beta),ncol=length(het)); # Pre-allocalte matrix to store power values
rownames(pow)=beta; colnames(pow)=het; # Assign row and column names

for (i in 1:length(beta)) {
	for (j in 1:length(het)) {
		q2=het[j]*(beta[i]^2); ncp=n*q2/(1-q2); # Calculate qsq and then NCP parameter
		if (is.finite(ncp) && ncp>=0) { pow[i,j]=pchisq(th,df=1,lower.tail=F,ncp=ncp); } # Calculate power when NCP is finite and >=0
}}

return(pow); # Return calculated power matrix
}






power_beta_maf <- function(beta, maf, n, pval = 5E-8) {
# Calculates power based on a range of values of beta and maf. n and pval has to be fixed. Could be used when HWE assumption holds.
# Parameters as defined above.
# beta and maf can be vectors.
# n and pval has to be single numeric values.
# Default for pval is the common genome-wide significance threshold 5E-8.
# Returns pow, a matrix containing power values, with beta along rows and maf along columns.
# Can contain NA if NCP is negative, NA or infinite, which may happen e.g. if beta is very large.
# Example call: pow = power_beta_maf(beta = (5:10)/50, maf = (1:5)/10, n = 5000, pval=5E-8)

if (missing(beta)) { stop("Parameter beta not found."); }
if (missing(maf)) { stop("Parameter maf not found."); }
if (missing(n)) { stop("Parameter n not found."); }
if ( !( is.vector(beta) && is.numeric(beta) ) ) { stop("Parameter beta not a numeric vector."); }
if ( !( is.vector(maf) && is.numeric(maf) ) ) { stop("Parameter maf not a numeric vector."); }
if ( !( is.scalar(n) ) ) { stop("Parameter n not a numeric scalar."); }
if ( !( is.scalar(pval) ) ) { stop("Parameter pval not a numeric scalar."); }
if ( any( !is.finite(beta) ) ) { stop("Parameter beta has unacceptable values."); }
if ( any( maf<0 | maf>0.5 | !is.finite(maf) ) ) { stop("Parameter maf has unacceptable values."); }
if ( n<0 | !is.finite(n) ) { stop("Parameter n has unacceptable value."); }
if ( pval<=0 | pval>=1 | !is.finite(pval) ) { stop("Parameter pval has unacceptable value."); }

th=qchisq(pval,df=1,lower.tail=F); # Significance threshold for chi-square, corresponding to P-value threshold

pow=matrix(NA,nrow=length(beta),ncol=length(maf)); # Pre-allocalte matrix to store power values
rownames(pow)=beta; colnames(pow)=maf; # Assign row and column names

for (i in 1:length(beta)) {
	for (j in 1:length(maf)) {
		q2=2*maf[j]*(1-maf[j])*(beta[i]^2); ncp=n*q2/(1-q2); # Calculate qsq and then NCP parameter
		if (is.finite(ncp) && ncp>=0) { pow[i,j]=pchisq(th,df=1,lower.tail=F,ncp=ncp); } # Calculate power when NCP is finite and >=0
}}

return(pow); # Return calculated power matrix
}





power_plot <- function(pow, xlabel, ylabel) {
# Plots a power matrix as heatmap
# Need to provide X and Y axis labels in addition to the matrix.
# Uses ggplot2 to plot.
# Also requires packages reshape2 and RColorBrewer.
# Example call (includes saving image as pdf) : pdf("power_heatmap.pdf",width=10,height=9); power_plot(pow, "n", "q-squared"); dev.off();

if (missing(pow)) { stop("Power matrix not supplied not found."); }
if (missing(xlabel)) { stop("X axis label not found."); }
if (missing(ylabel)) { stop("Y axis label not found."); }

require("ggplot2")
require("reshape2")
require("RColorBrewer")

powm=melt(pow); # Convert to data frame format, suitable for plotting. melt is from reshape2
colnames(powm)=c("x","y","Power"); # Set column names

cols=brewer.pal(9,"YlOrRd"); # Use a yellow-red color gradient

# make the plot
ggplot(data=powm, aes(x=x,y=y,z=Power)) + geom_raster(aes(fill=Power),interpolate=TRUE) + scale_fill_gradientn(colours=cols) + theme_classic() + xlab(xlabel) + ylab(ylabel);

}



