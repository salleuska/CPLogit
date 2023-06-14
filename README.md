# Supplementary code for Bayesian linear and logistic regression using Centered Partition processes

Code to reproduce simulation results presented in  
**Centered Partition Processes: Informative Priors for Clustering (with Discussion)**   
by Sally Paganin, Amy H. Herring, Andrew F. Olshan, David B. Dunson  
*Bayesian Analysis*. 16(1), 301-370, (March 2021)
\[ 
[Link](https://projecteuclid.org/journals/bayesian-analysis/volume-16/issue-1/Centered-Partition-Processes-Informative-Priors-for-Clustering-with-Discussion/10.1214/20-BA1197.full) | 
[arXiV preprint](https://arxiv.org/abs/1901.10225)
\]

## Installation

Be sure to have installed the following R packages

```r
	install.packages("Rcpp")
	install.packages("RcppArmadillo")
	install.packages("devtools")
	install.packages("ggplot2")
	install.packages("reshape2")
	install.packages("magrittr")
	install.packages("dplyr")
	install.packages("sdols")
	install.packages("mcclust")
	install.packages("mvtnorm")
```

Download and install the package `CPLogit` using `devtools` 

```r
library(devtools)
devtools::install_github("salleuska/CPLogit", subdir="CPLogit")
```
or if you download the repo
```r
## How to install 
library(devtools)
library(Rcpp)
clean_dll()

## to export Rcpp functions
compileAttributes()
build()
install()
document()
```

## Reproduce plots in Section 3.3

To reproduce code of the prior probability distribution over the set partition space induced by the CP process, use script `Sec3.3_Prior_graphs.R`.

## Reproduce simulation results

Simulate data used in the Sec 5.2 and in Sec 1 of the Supplementary Materials

```bash
Rscript simulation/0_generateData_linearSimulation.R
Rscript simulation/0_generateData_logisticSimulation.R

```
Run models and reproduce plots in Sec 5.2 and in Sec 1 of the Supplementary Materials

```bash
Rscript simulation/1_linearSimulation.R
Rscript simulation/1_logisticSimulation.R
```
