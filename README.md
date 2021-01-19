# CPLogit

Code to reproduce results presented in [Centered Partition Processes: Informative Priors for Clustering](https://projecteuclid.org/euclid.ba/1581584439#abstract).

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
```

Download and install the package `CPLogit` using `devtools` and `Rcppp` packages

```r
library(devtools)

clean_dll()
build()
install()
```

<!-- 
## How to install 
library(devtools)
library(Rcpp)
clean_dll()
## to export Rcpp functiond
compileAttributes()
build()
install()
document()
 -->

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
