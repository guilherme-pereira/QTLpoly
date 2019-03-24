---
output:
  html_document:
    highlight: tango
    theme: united
---

# QTLpoly

The R package `qtlpoly` (v. 0.1-0) is an under development software to map quantitative trait loci (QTL) in full-sib families of outcrossing autopolyploid species based on a random-effect multiple QTL model (Pereira et al. 2019).

Variance components associated with putative QTL  are tested using score statistics from the R package `varComp` (v. 0.2-0) (Qu et al. 2013). Final models are fitted using residual maximum likelihood (REML) from the R package `sommer` (v. 3.6) (Covarrubias-Pazaran 2016). Plots for visualizing the results are based on `ggplot2` (v. 3.1.0) (Wickham 2016). 

## Install and load the `qtlpoly` package and data

As mentioned, the package `qtlpoly` depends on a couple of functions from `sommer` (v. 3.6) and `varComp` (v. 0.2-0). `varComp` has been [archived from CRAN](https://cran.r-project.org/src/contrib/Archive/varComp/), while `sommer` has been constantly updated (currently, in its v. 3.7). In order to avoid conflit with updates in functions and object structures, we decided to stick with an earlier `sommer` version (v. 3.6), so you will need to downgrade if you have installed the most recent version of `sommer`.

`qtlpoly` package is available at [GitHub](https://github.com/guilherme-pereira/qtlpoly). You can install all needed packages within R using the functions from the R package `devtools`:

```r
> install.packages("devtools")
> devtools::install_url("https://cran.r-project.org/src/contrib/Archive/varComp/varComp_0.2-0.tar.gz")
> devtools::install_version("sommer", version = "3.6", repos = "http://cran.us.r-project.org")
> devtools::install_github("guilherme-pereira/qtlpoly") 
```

## Tutorial

Simulated and real data set analyses will be posted here accordingly in order to help users to get familiar with the software and allow them to perform their own analyses:

1. [Tutorial on Multiple QTL Mapping in Autopolyploids with QTLpoly](http)

## Acknowledgments

This package has been developed as part of the [Genomic Tools for Sweetpotato Improvement](https://sweetpotatogenomics.cals.ncsu.edu/) (GT4SP) project, funded by [Bill \& Melinda Gates Foundation](https://www.gatesfoundation.org/).

## References

Covarrubias-Pazaran, Giovanny. 2016. “Genome-assisted prediction of quantitative traits using the R package sommer.” PLoS ONE 11 (6): 1–15. [doi:10.1371/journal.pone.0156744](doi:10.1371/journal.pone.0156744).

Pereira, Guilherme Silva, Dorcus C Gemenet, Marcelo Mollinari, Bode A Olukolu, Joshua C Wood, Veronica Mosquera, Wolfgang J Gruneberg, et al. 2019. “Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population.” Submitted, 1–15.

Qu, Long, Tobias Guennel, and Scott L. Marshall. 2013. “Linear score tests for variance components in linear mixed models and applications to genetic association studies.” Biometrics 69 (4): 883–92. [doi:10.1111/biom.12095](doi:10.1111/biom.12095).

Wickham, Hadley. 2016. ggplot2: Elegant Graphics for Data Analysis. Springer.
