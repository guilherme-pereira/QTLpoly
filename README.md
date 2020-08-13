[![Build Status](https://travis-ci.org/guilherme-pereira/QTLpoly.svg?branch=master)](https://travis-ci.org/guilherme-pereira/QTLpoly) [![Development](https://img.shields.io/badge/development-active-blue.svg)](https://img.shields.io/badge/development-active-blue.svg)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# QTLpoly

The R package `qtlpoly` (v. 0.1-0) is an under development software to map quantitative trait loci (QTL) in full-sib families of outcrossing autopolyploid species based on a random-effect multiple QTL model (Pereira et al. 2019). 

In order to do so, you will need a genetic map from which conditional probabilities of putative QTL can be computed. We recommend [`mappoly`](https://github.com/mmollina/MAPpoly), a hidden Markov model-based R package to construct genetic maps in autopolyploids (Mollinari and Garcia 2019).

Variance components associated with putative QTL are tested using score statistics from the R package `varComp` (v. 0.2-0) (Qu et al. 2013). Final models are fitted using residual maximum likelihood (REML) from the R package `sommer` (v. 3.6) (Covarrubias-Pazaran 2016). Plots for visualizing the results are based on `ggplot2` (v. 3.3.2) (Wickham 2016). 

## Install `qtlpoly` package

As mentioned, the package `qtlpoly` depends on a couple of functions from `sommer` (v. 3.6) and `varComp` (v. 0.2-0). [`varComp`](https://cran.r-project.org/src/contrib/Archive/varComp/) and its dependency [SPA3G](https://cran.r-project.org/src/contrib/Archive/SPA3G/) have been archived from CRAN, while `sommer` has been [constantly updated](https://cran.r-project.org/web/packages/sommer/index.html). In order to avoid conflit with updates in functions and object structures, we decided to stick with an earlier `sommer` version (v. 3.6), stable to our needs, so you will need to downgrade if you have installed the most recent version of `sommer`.

`qtlpoly` package is available here on [GitHub](https://github.com/guilherme-pereira/qtlpoly). You can install all needed packages within R using the functions from the R package `devtools`:

```r
> install.packages("devtools")
> devtools::install_url("https://cran.r-project.org/src/contrib/Archive/SPA3G/SPA3G_1.0.tar.gz")
> devtools::install_url("https://cran.r-project.org/src/contrib/Archive/varComp/varComp_0.2-0.tar.gz")
> devtools::install_version("sommer", version = "3.6", repos = "http://cran.us.r-project.org")
> devtools::install_github("guilherme-pereira/qtlpoly", upgrade = FALSE) 
```

Do not update `sommer` to the most recent version even if R asks you to do so while installing `qtlpoly`. 

## Documents

Tutorials as well as simulated and real data set analyses will be listed here opportunately in order to help users to get familiar with the software and allow them to perform their own analyses:

1. [Tutorial on Multiple QTL Mapping in Autopolyploids with QTLpoly](https://guilherme-pereira.github.io/QTLpoly/1-tutorial)
2. Score-based resampling method to assess genome-wide significance

## Acknowledgments

This package has been developed as part of the [Genomic Tools for Sweetpotato Improvement](https://sweetpotatogenomics.cals.ncsu.edu/) (GT4SP) project, funded by [Bill \& Melinda Gates Foundation](https://www.gatesfoundation.org/).

## References

Covarrubias-Pazaran G. 2016. “Genome-assisted prediction of quantitative traits using the R package sommer.” PLoS ONE 11 (6): 1-15. [doi:10.1371/journal.pone.0156744](https://doi.org/10.1371/journal.pone.0156744).

Mollinari M, Garcia AAF. 2019. “Linkage analysis and haplotype phasing in experimental autopolyploid populations with high ploidy level using hidden Markov models.” G3: Genes, Genomes, Genetics 9 (10): 3297-3314. [doi:10.1534/g3.119.400378](https://doi.org/10.1534/g3.119.400378).

Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB. 2020. “Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population.” Genetics 215 (3): 579-595. [doi:10.1534/genetics.120.303080](https://doi.org/10.1534/genetics.120.303080).

Qu L, Guennel T, Marshall SL. 2013. “Linear score tests for variance components in linear mixed models and applications to genetic association studies.” Biometrics 69 (4): 883-892. [doi:10.1111/biom.12095](https://doi.org/10.1111/biom.12095).

Wickham H. 2016. “ggplot2: Elegant Graphics for Data Analysis.” Springer. [doi:10.1007/978-0-387-98141-3](https://www.springer.com/gp/book/9780387981413).
