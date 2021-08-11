#' @title Random-effect multiple interval mapping (REMIM) model score-based resampling method
#'
#' @description Stores P-values from maximum score statistics for a number of resampling of given linkage map.
#'
#' @param data an object of class \code{qtlpoly.data}.
#'
#' @param pheno.col a numeric vector with the phenotype columns to be analyzed; if \code{NULL} (default), all phenotypes from \code{'data'} will be included.
#'
#' @param n.sim a number of simulations, e.g. 1000 (default).
#'
#' @param n.clusters a number of parallel processes to spawn.
#'
#' @param seed an integer for the \code{set.seed()} function; if \code{NULL}, no reproducible seeds are set.
#'
#' @param verbose if \code{TRUE} (default), current progress is shown; if \code{FALSE}, no output is produced.
#'
#' @param x an object of class \code{qtlpoly.perm} to be printed or plotted.
#'
#' @param probs a vector of probability values in [0, 1] representing the quantiles, e.g. c(0.90, 0.95) for the 90\% and 95\% quantiles.
#'
#' @return An object of class \code{qtlpoly.perm} which contains a list of \code{results} for each trait with the maximum LOD score per permutation.
#'
#' @return LOD score thresholds for given quantiles for each trait.
#'
#' @return A \pkg{ggplot2} histogram with the distribution of ordered maximum LOD scores and thresholds for given quantiles for each trait.
#'
#' @seealso \code{\link[qtlpoly]{feim}}
#'
#' @examples
#'   \dontrun{
#'   # load raw data
#'   data(maps)
#'   data(pheno)
#'
#'   # estimate conditional probabilities using mappoly package
#'   library(mappoly)
#'   genoprob <- lapply(maps, calc_genoprob)
#'
#'   # prepare data
#'   data <- read_data(ploidy = 6, geno.prob = genoprob, pheno = pheno, step = 1)
#'
#'   # perform resampling
#'   score.rsmp <- resampling(data = data, n.sim = 1000, n.clusters = 4)
#'   }
#'
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}
#'
#' @references
#'     Zou F, Fine JP, Hu J, Lin DY (2004) An efficient resampling method for assessing genome-wide statistical significance in mapping quantitative trait loci, \emph{Genetics} 168 (4): 2307-2316. \url{https://doi.org/10.1534/genetics.104.031427}.
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2020) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{Genetics} 215 (3): 579-595. \url{http://doi.org/10.1534/genetics.120.303080}.
#'
#' @export resampling
#' @import parallel

resampling <- function(data, n.sim = 1000, alpha = c(0.20, 0.05), n.clusters = NULL, seed = 123, verbose = TRUE) {
  
  data.sim <- simulate_qtl(data = data, mu = 0, h2.qtl = NULL, var.error = 1, n.sim = n.sim, missing = TRUE, seed = seed) 
  score.null <- null_model(data = data.sim$results, n.clusters = n.clusters, plot = NULL)
  min.pvl <- numeric(length(score.null$results))
  for(p in 1:length(score.null$results)) {
    min.pvl[p] <- score.null$results[[p]]$pval[which.max(score.null$results[[p]]$stat)]
  }
  quantile(sort(min.pvl), alpha)

  structure(list(data=deparse(substitute(data)),
                 data.sim=data.sim,
                 n.sim=n.sim,
                 alpha=alpha,
                 seed=seed,
                 score.null=score.null,
                 min.pvl=min.pvl
  ),
  class=c("qtlpoly.perm"))
}

#' @rdname resampling
#' @export

print.qtlpoly.rsmp <- function(x, pheno.col=NULL, probs=c(0.90, 0.95)) {
  if(any(class(x) == "qtlpoly.rsmp")) cat("This is an object of class 'qtlpoly.rsmp'\n")
  quantile(sort(min.pvl), alpha)
}

#' @rdname resampling
#' @import ggplot2
#' @export

plot.qtlpoly.rsmp <- function(x, alpha=c(0.20, 0.05)) {
  if(!is.null(x$results[[p]])) {
    data <- x$min.pvl
    plot <- ggplot(data, aes(LOD)) +
      geom_histogram(bins = 15) +
      geom_vline(xintercept = quantile(data$LOD, probs), linetype="dashed") +
      annotate("text", x=quantile(data$LOD, probs), y=length(data$LOD)/6, angle=90, vjust=-0.1, hjust=0,
               label=paste(names(quantile(data$LOD, probs)), round(quantile(data$LOD, probs), digits = 2), sep=" = ")) +
      labs(title=names(x$results)[p], y="Count", x="Maximum LOD scores") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5), title=element_text(face="bold"))
    print(plot)
  }
}

