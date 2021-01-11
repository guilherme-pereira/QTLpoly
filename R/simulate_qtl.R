#' Simulations of multiple QTL
#'
#' Simulate new phenotypes with a given number of QTL and creates new object with the same structure of class \code{qtlpoly.data} from an existing genetic map.
#'
#' @param data an object of class \code{qtlpoly.data}.
#'
#' @param mu simulated phenotype mean, e.g. 0 (default).
#'
#' @param h2.qtl vector with QTL heritabilities, e.g. \code{c(0.3, 0.2, 0.1)} for three QTL (default); if \code{NULL}, only error is simulated.
#'
#' @param var.error simulated error variance, e.g. 1 (default).
#'
#' @param linked if \code{TRUE} (default), at least two QTL will be linked; if \code{FALSE}, QTL will be randomly assigned along the genetic map. Linkage is defined by a genetic distance smaller than the selected \code{w.size}.
#'
#' @param n.sim number of simulations, e.g. 1000 (default).
#'
#' @param missing if \code{TRUE} (default), phenotypes are simulated with the same number of missing data observed in \code{data$pheno}.
#'
#' @param w.size the window size (in centiMorgans) between two (linked) QTL, e.g. 20 (default).
#'
#' @param seed integer for the \code{set.seed()} function.
#'
#' @param verbose if \code{TRUE} (default), current progress is shown; if \code{FALSE}, no output is produced.
#'
#' @param x an object of class \code{qtlpoly.sim} to be printed.
#'
#' @param detailed if \code{TRUE}, detailed information on linkage groups and phenotypes in shown; if \code{FALSE}, no details are printed.
#'
#' @return An object of class \code{qtlpoly.sim} which contains a list of \code{results} with the same structure of class \code{qtlpoly.data}.
#'
#' @seealso \code{\link[qtlpoly]{read_data}}
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
#'   # simulate new phenotypes
#'   sim.dat <- simulate_qtl(data = data, n.sim = 1000)
#'   }
#'
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2020) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{Genetics} 215 (3): 579-595. \url{http://doi.org/10.1534/genetics.120.303080}.
#'     
#' @export simulate_qtl
#'
#' @import MASS

simulate_qtl <- function(data, mu = 0, h2.qtl = c(0.3, 0.2, 0.1), var.error = 1, linked = FALSE, n.sim = 1000, missing = TRUE, w.size = 20, seed = 123, verbose = TRUE) {
  set.seed(seed)
  progress <- c()
  if(n.sim >= 2) {
    sequence <- unique(round(c(seq(1, n.sim, (n.sim/10))-1), 0))
    progress <- c(sequence[-which(sequence < 2)], n.sim)  
  }
  n.ind <- dim(data$G)[1]
  R <- diag(n.ind) * var.error
  if (verbose) cat("Simulating phenotype number ... 1")
  if(!is.null(h2.qtl)) {
    n.qtl <- length(h2.qtl)
    A <- matrix(rep(h2.qtl, n.qtl), nrow = n.qtl, byrow = FALSE)
    diag(A) <- diag(A) - 1
    b <- -1 * var.error * h2.qtl
    var.qtl <- solve(A, b)
    n.phe <- 0
    sim.mrks <- c()
    sim.pheno <- c()
    sim.error <- c()
    while(n.phe < n.sim) {
      if(!linked) {
        G <- vector("list", n.qtl)
        pheno <- rep(0, n.ind)
        mrks <- sample(1:dim(data$G)[3], size = n.qtl)
        if(all(diff(sort(mrks)) > w.size/data$step)) {
          for(q in 1:n.qtl) {
            G[[q]] <- data$G[,,mrks[q]] * var.qtl[q]
            pheno <- pheno + MASS::mvrnorm(n = 1, mu = rep(mu, n.ind), Sigma = G[[q]])
          }
          error <- MASS::mvrnorm(1, mu = rep(mu, n.ind), Sigma = R)
          pheno <- pheno + error
          n.phe <- n.phe + 1
          sim.mrks <- cbind(sim.mrks, mrks)
          sim.pheno <- cbind(sim.pheno, pheno)
          sim.error <- cbind(sim.error, error)
          if (verbose) if(n.phe %in% progress) cat(" ...", n.phe)
        }
      } else if (linked) {
        G <- vector("list", n.qtl)
        pheno <- rep(0, n.ind)
        lgrs <- sample(1:data$nlgs, size = 1)
        mrks.lgr <- c((data$cum.nmrk[lgrs]+1):data$cum.nmrk[lgrs+1])
        mrks <- sample(mrks.lgr, size = 2)
        mrks <- sample(c(mrks, sample(c(1:dim(data$G)[3])[-mrks.lgr], size = n.qtl-2)))
        if(all(diff(sort(mrks)) > w.size/data$step)) {
          for(q in 1:n.qtl) {
            G[[q]] <- data$G[,,mrks[q]] * var.qtl[q]
            pheno <- pheno + MASS::mvrnorm(n = 1, mu = rep(mu, n.ind), Sigma = G[[q]])
          }
          error <- MASS::mvrnorm(1, mu = rep(mu, n.ind), Sigma = R)
          pheno <- pheno + error
          n.phe <- n.phe + 1
          sim.mrks <- cbind(sim.mrks, mrks)
          sim.pheno <- cbind(sim.pheno, pheno)
          sim.error <- cbind(sim.error, error)
          if (verbose) if(n.phe %in% progress) cat(" ...", n.phe)
        }
      }
    }
    colnames(sim.mrks) <- paste("T", formatC(1:n.phe, width = nchar(n.sim), flag = "0"), sep="")
    colnames(sim.pheno) <- paste(paste("T", formatC(1:n.phe, width = nchar(n.sim), flag = "0"), sep=""), apply(sim.mrks, 2, paste, sep="", collapse="_"), sep=".")
    colnames(sim.error) <- paste("E", formatC(1:n.phe, width = nchar(n.sim), flag = "0"), sep="")
    rownames(sim.error) <- rownames(sim.pheno) <- data$ind.names
  } else {
    sim.pheno <- sim.error <- t(MASS::mvrnorm(n.sim, mu = rep(mu, n.ind), Sigma = R))
    colnames(sim.pheno) <- colnames(sim.error) <- paste("E", formatC(1:n.sim, width = nchar(n.sim), flag = "0"), sep="")
    rownames(sim.error) <- rownames(sim.pheno) <- data$ind.names
    w.size <- sim.mrks <- var.qtl <- NULL
    # if(verbose) cat(" ...", n.sim)
    if (verbose) cat("", paste0("... ", progress))
  }
  if (missing) {
    n.miss <- numeric(data$nphe)
    for(p in 1:data$nphe) n.miss[p] <- sum(!is.na(data$pheno[,p]))
    n.miss <- data$nind - n.miss
    for(i in 1:1000) {
      n.miss.sample <- sample(n.miss, 1)
      if (n.miss.sample > 0) sim.pheno[sample(c(1:data$nind), n.miss.sample), i] <- NA
    }
  }
  if(verbose) cat(". Done!\n")
  data$pheno <- sim.pheno
  data$nphe <- n.sim
  
  structure(list(data = "data",
                 mu = mu,
                 h2.qtl = h2.qtl,
                 var.qtl = var.qtl,
                 var.error = var.error,
                 n.sim = n.sim,
                 missing = missing,
                 w.size = w.size,
                 seed = seed,
                 sim.mrks = sim.mrks,
                 sim.error = sim.error,
                 results = data),
            class=c("qtlpoly.simul"))
}

#' @rdname simulate_qtl
#' @export

print.qtlpoly.simul <- function(x, detailed = FALSE) {
  x <- x$results
  cat("This is an object of class 'qtlpoly.simul'\n")
  cat("  Ploidy level:       ", x$ploidy, "\n", sep="")
  cat("  No. individuals:    ", x$nind, "\n", sep="")
  cat("  No. linkage groups: ", x$nlgs, "\n", sep="")
  cat("  Step size:          ", x$step, " cM \n", sep="")
  cat("  Map size:           ", round(last(x$cum.size), 2), " cM (", last(x$cum.nmrk), " positions) \n", sep="")
  if(detailed) for(c in 1:x$nlgs) cat("    LG ", c, ": ", round(x$lgs.size[[c]], 2), " cM (", x$lgs.nmrk[[c]], " positions) \n", sep="")
  cat("  No. phenotypes:     ", x$nphe, "\n", sep="")
  if(detailed) for(p in 1:x$nphe) cat("    Trait ", p, ": ", sQuote(colnames(x$pheno)[p]), " (", sum(!is.na(x$pheno[,p])), " individuals) \n", sep="")
}
