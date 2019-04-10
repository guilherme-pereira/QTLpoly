#' Simulations of multiple QTL
#'
#' Simulate new phenotypes with a given number of QTL and creates new object with the same structure of class \code{qtlpoly.data} from an existing genetic map.
#'
#' @param data an object of class \code{qtlpoly.data}.
#'
#' @param mu simulated phenotype mean, e.g. 0 (default).
#'
#' @param h2.qtl vector with QTL heritabilities, e.g. c(0.3, 0.2, 0.1) for three QTL (default).
#'
#' @param var.error simulated error variance, e.g. 1 (default).
#'
#' @param n.sim number of simulations, e.g. 1000 (default).
#'
#' @param w.size the window size (in centiMorgans) to avoid on either side of QTL already in the model when looking selecting new QTL, e.g. 20 (default).
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
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2019) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{bioRxiv}. \url{doi:}.
#'
#' @export simulate_qtl
#'
#' @import MASS

simulate_qtl <- function(data, mu = 0, h2.qtl = c(0.3, 0.2, 0.1), var.error = 1, n.sim = 1000, w.size = 20, seed = 123, verbose = TRUE) {
  set.seed(seed)
  progress <- seq(0, n.sim, (n.sim/10))
  n.qtl <- length(h2.qtl)
  n.phe <- 0
  n.ind <- dim(data$G)[1]
  sim.mrks <- c()
  sim.pheno <- c()
  A <- matrix(rep(h2.qtl, n.qtl), nrow = n.qtl, byrow = FALSE)
  diag(A) <- diag(A)-1
  b <- -1 * h2.qtl
  var.qtl <- solve(A, b)
  var.phen <- sum(var.qtl, var.error)
  R <- diag(n.ind) * var.error
  if (verbose) cat("  Simulating phenotype number ... 1")
  while(n.phe < n.sim) {
    G <- vector("list", n.qtl)
    pheno <- rep(0, n.ind)
    mrks <- sample(1:dim(data$G)[3], size = n.qtl)
    if(all(diff(sort(mrks)) > w.size)) {
      for(q in 1:n.qtl) {
        G[[q]] <- data$G[,,mrks[q]] * var.qtl[q]
        pheno <- pheno + MASS::mvrnorm(n = 1, mu = rep(mu, n.ind), Sigma = G[[q]])
      }
      pheno <- pheno + MASS::mvrnorm(1, mu = rep(mu, n.ind), Sigma = R)
      n.phe <- n.phe + 1
      sim.mrks <- cbind(sim.mrks, mrks)
      sim.pheno <- cbind(sim.pheno, pheno)
      if (verbose) if(n.phe %in% progress) cat(" ...", n.phe)
    }
  }
  colnames(sim.mrks) <- NULL
  colnames(sim.pheno) <- paste(paste("T", formatC(1:n.phe, width = 4, flag = "0"), sep=""), apply(sim.mrks, 2, paste, sep="", collapse="_"), sep=".")
  if(verbose) cat(". Done!\n")
  data$pheno <- data.frame(sim.pheno)
  data$nphe <- n.phe

  structure(list(data = deparse(substitute(data)),
                 mu = mu,
                 h2.qtl = h2.qtl,
                 var.qtl = var.qtl,
                 var.error = var.error,
                 n.sim = n.sim,
                 w.size = w.size,
                 seed = seed,
                 sim.mrks = sim.mrks,
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
