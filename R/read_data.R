#' Read data
#'
#' Reads files in specific formats and creates a \code{qtlpoly.data} object to be used in subsequent analyses.
#'
#' @param ploidy a numeric value of ploidy level of the cross.
#'
#' @param geno.prob an object of class \code{mappoly.genoprob} from \pkg{mappoly}.
#'
#' @param pheno a data frame of phenotypes (columns) with individual names (rows) identical to individual names in \code{geno.prob} object.
#'
#' @param step a numeric value of step size (in centiMorgans) where tests will be performed, e.g. 1 (default); if \code{NULL}, tests will be performed at every marker.
#'
#' @param x an object of class \code{qtlpoly.data} to be printed.
#'
#' @param detailed if \code{TRUE}, detailed information on linkage groups and phenotypes in shown; if \code{FALSE}, no details are printed.
#'
#' @return An object of class \code{qtlpoly.data} which is a list containing the following components:
#'
#'     \item{ploidy}{a scalar with ploidy level.}
#'     \item{nlgs}{a scalar with the number of linkage groups.}
#'     \item{nind}{a scalar with the number of individuals.}
#'     \item{nmrk}{a scalar with the number of marker positions.}
#'     \item{nphe}{a scalar with the number of phenotypes.}
#'     \item{lgs.size}{a vector with linkage group sizes.}
#'     \item{cum.size}{a vector with cumulative linkage group sizes.}
#'     \item{lgs.nmrk}{a vector with number of marker positions per linkage group.}
#'     \item{cum.nmrk}{a vector with cumulative number of marker positions per linkage group.}
#'     \item{lgs}{a list with selected marker positions per linkage group.}
#'     \item{lgs.all}{a list with all marker positions per linkage group.}
#'     \item{step}{a scalar with the step size.}
#'     \item{pheno}{a data frame with phenotypes.}
#'     \item{G}{a list of relationship matrices for each marker position.}
#'     \item{Z}{a list of conditional probability matrices for each marker position for genotypes.}
#'     \item{X}{a list of conditional probability matrices for each marker position for alleles.}
#'     \item{Pi}{a matrix of identical-by-descent shared alleles among genotypes.}
#'
#' @seealso \code{\link[mappoly]{read_data}}, \code{\link[qtlpoly]{maps}}, \code{\link[qtlpoly]{pheno}}
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
#'   }
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2019) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{bioRxiv}. \url{doi:}.
#'
#' @export read_data
#' @importFrom abind abind

read_data <- function(ploidy = 6, geno.prob = genoprob, pheno = pheno, step = 1) {

  if (length(which(rownames(pheno) %in% dimnames(geno.prob[[1]]$probs)[[3]])) == 0)
    stop("Names of individuals from both genotype and phenotype data do not match. Please, check data \n")

  if(is.null(step)) step <- 1e-10
  nlgs <- length(geno.prob)
  probs <- vector("list", nlgs)
  lgs <- vector("list", nlgs)
  lgs.all <- vector("list", nlgs)
  lgs.size <- numeric(nlgs)
  lgs.nmrk <- numeric(nlgs)
  for(c in 1:nlgs) {
    lgs[[c]] <- unique(c(geno.prob[[c]]$map[!duplicated(floor(geno.prob[[c]]$map/step)*step)], last(geno.prob[[c]]$map)))
    lgs.all[[c]] <- geno.prob[[c]]$map
    lgs.size[c] <- last(lgs[[c]])
    names(lgs[[c]]) <- names(geno.prob[[c]]$map)[which(geno.prob[[c]]$map %in% lgs[[c]])]
    probs[[c]] <- geno.prob[[c]]$probs[,which(geno.prob[[c]]$map %in% lgs[[c]]),]
    lgs.nmrk[c] <- dim(probs[[c]])[2]
  }
  cum.size <- c(0, cumsum(lgs.size))
  cum.nmrk <- c(0, cumsum(lgs.nmrk))

  Z <- abind::abind(probs, along = 2); dim(Z)
  indnames <- dimnames(Z)[[3]]
  mrknames <- dimnames(Z)[[2]]
  nmrk <- dim(Z)[2]
  nind <- dim(Z)[3]

  Palleles <- letters[1:ploidy]
  Pgametes <- lapply(combn(Palleles, ploidy/2, simplify = FALSE), paste, collapse="")
  Qalleles <- letters[(ploidy+1):(2*ploidy)]
  Qgametes <- lapply(combn(Qalleles, ploidy/2, simplify = FALSE), paste, collapse="")
  genotypes <- as.vector(t(outer(Pgametes, Qgametes, paste, sep="")))
  sibs <- sapply( genotypes, FUN=function(x) paste(x, genotypes, sep="") )
  Pi <- matrix(data = NA, nrow = length(Pgametes)^2, ncol = length(Pgametes)^2)
  for(i in 1:ncol(sibs)) {
    for(j in 1:nrow(sibs)) {
      Pi[i,j] <- ((2*ploidy)-length(unique(strsplit(sibs[i,j], "")[[1]])))/ploidy
    }
  }
  colnames(Pi) <- rownames(Pi) <- as.vector(t(outer(Pgametes, Qgametes, paste, sep=":")))

  G <- array(data = NA, dim = c(nind, nind, nmrk), dimnames = list(c(indnames), c(indnames), c(mrknames)))
  for(m in 1:nmrk) {
    G[,,m] <- t(Z[,m,])%*%Pi%*%Z[,m,]
  }

  nphe <- dim(pheno)[2]
  pheno.new <- as.matrix(pheno[which(rownames(pheno) %in% dimnames(G)[[1]]),])
  rownames(pheno.new) <- rownames(pheno)[which(rownames(pheno) %in% dimnames(G)[[1]])]

  dimnames(Z)[[1]] <- as.vector(t(outer(Pgametes, Qgametes, paste, sep=":"))) #EXCLUIR

  alleles <- matrix(unlist(strsplit(dimnames(Z)[[1]], '')), ncol=(ploidy+1), byrow=TRUE)[,-c((ploidy/2)+1)]
  X <- array(data = NA, dim = c(nind, (ploidy*2), nmrk), dimnames = list(c(indnames), letters[1:(ploidy*2)], c(mrknames)))
  for(m in 1:nmrk) {
    for(i in 1:nind) {
      a <- vector("list", (ploidy*2))
      names(a) <- letters[1:(ploidy*2)]
      for(j in 1:(ploidy*2)) {
        a[[j]] <- which(alleles == letters[j], arr.ind = TRUE)[,1]
        a[[j]] <- sum(Z[,m,i][Reduce(intersect, list(a[[j]]))])
      }
      X[i,,m] <- unlist(a)
    }
  }

  cat("Reading the following data: \n")
  cat("  Ploidy level:       ", ploidy, "\n", sep="")
  cat("  No. individuals:    ", nind, "\n", sep="")
  cat("  No. linkage groups: ", nlgs, "\n", sep="")
  cat("  Step size:          ", step, " cM \n", sep="")
  cat("  Map size:           ", round(last(cum.size), 2), " cM (", last(cum.nmrk), " positions) \n", sep="")
  cat("  No. phenotypes:     ", nphe, "\n", sep="")

  structure(list(ploidy = ploidy,
                 nlgs = nlgs,
                 nind = nind,
                 nmrk = nmrk,
                 nphe = nphe,
                 lgs.size = lgs.size,
                 cum.size = cum.size,
                 lgs.nmrk = lgs.nmrk,
                 cum.nmrk = cum.nmrk,
                 lgs = lgs,
                 lgs.all = lgs.all,
                 step = step,
                 pheno = pheno.new,
                 G = G,
                 Z = Z,
                 X = X,
                 Pi = Pi
                 ),
            class="qtlpoly.data")

}

#' @rdname read_data
#' @export
print.qtlpoly.data <- function(x, detailed = FALSE) {
  cat("This is an object of class 'qtlpoly.data'\n")
  cat("  Ploidy level:       ", x$ploidy, "\n", sep="")
  cat("  No. individuals:    ", x$nind, "\n", sep="")
  cat("  No. linkage groups: ", x$nlgs, "\n", sep="")
  cat("  Step size:          ", x$step, " cM \n", sep="")
  cat("  Map size:           ", round(last(x$cum.size), 2), " cM (", last(x$cum.nmrk), " positions) \n", sep="")
  if(detailed) for(c in 1:x$nlgs) cat("    LG ", c, ": ", round(x$lgs.size[[c]], 2), " cM (", x$lgs.nmrk[[c]], " positions) \n", sep="")
  cat("  No. phenotypes:     ", x$nphe, "\n", sep="")
  if(detailed) for(p in 1:x$nphe) cat("    Trait ", p, ": ", sQuote(colnames(x$pheno)[p]), " (", sum(!is.na(x$pheno[,p])), " individuals) \n", sep="")
}
