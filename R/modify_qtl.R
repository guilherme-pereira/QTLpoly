#' Modify QTL model
#'
#' Adds or removes QTL manually from a given model.
#'
#' @param model an object of class \code{qtlpoly.model} containing the QTL to be modified.
#'
#' @param pheno.col a phenotype column number whose model will be modified or printed.
#'
#' @param add.qtl a marker position number to be added.
#'
#' @param drop.qtl a marker position number to be removed.
#'
#' @param verbose if \code{TRUE} (default), current progress is shown; if \code{FALSE}, no output is produced.
#'
#' @param x an object of class \code{qtlpoly.modify} to be printed.
#'
#' @return An object of class \code{qtlpoly.modify} which contains a list of \code{results} for each trait with the following components:
#'
#'     \item{pheno.col}{a phenotype column number.}
#'     \item{stat}{a vector containing values from score statistics.}
#'     \item{pval}{a vector containing \emph{p}-values from score statistics.}
#'     \item{qtls}{a data frame with information from the mapped QTL.}
#'
#' @seealso \code{\link[qtlpoly]{read_data}}, \code{\link[qtlpoly]{remim}}
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
#'   # perform remim
#'   remim.mod <- remim(data = data, w.size = 15, sig.fwd = 0.01, sig.bwd = 0.0001,
#'     d.sint = 1.5, n.clusters = 4, plot = "remim")
#'
#'   # modify model
#'   modified.mod <- modify_qtl(model = remim.mod, pheno.col = 3, add.qtl = 184)
#'   }
#'
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2019) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{bioRxiv}. \url{doi:}.
#'
#' @export modify_qtl

modify_qtl <- function(model, pheno.col = NULL, add.qtl = NULL, drop.qtl = NULL, verbose = TRUE) {

  p <- which(model$pheno.col == pheno.col)
  stat <- model$results[[p]]$stat
  pval <- model$results[[p]]$pval

  if(!is.null(model$results[[p]]$qtls)) {
    qtl.mrk <- model$results[[p]]$qtl[,"Nmrk"]
    qtl.lgr <- model$results[[p]]$qtl[,"LG"]
    qtl.pos <- model$results[[p]]$qtl[,"Pos"]
  } else {
    qtl.mrk <- c()
    qtl.lgr <- c()
    qtl.pos <- c()
  }

  if(verbose) {
    if(length(qtl.mrk) == 0) cat("Model modification for trait ", pheno.col, " ", sQuote(names(model$results)[[p]]), "; there are no QTL in the model \n", sep = "")
    if(length(qtl.mrk) == 1) cat("Model modification for trait ", pheno.col, " ", sQuote(names(model$results)[[p]]), "; there is ", length(qtl.mrk), " QTL in the model already \n", sep = "")
    if(length(qtl.mrk) >= 2) cat("Model modification for trait ", pheno.col, " ", sQuote(names(model$results)[[p]]), "; there are ", length(qtl.mrk), " QTL in the model already \n", sep = "")
  }

  if(!is.null(add.qtl)) {
    qtl.mrk <- c(qtl.mrk, add.qtl)
    qtl.lgr <- c(qtl.lgr, last(which(last(qtl.mrk) > data$cum.nmrk)))
    qtl.pos <- c(qtl.pos, round(unlist(data$lgs)[[last(qtl.mrk)]], digits = 2))
    if(verbose) cat("  QTL was added on LG ", last(qtl.lgr), " at ", last(qtl.pos), " cM (position number ", last(qtl.mrk), ") \n\n", sep = "")
  }

  if(!is.null(drop.qtl)) {
    q <- which(qtl.mrk == drop.qtl)
    if(verbose) cat("  QTL was dropped from LG ", last(qtl.lgr[q]), " at ", last(qtl.pos[q]), " cM (position number ", last(qtl.mrk[q]), ") \n\n", sep = "")
    qtl.mrk <- qtl.mrk[-q]
    qtl.lgr <- qtl.lgr[-q]
    qtl.pos <- qtl.pos[-q]
  }

  if(length(qtl.mrk) > 0) {
    nqtl <- length(qtl.mrk)
    qtl <- c()
    for(q in 1:nqtl) {
      qtl <- c(qtl, c(qtl.lgr[q],
                      qtl.pos[q],
                      qtl.mrk[q],
                      names(unlist(data$lgs))[qtl.mrk[q]],
                      stat[qtl.mrk[q]],
                      pval[qtl.mrk[q]]))
    }
    qtls <- as.data.frame(matrix(qtl, ncol=6, byrow=TRUE), stringsAsFactors=FALSE)
    colnames(qtls) <- c("LG", "Pos", "Nmrk", "Mrk", "Score", "Pval")
    qtls[, c(1,2,3,5,6)] <- sapply(qtls[, c(1,2,3,5,6)], as.numeric)
    qtls[, c(2,5)] <- round(qtls[, c(2,5)], digits = 2)
    qtls[, c(6)] <- formatC(qtls[, c(6)], format="e", digits = 2)
    if(any(qtls[, c(6)] == "0.00e+00")) qtls[which(qtls[,6] == "0.00e+00"), c(6)] <- "<1.00e-16"
  } else {
    qtls <- NULL
  } # output QTL

  model$results[[p]] <- list(
    pheno.col=pheno.col,
    stat=stat,
    pval=pval,
    qtls=qtls)

  structure(list(data=model$data,
                 pheno.col=model$pheno.col,
                 w.size=model$w.size,
                 sig.fwd=model$sig.fwd,
                 sig.bwd=model$sig.fwd,
                 polygenes=model$polygenes,
                 d.sint=model$d.sint,
                 results=model$results),
            class=c("qtlpoly.modify"))

}

#' @rdname modify_qtl
#' @export
#'
print.qtlpoly.modify <- function(x, pheno.col=NULL) {
  if(any(class(x) == "qtlpoly.modify")) cat("This is an object of class 'qtlpoly.modify' and it has not been profiled yet\n")
  if(is.null(pheno.col)) {
    pheno.col <- 1:length(x$results)
  } else {
    pheno.col <- which(x$pheno.col %in% pheno.col)
  }
  for(p in pheno.col) {
    cat("\n* Trait", x$results[[p]]$pheno.col, sQuote(names(x$results)[[p]]), "\n")
    if(!is.null(x$results[[p]]$qtls)) print(x$results[[p]]$qtls)
    else cat("There are no QTL in the model \n")
  }
}
