#' Fixed-effect interval mapping (FEIM)
#'
#' Performs interval mapping using the single-QTL, fixed-effect model proposed by Hackett et al. (2001).
#'
#' @param data an object of class \code{qtlpoly.data}.
#'
#' @param pheno.col a numeric vector with the phenotype columns to be analyzed; if \code{NULL} (default), all phenotypes from \code{'data'} will be included.
#'
#' @param w.size a number representing the window size (in centiMorgans) to be avoided on either side of QTL already in the model when looking for a new QTL, e.g. 15 (default).
#'
#' @param sig.lod the vector of desired significance LOD thresholds (usually permutation-based) for declaring a QTL for each trait, e.g. 5 (default); if a single value is provided, the same LOD threshold will be applied to all traits.
#'
#' @param d.sint a \eqn{d} value to subtract from logarithm of the odds (\eqn{LOD-d}) for support interval calculation, e.g. \eqn{d=1.5} (default) represents approximate 95\% support interval.
#'
#' @param plot a suffix for the file's name containing plots of every algorithm step, e.g. "remim" (default); if \code{NULL}, no file is produced.
#'
#' @param verbose if \code{TRUE} (default), current progress is shown; if \code{FALSE}, no output is produced.
#'
#' @param x an object of class \code{qtlpoly.feim} to be printed.
#'
#' @param sint whether \code{"upper"} or \code{"lower"} support intervals should be printed; if \code{NULL} (default), QTL peak information will be printed.
#'
#' @return An object of class \code{qtlpoly.feim} which contains a list of \code{results} for each trait with the following components:
#'
#'     \item{pheno.col}{a phenotype column number.}
#'     \item{LRT}{a vector containing LRT values.}
#'     \item{LOD}{a vector containing LOD scores.}
#'     \item{AdjR2}{a vector containing adjusted \eqn{R^2}.}
#'     \item{qtls}{a data frame with information from the mapped QTL.}
#'     \item{lower}{a data frame with information from the lower support interval of mapped QTL.}
#'     \item{upper}{a data frame with information from the upper support interval of mapped QTL.}
#'
#' @seealso \code{\link[qtlpoly]{permutations}}
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
#'   feim.mod <- feim(data = data, sig.lod = 7, plot = "feim")
#'   }
#'
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2019) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{bioRxiv}. \url{doi:}.
#'
#'     Hackett, C.A., Bradshaw, J.E., McNicol, J.W. (2001) Interval mapping of quantitative trait loci in autotetraploid species, \emph{Genetics} 159: 1819-1832. \url{http://www.genetics.org/content/159/4/1819}
#'
#' @export feim

feim <- function(data = data, pheno.col = NULL, w.size = 15, sig.lod = 7, d.sint = 1.5, plot = "feim", verbose = TRUE) {

  w.size <- w.size/data$step
  start <- proc.time()
  if(is.null(pheno.col)) pheno.col <- 1:dim(data$pheno)[2]
  if(!is.null(plot)) pdf(paste(plot, "pdf", sep = "."))
  if(length(sig.lod) == 1) sig.lod <- rep(sig.lod, length(pheno.col))
  results <- vector("list", length(pheno.col))
  names(results) <- colnames(data$pheno)[pheno.col]

  for(p in 1:length(results)) {

    if(verbose) cat("FEIM for trait", pheno.col[p], sQuote(colnames(data$pheno)[pheno.col[p]]), "\n")
    # ind <- rownames(data$pheno)[which(!is.na(data$pheno[,pheno.col[p]]))]
    ind <- rownames(data$pheno)[which(dimnames(data$pheno)[[1]] %in% dimnames(data$X)[[1]])]
    Y <- data$pheno[ind,pheno.col[p]]

    LRT <- AdjR2 <- numeric(data$nmrk)
    for(m in 1:data$nmrk) {
      full.mod <- lm(Y ~ 1 + data$X[ind,-c(1,(data$ploidy+1)),m])
      null.mod <- lm(full.mod$model$Y ~ 1)
      LRT[m] <- -2 * (logLik(null.mod) - logLik(full.mod))
      AdjR2[m] <- summary(full.mod)$adj.r.squared
    }

    LOD <- LRT/(2*log(10))

    qtl.mrk <- c()
    qtl.lgr <- c()
    qtl.pos <- c()
    if(any(LOD > sig.lod[p])) {
      for(c in 1:data$nlgs) {
        x = which(LOD[(data$cum.nmrk[c]+1):data$cum.nmrk[c+1]] >= sig.lod[p])
        LODx = LOD[x+data$cum.nmrk[c]]
        if(length(x) > 0) {
          intervals <- c(0, which(diff(x) > w.size), length(x))
          for(i in 1:(length(intervals)-1)) {
            interval <- c((intervals[i]+1):intervals[i+1])
            qtl.mrk <- c(qtl.mrk, (x[interval[which.max(LODx[interval])]])+data$cum.nmrk[c])
            qtl.lgr <- c(qtl.lgr, last(which(last(qtl.mrk) > data$cum.nmrk)))
            qtl.pos <- c(qtl.pos, round(unlist(data$lgs)[[last(qtl.mrk)]], digits = 2))
            if(verbose) cat("  QTL was found on LG ", last(qtl.lgr), " at ", last(qtl.pos), " cM (position number ", last(qtl.mrk), ")\n", sep="")
          }
        }
      }
      if(verbose) cat("\n")
    } else {
      qtl.mrk0 <- which.max(LOD)
      qtl.lgr0 <- last(which(last(qtl.mrk0) > data$cum.nmrk))
      qtl.pos0 <- round(unlist(data$lgs)[[last(qtl.mrk0)]], digits = 2)
      if(verbose) cat("  No QTL were found. A putative QTL on LG ", last(qtl.lgr0), " at ", last(qtl.pos0), " cM (position number ", last(qtl.mrk0), ") did not reach the threshold; its LOD was ", round(max(LOD), 5), "\n\n", sep="")
    }

    if(length(qtl.mrk) > 0) {
      nqtl <- length(qtl.mrk)
      qtl <- c()
      for(q in 1:nqtl) {
        qtl <- c(qtl, c(qtl.lgr[q],
                        qtl.pos[q],
                        qtl.mrk[q],
                        names(unlist(data$lgs))[qtl.mrk[q]],
                        LRT[qtl.mrk[q]],
                        LOD[qtl.mrk[q]],
                        AdjR2[qtl.mrk[q]]))
      }
      qtls <- as.data.frame(matrix(qtl, ncol=7, byrow=TRUE), stringsAsFactors=FALSE)
      colnames(qtls) <- c("LG", "Pos", "Nmrk", "Mrk", "LRT", "LOD", "AdjR2")
      qtls[, c(1,2,3,5,6,7)] <- sapply(qtls[, c(1,2,3,5,6,7)], as.numeric)
      qtls[, c(2,5,6)] <- round(qtls[, c(2,5,6)], digits = 2)

      if(!is.null(d.sint)) {
        low <- upp <- c()
        for(q in 1:nqtl) {
          same.lgr <- which(!is.na(match(qtl.lgr, qtl.lgr[q])))
          if(length(same.lgr) > 1) {
            diff.mrk <- sort(qtl.mrk[same.lgr])
            midpoint <- diff.mrk[-length(diff.mrk)] + diff(diff.mrk)/2
            if(which(diff.mrk == qtl.mrk[q]) == 1) { # supports 3 QTL max in the same LG
              markers.out <- (data$cum.nmrk[qtl.lgr[q]]+1):floor(midpoint[1])#; print(markers.out)
            } else if(diff.mrk[which(diff.mrk == qtl.mrk[q])] == last(diff.mrk)) {
              markers.out <- (floor(last(midpoint))+1):(data$cum.nmrk[qtl.lgr[q]+1])#; print(markers.out)
            } else {
              markers.out <- (floor(midpoint[1])+1):(floor(midpoint[2]))#; print(markers.out)
            }
          } else {
            markers.out <- (data$cum.nmrk[qtl.lgr[q]]+1):(data$cum.nmrk[qtl.lgr[q]+1])#; print(markers.out)
          }
          lod.out <- LOD[markers.out]
          lod.qtl <- LOD[qtl.mrk[q]]
          lower <- head(markers.out[which(lod.out >= (lod.qtl - d.sint))],1)
          upper <- tail(markers.out[which(lod.out >= (lod.qtl - d.sint))],1)
          low <- c(low, c(qtl.lgr[q],
                          round(unlist(data$lgs)[[lower]], digits = 2),
                          lower,
                          names(unlist(data$lgs))[lower],
                          LRT[lower],
                          LOD[lower],
                          AdjR2[lower]))
          upp <- c(upp, c(qtl.lgr[q],
                          round(unlist(data$lgs)[[upper]], digits = 2),
                          upper,
                          names(unlist(data$lgs))[upper],
                          LRT[upper],
                          LOD[upper],
                          AdjR2[upper]))
        }
        lower <- as.data.frame(matrix(low, ncol=7, byrow=TRUE), stringsAsFactors=FALSE)
        colnames(lower) <- c("LG", "Pos_lower", "Nmrk_lower", "Mrk_lower", "LRT_lower", "LOD_lower", "AdjR2_lower")
        lower[, c(1,2,3,5,6,7)] <- sapply(lower[, c(1,2,3,5,6,7)], as.numeric)
        lower[, c(2,5,6)] <- round(lower[, c(2,5)], digits = 2)

        upper <- as.data.frame(matrix(upp, ncol=7, byrow=TRUE), stringsAsFactors=FALSE)
        colnames(upper) <- c("LG", "Pos_upper", "Nmrk_upper", "Mrk_upper", "LRT_upper", "LOD_upper", "AdjR2_upper")
        upper[, c(1,2,3,5,6,7)] <- sapply(upper[, c(1,2,3,5,6,7)], as.numeric)
        upper[, c(2,5,6)] <- round(upper[, c(2,5)], digits = 2)

      } else {
        lower <- upper <- NULL
      }

    } else {
      qtls <- lower <- upper <- NULL
    }

    if(!is.null(plot)) {
      plot(LOD, xlab="Position number", main=paste("Trait ", pheno.col[p], " '", colnames(data$pheno)[pheno.col[p]], "'", sep=""), ylim=c(0,max(LOD)))
      abline(v=data$cum.nmrk, lty=3)
      if(!is.null(sig.lod[p])) abline(h=sig.lod[p], lty=5)
      if(!is.null(qtl.mrk)) points(x=qtl.mrk, y=rep(0, length(qtl.mrk)), pch=6, lwd=1.5, col="red")
    }

    results[[p]] <- list(
      pheno.col=pheno.col[p],
      LRT=LRT,
      LOD=LOD,
      AdjR2=AdjR2,
      qtls=qtls,
      lower=lower,
      upper=upper)

    end <- proc.time()
  }

  if(verbose) cat("Calculation took", round((end - start)[3], digits = 2), "seconds\n\n")
  if(!is.null(plot)) dev.off()

  structure(list(data=deparse(substitute(data)),
                 pheno.col=pheno.col,
                 w.size=w.size,
                 sig.lod=sig.lod,
                 results=results),
            class=c("qtlpoly.feim"))

}

#' @rdname feim
#' @export
print.qtlpoly.feim <- function(x, pheno.col = NULL, sint=NULL) {
  if(any(class(x) == "qtlpoly.feim")) cat("This is an object of class 'qtlpoly.feim'\n")
  if(is.null(pheno.col)) {
    pheno.col <- 1:length(x$results)
  } else {
    pheno.col <- which(x$pheno.col %in% pheno.col)
  }
  for(p in pheno.col) {
    cat("\n* Trait", x$results[[p]]$pheno.col, sQuote(names(x$results)[[p]]), "\n")
    if(is.null(sint)) {
      if(!is.null(x$results[[p]]$qtls)) print(x$results[[p]]$qtls)
      else cat("There are no QTL in the model \n")
    } else if(sint=="lower") {
      if(!is.null(x$results[[p]]$lower)) print(x$results[[p]]$lower)
      else cat("There are no QTL in the model \n")
    } else if(sint=="upper") {
      if(!is.null(x$results[[p]]$upper)) print(x$results[[p]]$upper)
      else cat("There are no QTL in the model \n")
    }
  }
}
