#' Model optimization
#'
#' Tests each QTL at a time and updates its position (if it changes) or drops the QTL (if non-significant).
#'
#' @param data an object of class \code{qtlpoly.data}.
#'
#' @param model an object of class \code{qtlpoly.model} containing the QTL to be optimized.
#'
#' @param sig.bwd the desired score-based \emph{p}-value threshold for backward elimination, e.g. 0.0001 (default).
#'
#' @param polygenes if \code{TRUE} all QTL but the one being tested are treated as a single polygenic effect, if \code{FALSE} (default) all QTL effect variances have to estimated.
#'
#' @param n.clusters number of parallel processes to spawn.
#'
#' @param plot a suffix for the file's name containing plots of every QTL optimization round, e.g. "optimize" (default); if \code{NULL}, no file is produced.
#'
#' @param verbose if \code{TRUE} (default), current progress is shown; if \code{FALSE}, no output is produced.
#'
#' @param x an object of class \code{qtlpoly.optimize} to be printed.
#'
#' @param pheno.col a numeric vector with the phenotype columns to be printed; if \code{NULL}, all phenotypes from \code{'data'} will be included.
#'
#' @return An object of class \code{qtlpoly.optimize} which contains a list of \code{results} for each trait with the following components:
#'
#'     \item{pheno.col}{a phenotype column number.}
#'     \item{stat}{a vector containing values from score statistics.}
#'     \item{pval}{a vector containing \emph{p}-values from score statistics.}
#'     \item{qtls}{a data frame with information from the mapped QTL.}
#'
#' @seealso \code{\link[qtlpoly]{read_data}}, \code{\link[qtlpoly]{null_model}}, \code{\link[qtlpoly]{search_qtl}}
#'
#' @examples
#'   \dontrun{
#'   # load raw data
#'   data(maps)
#'   data(pheno)
#'
#'   # estimate conditional probabilities using 'mappoly' package
#'   library(mappoly)
#'   genoprob <- lapply(maps, calc_genoprob)
#'
#'   # prepare data
#'   data <- read_data(ploidy = 6, geno.prob = genoprob, pheno = pheno, step = 1)
#'
#'   # build null models
#'   null.mod <- null_model(data = data, n.clusters = 4, plot = "null")
#'
#'   # perform forward search
#'   search.mod <- search_qtl(data = data, model = null.mod, w.size = 15, sig.fwd = 0.01,
#'     n.clusters = 4, plot = "search")
#'
#'   # optimize model
#'   optimize.mod <- optimize_qtl(data = data, model = search.mod, sig.bwd = 0.0001,
#'     n.clusters = 4, plot = "optimize")
#'   }
#'
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2019) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{bioRxiv}. \url{doi:}.
#'     
#'     Qu L, Guennel T, Marshall SL (2013) Linear score tests for variance components in linear mixed models and applications to genetic association studies. \emph{Biometrics} 69 (4): 883–92. \url{doi:10.1111/biom.12095}.
#'
#' @export optimize_qtl
#' @import varComp parallel

optimize_qtl <- function(data, model, sig.bwd = 0.0001, polygenes = FALSE, n.clusters = NULL, plot = "optimize", verbose = TRUE) {

  if(is.null(n.clusters)) n.clusters <- 1
  cat("INFO: Using", n.clusters, "CPUs for calculation\n\n")
  cl <- makeCluster(n.clusters)
  clusterEvalQ(cl, require(varComp))

  if(!is.null(plot)) plot <- paste(plot, "pdf", sep = ".")
  results <- vector("list", length(model$results))
  names(results) <- names(model$results)
  w.size <- model$w.size/data$step

  for(p in 1:length(results)) {

    start <- proc.time()
    pheno.col <- model$results[[p]]$pheno.col
    w.size <- w.size/data$step
    stat <- model$results[[p]]$stat
    pval <- model$results[[p]]$pval
    qtl.mrk <- model$results[[p]]$qtl[,"Nmrk"]
    qtl.lgr <- model$results[[p]]$qtl[,"LG"]
    qtl.pos <- model$results[[p]]$qtl[,"Pos"]
    if(verbose) {
      if(length(qtl.mrk) == 0) cat("Model optimization for trait ", pheno.col, " ", sQuote(colnames(data$pheno)[pheno.col]), "; there are no QTL in the model \n", sep="")
      if(length(qtl.mrk) == 1) cat("Model optimization for trait ", pheno.col, " ", sQuote(colnames(data$pheno)[pheno.col]), "; there is ", length(qtl.mrk), " QTL in the model already \n", sep="")
      if(length(qtl.mrk) >= 2) cat("Model optimization for trait ", pheno.col, " ", sQuote(colnames(data$pheno)[pheno.col]), "; there are ", length(qtl.mrk), " QTL in the model already \n", sep="")
    }
    if(!is.null(plot)) pdf(paste(colnames(data$pheno)[pheno.col], plot, sep = "_"))
    ind <- rownames(data$pheno)[which(!is.na(data$pheno[,pheno.col]))]
    Y <- data$pheno[ind,pheno.col]

    qtl.out <- c(1)
    while(!is.null(qtl.out)) {
      qtl.out <- c()
      if(length(qtl.mrk) == 1) {
        if(verbose) cat("  Refining QTL positions ...", qtl.mrk, "\n")
        markers.out <- c((data$cum.nmrk[qtl.lgr[1]]+1):(data$cum.nmrk[qtl.lgr[1]+1]))
        temp <- parSapply(cl, as.character(markers.out), function(x) { #like first search
          m <- as.numeric(x)
          full.mod <- varComp(Y ~ 1, varcov = list(data$G[ind,ind,m]))
          test <- varComp.test(full.mod, null=integer(0L))
          c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
        })
        stat[as.numeric(colnames(temp))] <- temp["st",]
        pval[as.numeric(colnames(temp))] <- temp["pv",]
        if(pval[markers.out[which.max(stat[markers.out])]] <= sig.bwd) { # updates position
          qtl.mrk[1] <- markers.out[which.max(stat[markers.out])]
        } else { # stores non-significant
          qtl.out <- c(1)
        }
        if(!is.null(plot)) {
          plot(-log10(pval), xlab="Marker number", ylab="-log10(p)", main="Refining round #1", ylim=c(0,10))
          abline(v=data$cum.nmrk, lty=3); abline(h=-log10(sig.bwd), lty=5); points(x=qtl.mrk, y=rep(-0.15, length(qtl.mrk)), pch=6, lwd=1.5, col="red")
          if(!is.null(qtl.out)) points(x=qtl.mrk[qtl.out], y=rep(-0.15, length(qtl.out)), pch=4, lwd=1.5, col="red")
        }
        if(!is.null(qtl.out)) {
          if(verbose) cat("  Excluding non-significant QTL", paste("...", qtl.mrk[qtl.out]), "\n")
          qtl.mrk <- qtl.mrk[-qtl.out]
          qtl.lgr <- qtl.lgr[-qtl.out]
        }
      }
      if(length(qtl.mrk) > 1) {
        if(verbose) cat("  Refining QTL positions ")
        for(q in 1:length(qtl.mrk)) {
          if(length(qtl.mrk) > 1 & (length(qtl.mrk)-length(qtl.out)) > 1) {
            same.lgr <- which(!is.na(match(qtl.lgr, qtl.lgr[q])))
            if(length(same.lgr) > 1) {
              diff.mrk <- sort(qtl.mrk[same.lgr])
              midpoint <- diff.mrk[-length(diff.mrk)] + diff(diff.mrk)/2
              if(which(diff.mrk == qtl.mrk[q])[1] == 1) { # supports up to 3 QTL in the same LG
                markers.out <- (data$cum.nmrk[qtl.lgr[q]]+1):floor(midpoint[1])
              } else if(diff.mrk[which(diff.mrk == qtl.mrk[q])] == last(diff.mrk)) {
                markers.out <- (floor(last(midpoint))+1):(data$cum.nmrk[qtl.lgr[q]+1])
              } else {
                markers.out <- (floor(midpoint[1])+1):(floor(midpoint[2]))
              }
            } else {
              markers.out <- (data$cum.nmrk[qtl.lgr[q]]+1):(data$cum.nmrk[qtl.lgr[q]+1])
            }
            qtl.vcv <- NULL
            qtl.mrk0 <- c()
            for(q0 in which(!(qtl.mrk %in% c(qtl.mrk[q], qtl.mrk[qtl.out])))) {
              qtl.vcv <- c(qtl.vcv, list(data$G[ind,ind,qtl.mrk[q0]]))
              qtl.mrk0 <- c(qtl.mrk0, qtl.mrk[q0])
            }
            if(polygenes) {
              Gstar <- apply(data$G[ind,ind,qtl.mrk0], MARGIN = c(1,2), sum)/length(qtl.mrk0); Gstar[1:5,1:5]
              full.mod0 <- varComp(Y ~ 1, varcov = list(Gstar))
              control <- varComp.control(start = c(coef(full.mod0, what = "var.ratio"),0))
              temp <- parSapply(cl, as.character(markers.out), function(x) {
                m <- as.numeric(x)
                full.mod <- varComp(Y ~ 1, varcov = list(Gstar, data$G[ind,ind,m]), control = control)
                test <- varComp.test(full.mod, null=1L)
                c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
              })
            } else {
              withCallingHandlers(full.mod0 <- varComp(Y ~ 1, varcov = c(qtl.vcv)), warning = h)
              control <- varComp.control(start = c(coef(full.mod0, what = "var.ratio"),0))
              temp <- parSapply(cl, as.character(markers.out), function(x) {
                m <- as.numeric(x)
                full.mod <- varComp(Y ~ 1, varcov = c(qtl.vcv, list(data$G[ind,ind,m])), control = control)
                test <- varComp.test(full.mod, null=c(1:length(qtl.vcv)))
                c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
              })
            }
            stat[as.numeric(colnames(temp))] <- temp["st",]
            pval[as.numeric(colnames(temp))] <- temp["pv",]
            if(pval[markers.out[which.max(stat[markers.out])]] <= sig.bwd) { # updates position
              qtl.mrk[q] <- markers.out[which.max(stat[markers.out])]
            } else { # stores non-significant
              qtl.out <- c(qtl.out, q)
            }
            if(!is.null(plot)) {
              plot(-log10(pval), xlab="Marker number", ylab="-log10(p)", main=paste("Refining round #", q, sep=""), ylim=c(0,10))
              abline(v=data$cum.nmrk, lty=3); abline(h=-log10(sig.bwd), lty=5); points(x=qtl.mrk, y=rep(-0.15, length(qtl.mrk)), pch=6, lwd=1.5, col="red")
              if(!is.null(qtl.out)) points(x=qtl.mrk[qtl.out], y=rep(-0.15, length(qtl.out)), pch=4, lwd=1.5, col="red")
            }
          }
          if(verbose) {
            if(length(qtl.mrk) > 1) cat("...", qtl.mrk[q], "")
            if(length(qtl.mrk) > 1 & q == length(qtl.mrk)) cat("\n")
            if(q > length(qtl.mrk) & !is.null(qtl.out)) cat(paste("...", qtl.mrk), "\n")
          }
        }
        if(!is.null(qtl.out) & q == (length(qtl.out)+1)) qtl.out <- unique(c(qtl.out, q))
        if(!is.null(qtl.out) & q <= length(qtl.mrk)) {
          if(verbose) cat("  Excluding non-significant QTL", paste("...", qtl.mrk[qtl.out]), "\n")
          # if(!is.null(plot)) {
          #   plot(-log10(pval), xlab="Marker number", ylab="-log10(p)", main=paste("Refining round #", q, sep=""), ylim=c(0,10))
          #   abline(v=data$cum.nmrk, lty=3); abline(h=-log10(sig.bwd), lty=5); points(x=qtl.mrk, y=rep(-0.15, length(qtl.mrk)), pch=6, lwd=1.5, col="red")
          #   if(!is.null(qtl.out)) points(x=qtl.mrk[qtl.out], y=rep(-0.15, length(qtl.out)), pch=4, lwd=1.5, col="red")
          # }
          qtl.mrk <- qtl.mrk[-qtl.out]
          qtl.lgr <- qtl.lgr[-qtl.out]
        }
      }
    } # keeps refining until all QTL are significant
    #end backward

    if(!is.null(plot)) dev.off()
    end <- proc.time()
    if(verbose) cat("  Calculation took", round((end - start)[3], digits = 2), "seconds\n\n")

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
      if(any(qtls[, c(6)] == "0.00e+00")) qtls[which(qtls[,6] == "0.00e+00"), c(6)] <- "<2.22e-16"
    } else {
      qtls <- NULL
    } # output QTL

    results[[p]] <- list(
      pheno.col=pheno.col,
      stat=stat,
      pval=pval,
      qtls=qtls)

  }

  stopCluster(cl)

  structure(list(data=deparse(substitute(data)),
                 pheno.col=model$pheno.col,
                 w.size=model$w.size,
                 sig.fwd=model$sig.fwd,
                 sig.bwd=sig.bwd,
                 polygenes=polygenes,
                 d.sint=NULL,
                 results=results),
            class=c("qtlpoly.model","qtlpoly.optimize"))

}

#' @rdname optimize_qtl
#' @export

print.qtlpoly.optimize <- function(x, pheno.col = NULL) {
  if(any(class(x) == "qtlpoly.optimize")) cat("This is an object of class 'qtlpoly.optimize'\n")
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
