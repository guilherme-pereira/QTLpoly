#' QTL profiling
#'
#' Generates the score-based genome-wide profile conditional to the selected QTL.
#'
#' @param data an object of class \code{qtlpoly.data}.
#'
#' @param model an object of class \code{qtlpoly.model} containing the QTL to be profiled.
#'
#' @param d.sint a \eqn{d} value to subtract from logarithm of \emph{p}-value (\eqn{LOP-d}) for support interval calculation, e.g. \eqn{d=1.5} (default) represents approximate 95\% support interval.
#'
#' @param polygenes if \code{TRUE} all QTL but the one being tested are treated as a single polygenic effect, if \code{FALSE} (default) all QTL effect variances have to estimated.
#'
#' @param n.clusters number of parallel processes to spawn.
#'
#' @param plot a suffix for the file's name containing plots of every QTL profiling round, e.g. "profile" (default); if \code{NULL}, no file is produced.
#'
#' @param verbose if \code{TRUE} (default), current progress is shown; if \code{FALSE}, no output is produced.
#'
#' @param x an object of class \code{qtlpoly.profile} to be printed.
#'
#' @param sint whether \code{"upper"} or \code{"lower"} support intervals should be printed; if \code{NULL} (default), only QTL peak information will be printed.
#'
#' @param pheno.col a numeric vector with the phenotype column numbers to be plotted; if \code{NULL}, all phenotypes from \code{'data'} will be included.
#'
#' @return An object of class \code{qtlpoly.profile} which contains a list of \code{results} for each trait with the following components:
#'
#'     \item{pheno.col}{a phenotype column number.}
#'     \item{stat}{a vector containing values from score statistics.}
#'     \item{pval}{a vector containing \emph{p}-values from score statistics.}
#'     \item{qtls}{a data frame with information from the mapped QTL.}
#'     \item{lower}{a data frame with information from the lower support interval of mapped QTL.}
#'     \item{upper}{a data frame with information from the upper support interval of mapped QTL.}
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
#'
#'   # profile model
#'   profile.mod <- profile_qtl(data = data, model = optimize.mod, d.sint = 1.5,
#'     polygenes = FALSE, n.clusters = 4, plot = "profile")
#'   }
#'
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2019) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{bioRxiv}. \url{doi:}.
#'     
#'     Qu L, Guennel T, Marshall SL (2013) Linear score tests for variance components in linear mixed models and applications to genetic association studies. \emph{Biometrics} 69 (4): 883–92. \url{doi:10.1111/biom.12095}.
#'
#' @export profile_qtl
#' @import varComp parallel

profile_qtl <- function(data, model, d.sint = 1.5, polygenes = FALSE, n.clusters = NULL, plot = "profile", verbose = TRUE) {

  if(is.null(n.clusters)) n.clusters <- 1
  cat("INFO: Using", n.clusters, "CPUs for calculation\n\n")
  cl <- makeCluster(n.clusters)
  clusterEvalQ(cl, require(varComp))

  if(!is.null(plot)) plot <- paste(plot, "pdf", sep = ".")
  sig.bwd <- model$sig.bwd
  results <- vector("list", length(model$results))
  names(results) <- names(model$results)
  w.size <- model$w.size/data$step

  for(p in 1:length(results)) {

    start <- proc.time()
    pheno.col <- model$results[[p]]$pheno.col
    qtl.mrk <- model$results[[p]]$qtl[,"Nmrk"]
    qtl.lgr <- model$results[[p]]$qtl[,"LG"]
    qtl.pos <- model$results[[p]]$qtl[,"Pos"]
    stat <- model$results[[p]]$stat
    pval <- model$results[[p]]$pval
    lower <- upper <- qtl.mrk
    if(verbose) {
      if(length(qtl.mrk) == 0) cat("QTL profile for trait ", pheno.col, " ", sQuote(colnames(data$pheno)[pheno.col]), "; there are no QTL in the model \n", sep="")
      if(length(qtl.mrk) == 1) cat("QTL profile for trait ", pheno.col, " ", sQuote(colnames(data$pheno)[pheno.col]), "; there is ", length(qtl.mrk), " QTL in the model \n", sep="")
      if(length(qtl.mrk) >= 2) cat("QTL profile for trait ", pheno.col, " ", sQuote(colnames(data$pheno)[pheno.col]), "; there are ", length(qtl.mrk), " QTL in the model \n", sep="")
    }
    if(!is.null(plot)) pdf(paste(colnames(data$pheno)[pheno.col], plot, sep = "_"))
    ind <- rownames(data$pheno)[which(!is.na(data$pheno[,pheno.col]))]
    Y <- data$pheno[ind,pheno.col]

    #begin profile
    if(length(qtl.mrk) == 0) {
      markers <- c(1:data$nmrk)
      temp <- parSapply(cl, as.character(markers), function(x) { #like first search
        m <- as.numeric(x)
        full.mod <- varComp(Y ~ 1, varcov = list(data$G[ind,ind,m]))
        test <- varComp.test(full.mod, null=integer(0L))
        c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
      })
      stat[as.numeric(colnames(temp))] <- temp["st",]
      pval[as.numeric(colnames(temp))] <- temp["pv",]
      if(!is.null(plot)) {
        plot(-log10(pval), xlab="Marker number", ylab="-log10(p)", main="Completing genome", ylim=c(0,10))
        abline(v=data$cum.nmrk, lty=3); abline(h=-log10(sig.bwd), lty=5)
      }
    } # completing genome when there's no QTL

    if(length(qtl.mrk) == 1) {
      if(verbose) cat("  Profiling QTL ...", qtl.mrk, "\n")
      markers.out <- c((data$cum.nmrk[qtl.lgr[1]]+1):(data$cum.nmrk[qtl.lgr[1]+1]))
      markers <- c(1:data$nmrk)[-markers.out]
      temp <- parSapply(cl, as.character(markers.out), function(x) { #like first search
        m <- as.numeric(x)
        full.mod <- varComp(Y ~ 1, varcov = list(data$G[ind,ind,m]))
        test <- varComp.test(full.mod, null=integer(0L))
        c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
      })
      stat[as.numeric(colnames(temp))] <- temp["st",]
      pval[as.numeric(colnames(temp))] <- temp["pv",]
      qtl.vcv <- list(data$G[ind,ind,qtl.mrk[1]])
      withCallingHandlers(full.mod0 <- varComp(Y ~ 1, varcov = c(qtl.vcv)), warning = h)
      control <- varComp.control(start = c(coef(full.mod0, what = "var.ratio"),0))
      temp <- parSapply(cl, as.character(markers), function(x) {
        m <- as.numeric(x)
        full.mod <- varComp(Y ~ 1, varcov = c(qtl.vcv, list(data$G[ind,ind,m])), control = control)
        test <- varComp.test(full.mod, null=1L)
        c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
      })
      stat[as.numeric(colnames(temp))] <- temp["st",]
      pval[as.numeric(colnames(temp))] <- temp["pv",]
      if(!is.null(plot)) {
        plot(-log10(pval), xlab="Marker number", ylab="-log10(p)", main="Profiling round #1", ylim=c(0,10))
        abline(v=data$cum.nmrk, lty=3); points(x=qtl.mrk, y=rep(-0.15, length(qtl.mrk)), pch=6, lwd=1.5, col="red"); abline(h=-log10(sig.bwd), lty=5)
      }
      if(!is.null(d.sint)) {
        lop.out <- -log10(pval[markers.out])
        lop.qtl <- -log10(pval[qtl.mrk[1]])
        lower[1] <- head(markers.out[which(lop.out >= (lop.qtl - d.sint))],1)
        upper[1] <- tail(markers.out[which(lop.out >= (lop.qtl - d.sint))],1)
      } # lower and upper markers for support interval
    } # profiling 1 QTL

    if(length(qtl.mrk) > 1) {
      qtl.lgr <- qtl.lgr[order(qtl.mrk)]
      qtl.mrk <- sort(qtl.mrk)
      markers <- c()
      if(verbose) cat("  Profiling QTL ")
      for(q in 1:length(qtl.mrk)) {
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
        markers <- c(markers, markers.out)
        qtl.vcv <- NULL
        qtl.mrk0 <- c()
        for(q0 in which(qtl.mrk != qtl.mrk[q])) {
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

        qtl.mrk[q] <- markers.out[which.max(stat[markers.out])]
        qtl.pos[q] <- round(unlist(data$lgs)[qtl.mrk[q]], digits = 2)

        if(!is.null(d.sint)) {
          lop.out <- -log10(pval[markers.out])
          lop.qtl <- -log10(pval[qtl.mrk[q]])
          lower[q] <- head(markers.out[which(lop.out >= (lop.qtl - d.sint))],1)
          upper[q] <- tail(markers.out[which(lop.out >= (lop.qtl - d.sint))],1)
        } # lower and upper markers for support interval

        if(verbose) {
          cat("...", qtl.mrk[q], "")
          if(q == length(qtl.mrk)) cat("\n")
        }
        if(!is.null(plot)) {
          plot(-log10(pval), xlab="Marker number", ylab="-log10(p)", main=paste("Profiling round #", q, sep=""), ylim=c(0,10))
          abline(v=data$cum.nmrk, lty=3); points(x=qtl.mrk, y=rep(-0.15, length(qtl.mrk)), pch=6, lwd=1.5, col="red"); abline(h=-log10(sig.bwd), lty=5)
        }
      }
      qtl.vcv <- NULL
      markers.out <- c(1:data$nmrk)[-markers]
      if(length(markers.out) > 0) {
        for(q in 1:length(qtl.mrk)) { # completando todo genoma com os qtls finais
          qtl.vcv <- c(qtl.vcv, list(data$G[ind,ind,qtl.mrk[q]]))
        }
        if(polygenes) {
          Gstar <- apply(data$G[ind,ind,qtl.mrk], MARGIN = c(1,2), sum)/length(qtl.mrk); Gstar[1:5,1:5]
          full.mod0 <- varComp(Y ~ 1, varcov = list(Gstar))
          control <- varComp.control(start = c(coef(full.mod0, what = "var.ratio"),0))
          temp <- parSapply(cl, as.character(markers.out), function(x) {
            m <- as.numeric(x)
            full.mod <- varComp(Y ~ 1, varcov = list(Gstar, data$G[ind,ind,m]), control = control)
            test <- varComp.test(full.mod, null=1L)
            c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
          })
        } else {
          withCallingHandlers(full.mod <- varComp(Y ~ 1, varcov = c(qtl.vcv)), warning = h)
          control <- varComp.control(start = c(coef(full.mod, what = "var.ratio"),0))
          temp <- parSapply(cl, as.character(markers.out), function(x) {
            m <- as.numeric(x)
            full.mod <- varComp(Y ~ 1, varcov = c(qtl.vcv, list(data$G[ind,ind,m])), control = control)
            test <- varComp.test(full.mod, null=c(1:length(qtl.vcv)))
            c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
          })
        }
        stat[as.numeric(colnames(temp))] <- temp["st",]
        pval[as.numeric(colnames(temp))] <- temp["pv",]
        if(!is.null(plot)) {
          plot(-log10(pval), xlab="Marker number", ylab="-log10(p)", main="Completing genome", ylim=c(0,10))
          abline(v=data$cum.nmrk, lty=3); points(x=qtl.mrk, y=rep(-0.15, length(qtl.mrk)), pch=6, lwd=1.5, col="red"); abline(h=-log10(sig.bwd), lty=5)
        }
      }
    } # profiling 1+ QTL
    #end profile
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

      if(!is.null(d.sint)) {
        low <- upp <- c()
        for(q in 1:nqtl) {
          low <- c(low, c(qtl.lgr[q],
                          round(unlist(data$lgs)[[lower[q]]], digits = 2),
                          lower[q],
                          names(unlist(data$lgs))[lower[q]],
                          stat[lower[q]],
                          pval[lower[q]]))
          upp <- c(upp, c(qtl.lgr[q],
                          round(unlist(data$lgs)[[upper[q]]], digits = 2),
                          upper[q],
                          names(unlist(data$lgs))[upper[q]],
                          stat[upper[q]],
                          pval[upper[q]]))
        }
        lower <- as.data.frame(matrix(low, ncol=6, byrow=TRUE), stringsAsFactors=FALSE)
        colnames(lower) <- c("LG", "Pos_lower", "Nmrk_lower", "Mrk_lower", "Score_lower", "Pval_lower")
        lower[, c(1,2,3,5,6)] <- sapply(lower[, c(1,2,3,5,6)], as.numeric)
        lower[, c(2,5)] <- round(lower[, c(2,5)], digits = 2)
        lower[, c(6)] <- formatC(lower[, c(6)], format="e", digits = 2)
        if(any(lower[, c(6)] == "0.00e+00")) lower[which(lower[,6] == "0.00e+00"), c(6)] <- "<2.22e-16"

        upper <- as.data.frame(matrix(upp, ncol=6, byrow=TRUE), stringsAsFactors=FALSE)
        colnames(upper) <- c("LG", "Pos_upper", "Nmrk_upper", "Mrk_upper", "Score_upper", "Pval_upper")
        upper[, c(1,2,3,5,6)] <- sapply(upper[, c(1,2,3,5,6)], as.numeric)
        upper[, c(2,5)] <- round(upper[, c(2,5)], digits = 2)
        upper[, c(6)] <- formatC(upper[, c(6)], format="e", digits = 2)
        if(any(upper[, c(6)] == "0.00e+00")) upper[which(upper[,6] == "0.00e+00"), c(6)] <- "<2.22e-16"

      } else {
        lower <- upper <- NULL
      }

    } else {
      qtls <- lower <- upper <- NULL
    } # output QTL plus support interval lower and upper bounds

    results[[p]] <- list(
      pheno.col=pheno.col,
      stat=stat,
      pval=pval,
      qtls=qtls,
      lower=lower,
      upper=upper)

  }

  stopCluster(cl)

  structure(list(data=deparse(substitute(data)),
                 pheno.col=model$pheno.col,
                 w.size=model$w.size,
                 sig.fwd=model$sig.fwd,
                 sig.bwd=model$sig.bwd,
                 polygenes=polygenes,
                 d.sint=d.sint,
                 results=results),
            class=c("qtlpoly.model","qtlpoly.profile"))

}

#' @rdname profile_qtl
#' @export

print.qtlpoly.profile <- function(x, pheno.col = NULL, sint=NULL) {
  if(any(class(x) == "qtlpoly.profile")) cat("This is an object of class 'qtlpoly.profile'\n")
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
