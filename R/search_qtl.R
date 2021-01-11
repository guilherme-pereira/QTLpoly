#' QTL forward search
#'
#' Searches for QTL and adds them one at a time to a multiple random-effect QTL model based on score statistics.
#'
#' @param data an object of class \code{qtlpoly.data}.
#'
#' @param offset.data a data frame with the same dimensions of \code{data$pheno} containing offset variables; if \code{NULL} (default), no offset variables are considered.
#' 
#' @param model an object of class \code{qtlpoly.model} from which a forward search will start.
#'
#' @param w.size the window size (in cM) to avoid on either side of QTL already in the model when looking for a new QTL.
#'
#' @param sig.fwd the desired score-based \emph{p}-value threshold for forward search, e.g. 0.01 (default).
#'
#' @param score.null an object of class \code{qtlpoly.null} with results of score statistics from resampling.
#'
#' @param polygenes if \code{TRUE} all QTL but the one being tested are treated as a single polygenic effect; if \code{FALSE} (default) all QTL effect variances have to estimated.
#'
#' @param n.rounds number of search rounds; if \code{Inf} (default) forward search will stop when no more significant positions can be found.
#'
#' @param n.clusters number of parallel processes to spawn.
#'
#' @param plot a suffix for the file's name containing plots of every QTL search round, e.g. "search" (default); if \code{NULL}, no file is produced.
#'
#' @param verbose if \code{TRUE} (default), current progress is shown; if \code{FALSE}, no output is produced.
#'
#' @param x an object of class \code{qtlpoly.search} to be printed.
#'
#' @param pheno.col a numeric vector with the phenotype column numbers to be printed; if \code{NULL}, all phenotypes from \code{'data'} will be included.
#'
#' @return An object of class \code{qtlpoly.search} which contains a list of \code{results} for each trait with the following components:
#'
#'     \item{pheno.col}{a phenotype column number.}
#'     \item{stat}{a vector containing values from score statistics.}
#'     \item{pval}{a vector containing \emph{p}-values from score statistics.}
#'     \item{qtls}{a data frame with information from the mapped QTL.}
#'
#' @seealso \code{\link[qtlpoly]{read_data}}, \code{\link[qtlpoly]{null_model}}
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
#'   search.mod <- search(data = data, model = null.mod, w.size = 15, sig.fwd = 0.01,
#'     n.clusters = 4, plot = "search")
#'   }
#'
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2020) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{Genetics} 215 (3): 579-595. \url{http://doi.org/10.1534/genetics.120.303080}.
#'     
#'     Qu L, Guennel T, Marshall SL (2013) Linear score tests for variance components in linear mixed models and applications to genetic association studies. \emph{Biometrics} 69 (4): 883–92. \url{doi.org/10.1111/biom.12095}.
#'
#'     Zou F, Fine JP, Hu J, Lin DY (2004) An efficient resampling method for assessing genome-wide statistical significance in mapping quantitative trait loci. \emph{Genetics} 168 (4): 2307-16. \url{doi.org/10.1534/genetics.104.031427}
#'
#' @export search_qtl

search_qtl <- function(data, offset.data = NULL, model, w.size = 15, sig.fwd = 0.20, score.null = NULL, polygenes = FALSE, n.rounds = Inf, n.clusters = NULL, plot = "search", verbose = TRUE) {
  
  if(is.null(n.clusters)) n.clusters <- 1
  cat("INFO: Using", n.clusters, "CPUs for calculation\n\n")
  cl <- makeCluster(n.clusters)
  clusterEvalQ(cl, require(varComp))
  
  sig.fwd0 <- sig.fwd
  
  min.pvl <- NULL
  if(!is.null(score.null)) {
    min.pvl <- numeric(length(score.null$results))
    for(p in 1:length(score.null$results)) {
      min.pvl[p] <- min(score.null$results[[p]]$pval)
    }        
  } # else if(!is.null(model$min.pvl)) {
  #   min.pvl <- model$min.pvl
  # }
  
  if(!is.null(plot)) plot <- paste(plot, "pdf", sep = ".")
  results <- vector("list", length(model$results))
  names(results) <- names(model$results)
  if(data$step >= 1) w.size <- w.size/data$step
  
  for(p in 1:length(results)) {
    
    if(!is.null(min.pvl)) {
      sig.fwd <- quantile(sort(min.pvl), sig.fwd0)#; cat(sig.fwd, "\n")
    } else {
      sig.fwd <- sig.fwd0
    }
    
    round <- 1
    start <- proc.time()
    pheno.col <- model$results[[p]]$pheno.col
    ind <- rownames(data$pheno)[which(!is.na(data$pheno[,pheno.col]))]
    Y <- data$pheno[ind,pheno.col]
    if(is.null(offset.data)) {
      offset <- NULL
    } else {
      offset <- offset.data[ind,pheno.col]
    }
    stat <- model$results[[p]]$stat
    pval <- model$results[[p]]$pval
    temp <- rbind(st = stat, pv = pval)
    colnames(temp) <- as.character(1:data$nmrk)
    if(is.null(model$results[[p]]$qtls)) {
      if(verbose) cat("Forward search for trait ", pheno.col, " ", sQuote(colnames(data$pheno)[pheno.col]), "; there are no QTL in the model \n", sep="")
      if(!is.null(plot)) pdf(paste(colnames(data$pheno)[pheno.col], plot, sep = "_"))
      round <- 0
      qtl.mrk <- c()
      qtl.lgr <- c()
      qtl.pos <- c()
    } else {
      qtl.mrk <- model$results[[p]]$qtl[,"Nmrk"]
      qtl.lgr <- model$results[[p]]$qtl[,"LG"]
      qtl.pos <- model$results[[p]]$qtl[,"Pos"]
      if(verbose) {
        if(length(qtl.mrk) == 1) cat("Forward search for trait ", pheno.col, " ", sQuote(colnames(data$pheno)[pheno.col]), "; there is ", length(qtl.mrk), " QTL in the model already \n", sep="")
        if(length(qtl.mrk) >= 2) cat("Forward search for trait ", pheno.col, " ", sQuote(colnames(data$pheno)[pheno.col]), "; there are ", length(qtl.mrk), " QTL in the model already \n", sep="")
      }
      if(!is.null(plot)) pdf(paste(colnames(data$pheno)[pheno.col], plot, sep = "_"))
      qtl.vcv <- NULL
      markers.out <- c()
      interval <- c()
      for(q in 1:length(qtl.mrk)) {
        qtl.vcv <- c(qtl.vcv, list(data$G[ind,ind,qtl.mrk[q]]))
        interval <- c((qtl.mrk[q]-w.size):(qtl.mrk[q]+w.size))
        markers.out <- c(markers.out, interval[c(which(interval >= (data$cum.nmrk[qtl.lgr[q]]+1) & interval <= (data$cum.nmrk[qtl.lgr[q]+1])))])
      }
      markers <- c(1:data$nmrk)[-markers.out]
      if(polygenes) {
        Gstar <- apply(data$G[ind,ind,qtl.mrk], MARGIN = c(1,2), sum)/length(qtl.mrk); Gstar[1:5,1:5]
        full.mod0 <- varComp(Y ~ 1, varcov = list(Gstar), offset = offset)
        control <- varComp.control(start = c(coef(full.mod0, what = "var.ratio"),0))
        temp <- parSapply(cl, as.character(markers), function(x) {
          m <- as.numeric(x)
          full.mod <- varComp(Y ~ 1, varcov = list(Gstar, data$G[ind,ind,m]), control = control, offset = offset)
          test <- varComp.test(full.mod, null=1L)
          c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
        })
      } else {
        withCallingHandlers(full.mod0 <- varComp(Y ~ 1, varcov = c(qtl.vcv), offset = offset), warning = h)
        control <- varComp.control(start = c(coef(full.mod0, what = "var.ratio"),0))
        temp <- parSapply(cl, as.character(markers), function(x) {
          m <- as.numeric(x)
          full.mod <- varComp(Y ~ 1, varcov = c(qtl.vcv, list(data$G[ind,ind,m])), control = control, offset = offset)
          test <- varComp.test(full.mod, null=c(1:length(qtl.vcv)))
          c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
        })
      }
      if(!is.null(plot)) {
        plot(x=as.numeric(names(temp["st",])), y=-log10(temp["pv",]), xlab="Position number", ylab="-log10(p)", main=paste("Search round #", round, sep=""), ylim=c(0,10))
        abline(v=data$cum.nmrk, lty=3); abline(h=-log10(sig.fwd), lty=5); points(x=qtl.mrk, y=rep(-0.15, length(qtl.mrk)), pch=6, lwd=1.5, col="red")
      }
    }
    
    while(temp["pv",which.max(temp["st",])] <= sig.fwd & round < n.rounds) {
      round <- round + 1
      qtl.mrk <- c(qtl.mrk,  as.numeric(names(which.max(temp["st",]))))
      qtl.lgr <- c(qtl.lgr, last(which(last(qtl.mrk) > data$cum.nmrk)))
      qtl.pos <- c(qtl.pos, round(unlist(data$lgs)[[last(qtl.mrk)]], digits = 2))
      if(verbose) cat("  QTL was found on LG ", last(qtl.lgr), " at ", last(qtl.pos), " cM (position number ", last(qtl.mrk), ")\n", sep="")
      qtl.vcv <- NULL
      markers.out <- c()
      interval <- c()
      for(q in 1:length(qtl.mrk)) {
        qtl.vcv <- c(qtl.vcv, list(data$G[ind,ind,qtl.mrk[q]]))
        interval <- c((qtl.mrk[q]-w.size):(qtl.mrk[q]+w.size))
        markers.out <- c(markers.out, interval[c(which(interval >= (data$cum.nmrk[qtl.lgr[q]]+1) & interval <= (data$cum.nmrk[qtl.lgr[q]+1])))])
      }
      markers <- c(1:data$nmrk)[-markers.out]
      if(polygenes) {
        Gstar <- apply(data$G[ind,ind,qtl.mrk], MARGIN = c(1,2), sum)/length(qtl.mrk); Gstar[1:5,1:5]
        full.mod0 <- varComp(Y ~ 1, varcov = list(Gstar), offset = offset)
        control <- varComp.control(start = c(coef(full.mod0, what = "var.ratio"),0))
        temp <- parSapply(cl, as.character(markers), function(x) {
          m <- as.numeric(x)
          full.mod <- varComp(Y ~ 1, varcov = list(Gstar, data$G[ind,ind,m]), control = control, offset = offset)
          test <- varComp.test(full.mod, null=1L)
          c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
        })
      } else {
        withCallingHandlers(full.mod0 <- varComp(Y ~ 1, varcov = c(qtl.vcv), offset = offset), warning = h)
        control <- varComp.control(start = c(coef(full.mod0, what = "var.ratio"),0))
        temp <- parSapply(cl, as.character(markers), function(x) {
          m <- as.numeric(x)
          full.mod <- varComp(Y ~ 1, varcov = c(qtl.vcv, list(data$G[ind,ind,m])), control = control, offset = offset)
          test <- varComp.test(full.mod, null=c(1:length(qtl.vcv)))
          c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
        })
      }
      if(!is.null(plot)) {
        round <- 1+round
        plot(x=as.numeric(names(temp["st",])), y=-log10(temp["pv",]), xlab="Position number", ylab="-log10(p)", main=paste("Search round #", round, sep=""), ylim=c(0,10))
        abline(v=data$cum.nmrk, lty=3); abline(h=-log10(sig.fwd), lty=5); points(x=qtl.mrk, y=rep(-0.15, length(qtl.mrk)), pch=6, lwd=1.5, col="red")
      }
    } # QTL search while significant p-values: forward search
    
    qtl.mrk0 <- as.numeric(names(which.max(temp["st",])))
    qtl.lgr0 <- last(which(last(qtl.mrk0) > data$cum.nmrk))
    qtl.pos0 <- round(unlist(data$lgs)[[qtl.mrk0]], digits = 2)
    if(!is.null(qtl.mrk) & verbose) cat("  No more QTL were found. A putative QTL on LG ", last(qtl.lgr0), " at ", last(qtl.pos0), " cM (position number ", last(qtl.mrk0), ") did not reach the threshold; its p-value was ", round(temp["pv",][which.max(temp["st",])], 5), "\n", sep="")
    stat[as.numeric(colnames(temp))] <- temp["st",]
    pval[as.numeric(colnames(temp))] <- temp["pv",]
    
    if(is.null(qtl.mrk)) {
      qtl.mrk0 <- which.max(stat)
      qtl.lgr0 <- last(which(last(qtl.mrk0) > data$cum.nmrk))
      qtl.pos0 <- round(unlist(data$lgs)[[last(qtl.mrk0)]], digits = 2)
      if(verbose) cat("  No QTL were found. A putative QTL on LG ", last(qtl.lgr0), " at ", last(qtl.pos0), " cM (position number ", last(qtl.mrk0), ") did not reach the threshold; its p-value was ", round(temp["pv",][which.max(temp["st",])], 5), "\n", sep="")
    }
    
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
    }
    
    results[[p]] <- list(
      pheno.col=pheno.col,
      stat=stat,
      pval=pval,
      qtls=qtls)
    
  }
  
  stopCluster(cl)
  
  structure(list(data=deparse(substitute(data)),
                 offset.data=deparse(substitute(offset.data)),
                 pheno.col=model$pheno.col,
                 w.size=w.size,
                 sig.fwd=sig.fwd0,
                 sig.bwd=NULL,
                 min.pvl=min.pvl,
                 polygenes=polygenes,
                 d.sint=NULL,
                 results=results),
            class=c("qtlpoly.model","qtlpoly.search"))
  
}

#' @rdname search_qtl
#' @export

print.qtlpoly.search <- function(x, pheno.col = NULL) {
  if(any(class(x) == "qtlpoly.search")) cat("This is an object of class 'qtlpoly.search'\n")
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
