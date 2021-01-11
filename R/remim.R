#' Random-effect multiple interval mapping (REMIM)
#'
#' Automatic function that performs REMIM algorithm using score statistics.
#'
#' @param data an object of class \code{qtlpoly.data}.
#'
#' @param pheno.col a numeric vector with the phenotype columns to be analyzed or printed; if \code{NULL} (default), all phenotypes from \code{'data'} will be included.
#'
#' @param w.size the window size (in centiMorgans) to avoid on either side of QTL already in the model when looking for a new QTL, e.g. 15 (default).
#'
#' @param sig.fwd the desired score-based significance level for forward search, e.g. 0.01 (default).
#'
#' @param sig.bwd the desired score-based significance level for backward elimination, e.g. 0.001 (default).
#'
#' @param score.null an object of class \code{qtlpoly.null} with results of score statistics from resampling.
#'
#' @param d.sint a \eqn{d} value to subtract from logarithm of \emph{p}-value (\eqn{LOP-d}) for support interval calculation, e.g. \eqn{d=1.5} (default) represents approximate 95\% support interval.
#'
#' @param polygenes if \code{TRUE} all QTL already in the model are treated as a single polygenic effect; if \code{FALSE} (default) all QTL effect variances have to estimated.
#'
#' @param n.clusters number of parallel processes to spawn.
#'
#' @param n.rounds number of search rounds; if \code{Inf} (default) forward search will stop when no more significant positions can be found.
#'
#' @param plot a suffix for the file's name containing plots of every algorithm step, e.g. "remim" (default); if \code{NULL}, no file is produced.
#'
#' @param verbose if \code{TRUE} (default), current progress is shown; if \code{FALSE}, no output is produced.
#'
#' @param x an object of class \code{qtlpoly.remim} to be printed.
#'
#' @param sint whether \code{"upper"} or \code{"lower"} support intervals should be printed; if \code{NULL} (default), only QTL peak information will be printed.
#'
#' @return An object of class \code{qtlpoly.remim} which contains a list of \code{results} for each trait with the following components:
#'
#'     \item{pheno.col}{a phenotype column number.}
#'     \item{stat}{a vector containing values from score statistics.}
#'     \item{pval}{a vector containing \emph{p}-values from score statistics.}
#'     \item{qtls}{a data frame with information from the mapped QTL.}
#'     \item{lower}{a data frame with information from the lower support interval of mapped QTL.}
#'     \item{upper}{a data frame with information from the upper support interval of mapped QTL.}
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
#'   # perform remim
#'   remim.mod <- remim(data = data, w.size = 15, sig.fwd = 0.01, sig.bwd = 0.0001,
#'     d.sint = 1.5, n.clusters = 4, plot = "remim")
#'   }
#'
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}
#'
#' @references
#'     Kao CH, Zeng ZB, Teasdale RD (1999) Multiple interval mapping for quantitative trait loci. \emph{Genetics} 152 (3): 1203–16. \url{www.genetics.org/content/152/3/1203}. 
#' 
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2020) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{Genetics} 215 (3): 579-595. \url{http://doi.org/10.1534/genetics.120.303080}.
#'     
#'     Qu L, Guennel T, Marshall SL (2013) Linear score tests for variance components in linear mixed models and applications to genetic association studies. \emph{Biometrics} 69 (4): 883–92. \url{doi.org/10.1111/biom.12095}.
#'
#'     Zou F, Fine JP, Hu J, Lin DY (2004) An efficient resampling method for assessing genome-wide statistical significance in mapping quantitative trait loci. \emph{Genetics} 168 (4): 2307-16. \url{doi.org/10.1534/genetics.104.031427}
#'
#' @export remim

remim <- function(data, pheno.col = NULL, w.size = 15, sig.fwd = 0.01, sig.bwd = 0.0001, score.null = NULL, d.sint = 1.5, polygenes = FALSE, n.clusters = NULL, n.rounds = Inf, plot = "remim", verbose = TRUE) {
  
  if(is.null(n.clusters)) n.clusters <- 1
  cat("INFO: Using", n.clusters, "CPUs for calculation\n\n")
  cl <- makeCluster(n.clusters)
  clusterEvalQ(cl, require(varComp))
  sig.fwd0 <- sig.fwd
  sig.bwd0 <- sig.bwd
  
  min.pvl <- NULL
  if(!is.null(score.null)) {
    min.pvl <- numeric(length(score.null$results))
    for(p in 1:length(score.null$results)) {
      min.pvl[p] <- score.null$results[[p]]$pval[which.max(score.null$results[[p]]$stat)]
    }        
  } 
  
  if(is.null(pheno.col)) pheno.col <- 1:dim(data$pheno)[2]
  if(!is.null(plot)) plot <- paste(plot, "pdf", sep = ".")
  results <- vector("list", length(pheno.col))
  names(results) <- colnames(data$pheno)[pheno.col]
  if(data$step > 1) w.size <- w.size/data$step
  
  for(p in 1:length(results)) {
    
    round <- 1
    stat <- numeric(data$nmrk)
    pval <- numeric(data$nmrk)
    start <- proc.time()
    if(verbose) cat("REMIM for trait", pheno.col[p], sQuote(colnames(data$pheno)[pheno.col[p]]), "\n")
    # if(!is.null(plot)) pdf(paste(colnames(data$pheno)[pheno.col], ".pdf", sep = ""))
    if(!is.null(plot)) pdf(paste(colnames(data$pheno)[pheno.col[p]], plot, sep = "_"))
    if(!is.null(min.pvl)) {
      sig.fwd <- quantile(sort(min.pvl), sig.fwd0); # cat(sig.fwd, "\n")
      sig.bwd <- quantile(sort(min.pvl), sig.bwd0); # cat(sig.bwd, "\n")
    } else {
      sig.fwd <- sig.fwd0
      sig.bwd <- sig.bwd0
    }
    
    ind <- rownames(data$pheno)[which(!is.na(data$pheno[,pheno.col[p]]))]
    Y <- data$pheno[ind,pheno.col[p]]
    if(!is.null(data$weights)) weight <- data$weights[ind,pheno.col[p]] else weight <- rep(1,length(ind))
    markers <- c(1:data$nmrk)
    temp <- parSapply(cl, as.character(markers), function(x) {
      m <- as.numeric(x)
      full.mod <- varComp(Y ~ 1, varcov = list(data$G[ind,ind,m]), weights = weight/max(weight))
      test <- varComp.test(full.mod, null=integer(0L))
      c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
    }) #first search
    if(!is.null(plot)) {
      search <- 1
      plot(x=as.numeric(names(temp["st",])), y=-log10(temp["pv",]), xlab="Position number", ylab="-log10(p)", main=paste("Search round #", search, sep=""), ylim=c(0,10))
      abline(v=data$cum.nmrk, lty=3); abline(h=-log10(sig.fwd), lty=5)
    }
    qtl.mrk <- c()
    qtl.lgr <- c()
    qtl.pos <- c()
    qtl.out0 <- c()
    
    if(temp["pv",which.max(temp["st",])] > sig.fwd) {
      qtl.mrk0 <- as.numeric(names(which.max(temp["st",])))
      qtl.lgr0 <- last(which(last(qtl.mrk0) > data$cum.nmrk))
      qtl.pos0 <- round(unlist(data$lgs)[[qtl.mrk0]], digits = 2)
      if(verbose) cat("  No QTL were found. A putative QTL on LG ", last(qtl.lgr0), " at ", last(qtl.pos0), " cM (position number ", last(qtl.mrk0), ") did not reach the threshold; its p-value was ", round(temp["pv",][which.max(temp["st",])], 5), "\n", sep="")
      stat[as.numeric(colnames(temp))] <- temp["st",]
      pval[as.numeric(colnames(temp))] <- temp["pv",]
    }
    
    while(temp["pv",which.max(temp["st",])] <= sig.fwd & round < n.rounds & !any(as.numeric(names(which.max(temp["st",]))) == qtl.mrk)) { # QTL search while significant p-values
      while(temp["pv",which.max(temp["st",])] <= sig.fwd) {
        qtl.mrk <- c(qtl.mrk, as.numeric(names(which.max(temp["st",]))))
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
          full.mod0 <- varComp(Y ~ 1, varcov = list(Gstar), weights = weight/max(weight))
          control <- varComp.control(start = c(coef(full.mod0, what = "var.ratio"),0))
          temp <- parSapply(cl, as.character(markers), function(x) {
            m <- as.numeric(x)
            full.mod <- varComp(Y ~ 1, varcov = list(Gstar, data$G[ind,ind,m]), control = control, weights = weight/max(weight))
            test <- varComp.test(full.mod, null=1L)
            c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
          })
        } else {
          withCallingHandlers(full.mod0 <- varComp(Y ~ 1, varcov = c(qtl.vcv), weights = weight/max(weight)), warning = h)
          control <- varComp.control(start = c(coef(full.mod0, what = "var.ratio"),0))
          temp <- parSapply(cl, as.character(markers), function(x) {
            m <- as.numeric(x)
            full.mod <- varComp(Y ~ 1, varcov = c(qtl.vcv, list(data$G[ind,ind,m])), control = control, weights = weight/max(weight))
            test <- varComp.test(full.mod, null=c(1:length(qtl.vcv)))
            c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
          })
        }
        if(!is.null(plot)) {
          search <- 1+search
          plot(x=as.numeric(names(temp["st",])), y=-log10(temp["pv",]), xlab="Position number", ylab="-log10(p)", main=paste("Search round #", search, sep=""), ylim=c(0,10), xlim = c(0,last(data$cum.nmrk)))
          abline(v=data$cum.nmrk, lty=3); abline(h=-log10(sig.fwd), lty=5); points(x=qtl.mrk, y=rep(-0.15, length(qtl.mrk)), pch=6, lwd=1.5, col="red")
        }
      } # QTL search while significant p-values: forward serch
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
      qtl.out <- c(1)
      qtl.out1 <- c()
      while(!is.null(qtl.out)) {
        qtl.out <- c()
        qtl.mrk1 <- qtl.mrk
        if(length(qtl.mrk) == 1) {
          if(verbose) cat("  Refining QTL positions ...", qtl.mrk, "\n")
          markers.out <- c((data$cum.nmrk[qtl.lgr[1]]+1):(data$cum.nmrk[qtl.lgr[1]+1]))
          temp <- parSapply(cl, as.character(markers.out), function(x) { #like first search
            m <- as.numeric(x)
            full.mod <- varComp(Y ~ 1, varcov = list(data$G[ind,ind,m]), weights = weight/max(weight))
            test <- varComp.test(full.mod, null=integer(0L))
            c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
          })
          stat[as.numeric(colnames(temp))] <- temp["st",]
          pval[as.numeric(colnames(temp))] <- temp["pv",]
          if(pval[markers.out[which.max(stat[markers.out])]] <= sig.bwd) { # updates position
            qtl.mrk[1] <- markers.out[which.max(stat[markers.out])]
            qtl.pos[1] <- round(unlist(data$lgs)[[qtl.mrk[1]]], digits = 2)
          } else { # stores non-significant
            qtl.out <- c(1)
          }
          if(!is.null(plot)) {
            plot(-log10(pval), xlab="Position number", ylab="-log10(p)", main="Refining round #1", ylim=c(0,10))
            abline(v=data$cum.nmrk, lty=3); abline(h=-log10(sig.bwd), lty=5); points(x=qtl.mrk, y=rep(-0.15, length(qtl.mrk)), pch=6, lwd=1.5, col="red")
            if(!is.null(qtl.out)) points(x=qtl.mrk[qtl.out], y=rep(-0.15, length(qtl.out)), pch=4, lwd=1.5, col="red")
          }
          if(!is.null(qtl.out)) {
            if(verbose) cat("  Excluding non-significant QTL", paste("...", qtl.mrk[qtl.out]), "\n")
            qtl.mrk <- qtl.mrk[-qtl.out]
            qtl.lgr <- qtl.lgr[-qtl.out]
            qtl.pos <- qtl.pos[-qtl.out]
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
                full.mod0 <- varComp(Y ~ 1, varcov = list(Gstar), weights = weight/max(weight))
                control <- varComp.control(start = c(coef(full.mod0, what = "var.ratio"),0))
                temp <- parSapply(cl, as.character(markers.out), function(x) {
                  m <- as.numeric(x)
                  full.mod <- varComp(Y ~ 1, varcov = list(Gstar, data$G[ind,ind,m]), control = control, weights = weight/max(weight))
                  test <- varComp.test(full.mod, null=1L)
                  c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
                })
              } else {
                withCallingHandlers(full.mod0 <- varComp(Y ~ 1, varcov = c(qtl.vcv), weights = weight/max(weight)), warning = h)
                control <- varComp.control(start = c(coef(full.mod0, what = "var.ratio"),0))
                temp <- parSapply(cl, as.character(markers.out), function(x) {
                  m <- as.numeric(x)
                  full.mod <- varComp(Y ~ 1, varcov = c(qtl.vcv, list(data$G[ind,ind,m])), control = control, weights = weight/max(weight))
                  test <- varComp.test(full.mod, null=c(1:length(qtl.vcv)))
                  c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
                })
              }
              stat[as.numeric(colnames(temp))] <- temp["st",]
              pval[as.numeric(colnames(temp))] <- temp["pv",]
              if(pval[markers.out[which.max(stat[markers.out])]] <= sig.bwd) { # updates position
                qtl.mrk[q] <- markers.out[which.max(stat[markers.out])]
                qtl.pos[q] <- round(unlist(data$lgs)[[qtl.mrk[q]]], digits = 2)
              } else { # stores non-significant
                qtl.out <- c(qtl.out, q)
              }
              if(!is.null(plot)) {
                plot(-log10(pval), xlab="Position number", ylab="-log10(p)", main=paste("Refining round #", q, sep=""), ylim=c(0,10))
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
            qtl.out1 <- qtl.mrk[qtl.out]
            if(verbose) cat("  Excluding non-significant QTL", paste("...", qtl.out1), "\n")
            # if(!is.null(plot)) {
            #   plot(-log10(pval), xlab="Position number", ylab="-log10(p)", main=paste("Refining round #", q, sep=""), ylim=c(0,10))
            #   abline(v=data$cum.nmrk, lty=3); abline(h=-log10(sig.bwd), lty=5); points(x=qtl.mrk, y=rep(-0.15, length(qtl.mrk)), pch=6, lwd=1.5, col="red")
            #   if(!is.null(qtl.out)) points(x=qtl.mrk[qtl.out], y=rep(-0.15, length(qtl.out)), pch=4, lwd=1.5, col="red")
            # }
            qtl.mrk <- qtl.mrk[-qtl.out]
            qtl.lgr <- qtl.lgr[-qtl.out]
            qtl.pos <- qtl.pos[-qtl.out]
          }
        }
        if(is.null(qtl.out) & length(qtl.mrk1) == length(qtl.mrk)) { # if QTL position changes a lot, significance may change too
          if(any(abs(qtl.mrk1 - qtl.mrk) > w.size)) qtl.out <- c(1)
        }
      } # keeps refining until all QTL are significant
      
      lower <- upper <- qtl.mrk
      
      if(length(qtl.mrk) == 0) {
        markers <- c(1:data$nmrk)
        temp <- parSapply(cl, as.character(markers), function(x) { #like first search
          m <- as.numeric(x)
          full.mod <- varComp(Y ~ 1, varcov = list(data$G[ind,ind,m]), weights = weight/max(weight))
          test <- varComp.test(full.mod, null=integer(0L))
          c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
        })
        stat[as.numeric(colnames(temp))] <- temp["st",]
        pval[as.numeric(colnames(temp))] <- temp["pv",]
        if(!is.null(plot)) {
          plot(-log10(pval), xlab="Position number", ylab="-log10(p)", main="Completing genome", ylim=c(0,10))
          abline(v=data$cum.nmrk, lty=3); abline(h=-log10(sig.bwd), lty=5)
        }
      } # completing genome when there's no QTL
      
      if(length(qtl.mrk) == 1) {
        if(verbose) cat("  Profiling QTL ...", qtl.mrk, "\n")
        markers.out <- c((data$cum.nmrk[qtl.lgr[1]]+1):(data$cum.nmrk[qtl.lgr[1]+1]))
        markers <- c(1:data$nmrk)[-markers.out]
        temp <- parSapply(cl, as.character(markers.out), function(x) { #like first search
          m <- as.numeric(x)
          full.mod <- varComp(Y ~ 1, varcov = list(data$G[ind,ind,m]), weights = weight/max(weight))
          test <- varComp.test(full.mod, null=integer(0L))
          c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
        })
        stat[as.numeric(colnames(temp))] <- temp["st",]
        pval[as.numeric(colnames(temp))] <- temp["pv",]
        qtl.vcv <- list(data$G[ind,ind,qtl.mrk[1]])
        withCallingHandlers(full.mod0 <- varComp(Y ~ 1, varcov = c(qtl.vcv), weights = weight/max(weight)), warning = h)
        control <- varComp.control(start = c(coef(full.mod0, what = "var.ratio"),0))
        temp <- parSapply(cl, as.character(markers), function(x) { #markers outside chr with QTL (second search)
          m <- as.numeric(x)
          full.mod <- varComp(Y ~ 1, varcov = c(qtl.vcv, list(data$G[ind,ind,m])), control = control, weights = weight/max(weight))
          test <- varComp.test(full.mod, null=1L)
          c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
        })
        stat[as.numeric(colnames(temp))] <- temp["st",]
        pval[as.numeric(colnames(temp))] <- temp["pv",]
        if(!is.null(plot)) {
          plot(-log10(pval), xlab="Position number", ylab="-log10(p)", main="Profiling round #1", ylim=c(0,10))
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
            full.mod0 <- varComp(Y ~ 1, varcov = list(Gstar), weights = weight/max(weight))
            control <- varComp.control(start = c(coef(full.mod0, what = "var.ratio"),0))
            temp <- parSapply(cl, as.character(markers.out), function(x) {
              m <- as.numeric(x)
              full.mod <- varComp(Y ~ 1, varcov = list(Gstar, data$G[ind,ind,m]), control = control, weights = weight/max(weight))
              test <- varComp.test(full.mod, null=1L)
              c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
            })
          } else {
            withCallingHandlers(full.mod0 <- varComp(Y ~ 1, varcov = c(qtl.vcv), weights = weight/max(weight)), warning = h)
            control <- varComp.control(start = c(coef(full.mod0, what = "var.ratio"),0))
            temp <- parSapply(cl, as.character(markers.out), function(x) {
              m <- as.numeric(x)
              full.mod <- varComp(Y ~ 1, varcov = c(qtl.vcv, list(data$G[ind,ind,m])), control = control, weights = weight/max(weight))
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
            plot(-log10(pval), xlab="Position number", ylab="-log10(p)", main=paste("Profiling round #", q, sep=""), ylim=c(0,10))
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
            full.mod0 <- varComp(Y ~ 1, varcov = list(Gstar), weights = weight/max(weight))
            control <- varComp.control(start = c(coef(full.mod0, what = "var.ratio"),0))
            temp <- parSapply(cl, as.character(markers.out), function(x) {
              m <- as.numeric(x)
              full.mod <- varComp(Y ~ 1, varcov = list(Gstar, data$G[ind,ind,m]), control = control, weights = weight/max(weight))
              test <- varComp.test(full.mod, null=1L)
              c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
            })
          } else {
            withCallingHandlers(full.mod <- varComp(Y ~ 1, varcov = c(qtl.vcv), weights = weight/max(weight)), warning = h)
            control <- varComp.control(start = c(coef(full.mod, what = "var.ratio"),0))
            temp <- parSapply(cl, as.character(markers.out), function(x) {
              m <- as.numeric(x)
              full.mod <- varComp(Y ~ 1, varcov = c(qtl.vcv, list(data$G[ind,ind,m])), control = control, weights = weight/max(weight))
              test <- varComp.test(full.mod, null=c(1:length(qtl.vcv)))
              c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
            })
          }
          stat[as.numeric(colnames(temp))] <- temp["st",]
          pval[as.numeric(colnames(temp))] <- temp["pv",]
          if(!is.null(plot)) {
            plot(-log10(pval), xlab="Position number", ylab="-log10(p)", main="Completing genome", ylim=c(0,10))
            abline(v=data$cum.nmrk, lty=3); points(x=qtl.mrk, y=rep(-0.15, length(qtl.mrk)), pch=6, lwd=1.5, col="red"); abline(h=-log10(sig.bwd), lty=5)
          }
        }
      } # profiling 1+ QTL
      
      sig.fwd <- sig.bwd # sig.fwd is updated to the last sig.bwd
      
      if(!is.null(qtl.out0) & !is.null(qtl.out1)) {
        if(qtl.out0 == qtl.out1) break
      } # breaks while if adds the same one excluded
      qtl.out0 <- qtl.out1
      round <- round + 1
    } # keeps doing forward, refinement and backward
    
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
    
    # sig.fwd <- sig.fwd0
    
    results[[p]] <- list(
      pheno.col=pheno.col[p],
      stat=stat,
      pval=pval,
      qtls=qtls,
      lower=lower,
      upper=upper)
    
  }
  
  stopCluster(cl)
  
  structure(list(data=deparse(substitute(data)),
                 pheno.col=pheno.col,
                 w.size=w.size*data$step,
                 sig.fwd=sig.fwd0,
                 sig.bwd=sig.bwd0,
                 min.pvl=min.pvl,
                 polygenes=polygenes,
                 d.sint=d.sint,
                 results=results),
            class=c("qtlpoly.model","qtlpoly.remim"))
  
}

#' @rdname remim
#' @export
print.qtlpoly.remim <- function(x, pheno.col = NULL, sint=NULL) {
  if(any(class(x) == "qtlpoly.remim")) cat("This is an object of class 'qtlpoly.remim'\n")
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
