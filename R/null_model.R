#' Null model
#'
#' Creates a null model (with no QTL) for each trait.
#'
#' @param data an object of class \code{qtlpoly.data}.
#'
#' @param offset.data a data frame with the same dimensions of \code{data$pheno} containing offset variables; if \code{NULL} (default), no offset variables are considered.
#' 
#' @param pheno.col a numeric vector with the phenotype columns to be analyzed; if \code{NULL}, all phenotypes from \code{'data'} will be included.
#'
#' @param n.clusters number of parallel processes to spawn.
#'
#' @param plot a suffix for the file's name containing simple plots of every QTL search round, e.g. "null" (default); if \code{NULL}, no file is produced.
#'
#' @param verbose if \code{TRUE} (default), current progress is shown; if \code{FALSE}, no output is produced.
#'
#' @param x an object of class \code{qtlpoly.null} to be printed.
#'
#' @return An object of class \code{qtlpoly.null} which contains a list of \code{results} for each trait with the following components:
#'
#'     \item{pheno.col}{a phenotype column number.}
#'     \item{stat}{a vector containing values from score statistics.}
#'     \item{pval}{a vector containing \emph{p}-values from score statistics.}
#'     \item{qtls}{a data frame with information from the mapped QTL (\code{NULL} at this point).}
#'
#' @seealso \code{\link[qtlpoly]{read_data}}
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
#'   }
#'
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2019) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{bioRxiv}. \url{doi.org/10.1101/622951}.
#'     
#'     Qu L, Guennel T, Marshall SL (2013) Linear score tests for variance components in linear mixed models and applications to genetic association studies. \emph{Biometrics} 69 (4): 883–92. \url{doi.org/10.1111/biom.12095}.
#'
#' @export null_model
#' @import varComp parallel

null_model <- function(data, offset.data = NULL, pheno.col = NULL, n.clusters = NULL, plot = "null", verbose = TRUE) {
  
  if(is.null(n.clusters)) n.clusters <- 1
  cat("INFO: Using", n.clusters, "CPUs for calculation\n\n")
  cl <- makeCluster(n.clusters)
  clusterEvalQ(cl, require(varComp))
  
  if(is.null(pheno.col)) pheno.col <- 1:dim(data$pheno)[2]
  if(!is.null(plot)) plot <- paste(plot, "pdf", sep = ".")
  results <- vector("list", length(pheno.col))
  names(results) <- colnames(data$pheno)[pheno.col]
  
  for(p in 1:length(results)) {
    
    start <- proc.time()
    stat <- numeric(data$nmrk)
    pval <- numeric(data$nmrk)
    ind <- rownames(data$pheno)[which(!is.na(data$pheno[,pheno.col[p]]))]
    Y <- data$pheno[ind,pheno.col[p]]
    if(is.null(offset.data)) {
      offset <- NULL
    } else {
      offset <- offset.data[ind,pheno.col[p]]
    }
    if(verbose) cat("Null model for trait", pheno.col[p], sQuote(colnames(data$pheno)[pheno.col[p]]), "\n")
    if(!is.null(plot)) pdf(paste(colnames(data$pheno)[pheno.col[p]], plot, sep = "_"))
    
    markers <- c(1:data$nmrk)
    temp <- parSapply(cl, as.character(markers), function(x) {
      m <- as.numeric(x)
      full.mod <- varComp(Y ~ 1, varcov = list(data$G[ind,ind,m]), offset = offset)
      test <- varComp.test(full.mod, null=integer(0L))
      c(st=as.numeric(test[[1]][[1]][[1]]$statistic), pv=as.numeric(test[[1]][[1]][[1]]$p.value))
    })
    if(!is.null(plot)) {
      round <- 1
      plot(x=as.numeric(names(temp["st",])), y=-log10(temp["pv",]), xlab="Marker number", ylab="-log10(p)", main="Null model", ylim=c(0,10))
      abline(v=data$cum.nmrk, lty=3)
    }
    stat[as.numeric(colnames(temp))] <- temp["st",]
    pval[as.numeric(colnames(temp))] <- temp["pv",]
    
    results[[p]] <- list(
      pheno.col=pheno.col[p],
      stat=stat,
      pval=pval,
      qtls=NULL)
    
    if(!is.null(plot)) dev.off()
    end <- proc.time()
    if(verbose) cat("  Calculation took", round((end - start)[3], digits = 2), "seconds\n\n")
    
  }
  
  stopCluster(cl)
  
  structure(list(data=deparse(substitute(data)),
                 offset.data=deparse(substitute(offset.data)),
                 pheno.col=pheno.col,
                 w.size=NULL,
                 sig.fwd=NULL,
                 sig.bwd=NULL,
                 polygenes=NULL,
                 d.sint=NULL,
                 results=results),
            class=c("qtlpoly.model","qtlpoly.null"))
  
}

#' @rdname null_model
#' @export

print.qtlpoly.null <- function(x, pheno.col = NULL) {
  if(any(class(x) == "qtlpoly.null")) cat("This is an object of class 'qtlpoly.null'\n")
  if(is.null(pheno.col)) {
    pheno.col <- 1:length(x$results)
  } else {
    pheno.col <- which(x$pheno.col %in% pheno.col)
  }
  for(p in pheno.col) {
    cat("\n* Trait", x$results[[p]]$pheno.col, sQuote(names(x$results)[[p]]), "\n")
    cat("There are no QTL in the model \n")
  }
}
