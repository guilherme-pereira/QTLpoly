#' QTL allele effect estimation
#'
#' Computes allele specific and allele combination additive effects from multiple QTL models.
#'
#' @param ploidy a numeric value of ploidy level of the cross (currently, only 4 or 6).
#'
#' @param fitted a fitted multiple QTL model of class \code{qtlpoly.fitted}.
#'
#' @param x an object of class \code{qtlpoly.effects} to be plotted.
#'
#' @param pheno.col a numeric vector with the phenotype column numbers to be plotted; if \code{NULL}, all phenotypes from \code{'fitted'} will be included.
#'
#' @param p1 a character string with the first parent name, e.g. \code{"P1"} (default).
#'
#' @param p2 a character string with the second parent name, e.g. \code{"P2"} (default).
#'
#' @return An object of class \code{qtlpoly.effects} which is a list of \code{results} for each containing the following components:
#'
#'     \item{pheno.col}{a phenotype column number.}
#'     \item{y.hat}{a vector with the predicted values.}
#'
#' @return A \pkg{ggplot2} barplot with parental allele and allele combination effects.
#'
#' @seealso \code{\link[qtlpoly]{read_data}}, \code{\link[qtlpoly]{remim}}, \code{\link[qtlpoly]{fit_model}}
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
#'   # fit model
#'   fitted.mod <- fit_model(data=data, model=remim.mod, probs="joint", polygenes="none")
#'
#'   # estimate effects
#'   est.effects <- qtl_effects(ploidy = 6, fitted = fitted.mod)
#'   plot(est.effects)
#'   }
#'   
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2019) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{bioRxiv}. \url{doi:}.
#'
#' @export qtl_effects

qtl_effects <- function(ploidy = 6, fitted) {

  results <- vector("list", length(fitted$results))
  names(results) <- names(fitted$results)

  for(p in 1:length(results)) {

    if(!is.null(fitted$results[[p]]$qtls)) {

      nqtl <- dim(fitted$results[[p]]$qtls)[1]
      if(nqtl > 1) nqtl <- nqtl - 1
      effects <- vector("list", nqtl)

      if(ploidy == 6) {

        for(q in 1:nqtl) {

          blups <- fitted$results[[p]]$fitted$u.hat[[q]]
          alleles <- matrix(unlist(strsplit(rownames(blups), '')), ncol=7, byrow=TRUE)[,-4]

          A <- t(combn(letters[1:12],1))
          D <- t(combn(letters[1:12],2))
          T <- t(combn(letters[1:12],3))
          F <- t(combn(letters[1:12],4))
          G <- t(combn(letters[1:12],5))
          S <- t(combn(letters[1:12],6))

          a <- vector("list", dim(A)[1])
          d <- vector("list", dim(D)[1])
          t <- vector("list", dim(T)[1])
          f <- vector("list", dim(F)[1])
          g <- vector("list", dim(G)[1])
          s <- vector("list", dim(S)[1])

          for(i in 1:dim(A)[1]) {
            a[[i]] <- which(alleles == as.character(A[i,]), arr.ind = TRUE)[,1]
            a[[i]] <- mean(blups[Reduce(intersect, list(a[[i]]))])
          }
          names(a) <- as.character(A)
          for(i in 1:dim(D)[1]) {
            d[[i]] <- which(apply(alleles == as.character(D[i,1]) | alleles == as.character(D[i,2]), 1, sum) == 2)
            d[[i]] <- mean(blups[Reduce(intersect, list(d[[i]]))])
          }
          names(d) <- apply(D, 1, paste, collapse="")
          for(i in 1:dim(T)[1]) {
            t[[i]] <- which(apply(alleles == as.character(T[i,1]) | alleles == as.character(T[i,2]) | alleles == as.character(T[i,3]), 1, sum) == 3)
            t[[i]] <- mean(blups[Reduce(intersect, list(t[[i]]))])
          }
          names(t) <- apply(T, 1, paste, collapse="")
          for(i in 1:dim(F)[1]) {
            f[[i]] <- which(apply(alleles == as.character(F[i,1]) | alleles == as.character(F[i,2]) | alleles == as.character(F[i,3]) | alleles == as.character(F[i,4]), 1, sum) == 4)
            f[[i]] <- mean(blups[Reduce(intersect, list(f[[i]]))])
          }
          names(f) <- apply(F, 1, paste, collapse="")
          for(i in 1:dim(G)[1]) {
            g[[i]] <- which(apply(alleles == as.character(G[i,1]) | alleles == as.character(G[i,2]) | alleles == as.character(G[i,3]) | alleles == as.character(G[i,4]) | alleles == as.character(G[i,5]), 1, sum) == 5)
            g[[i]] <- mean(blups[Reduce(intersect, list(g[[i]]))])
          }
          names(g) <- apply(G, 1, paste, collapse="")
          for(i in 1:dim(S)[1]) {
            s[[i]] <- which(apply(alleles == as.character(S[i,1]) | alleles == as.character(S[i,2]) | alleles == as.character(S[i,3]) | alleles == as.character(S[i,4]) | alleles == as.character(S[i,5]) | alleles == as.character(S[i,6]), 1, sum) == 6)
            s[[i]] <- mean(blups[Reduce(intersect, list(s[[i]]))])
          }
          names(s) <- apply(S, 1, paste, collapse="")

          a <- a[!is.nan(unlist(a))]
          d <- d[!is.nan(unlist(d))]
          t <- t[!is.nan(unlist(t))]
          f <- f[!is.nan(unlist(f))]
          g <- g[!is.nan(unlist(g))]
          s <- s[!is.nan(unlist(s))]

          for(i in 1:length(d)) {
            d[[i]] <- d[[i]] -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(d), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }

          for(i in 1:length(t)) {
            t[[i]] <- t[[i]] -
              sum(unlist(d[which(lapply(lapply(strsplit(names(d), split = ""), function(x) intersect(strsplit(names(t), split = "")[[i]], x)), function(x) length(x) == 2) == TRUE)])) -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(t), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }

          for(i in 1:length(f)) {
            f[[i]] <- f[[i]] -
              sum(unlist(t[which(lapply(lapply(strsplit(names(t), split = ""), function(x) intersect(strsplit(names(f), split = "")[[i]], x)), function(x) length(x) == 3) == TRUE)])) -
              sum(unlist(d[which(lapply(lapply(strsplit(names(d), split = ""), function(x) intersect(strsplit(names(f), split = "")[[i]], x)), function(x) length(x) == 2) == TRUE)])) -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(f), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }

          for(i in 1:length(g)) {
            g[[i]] <- g[[i]] -
              sum(unlist(f[which(lapply(lapply(strsplit(names(f), split = ""), function(x) intersect(strsplit(names(g), split = "")[[i]], x)), function(x) length(x) == 4) == TRUE)])) -
              sum(unlist(t[which(lapply(lapply(strsplit(names(t), split = ""), function(x) intersect(strsplit(names(g), split = "")[[i]], x)), function(x) length(x) == 3) == TRUE)])) -
              sum(unlist(d[which(lapply(lapply(strsplit(names(d), split = ""), function(x) intersect(strsplit(names(g), split = "")[[i]], x)), function(x) length(x) == 2) == TRUE)])) -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(g), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }

          for(i in 1:length(s)) {
            s[[i]] <- s[[i]] -
              sum(unlist(g[which(lapply(lapply(strsplit(names(g), split = ""), function(x) intersect(strsplit(names(s), split = "")[[i]], x)), function(x) length(x) == 5) == TRUE)])) -
              sum(unlist(f[which(lapply(lapply(strsplit(names(f), split = ""), function(x) intersect(strsplit(names(s), split = "")[[i]], x)), function(x) length(x) == 4) == TRUE)])) -
              sum(unlist(t[which(lapply(lapply(strsplit(names(t), split = ""), function(x) intersect(strsplit(names(s), split = "")[[i]], x)), function(x) length(x) == 3) == TRUE)])) -
              sum(unlist(d[which(lapply(lapply(strsplit(names(d), split = ""), function(x) intersect(strsplit(names(s), split = "")[[i]], x)), function(x) length(x) == 2) == TRUE)])) -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(s), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }

          effects[[q]] <- list(unlist(a), unlist(d), unlist(t), unlist(f), unlist(g), unlist(s))

        }
      }

      if(ploidy == 4) {

        for(q in 1:nqtl) {

          blups <- fitted$results[[p]]$fitted$u.hat[[q]]
          alleles <- matrix(unlist(strsplit(rownames(blups), '')), ncol=5, byrow=TRUE)[,-3]

          A <- t(combn(letters[1:8],1))
          D <- t(combn(letters[1:8],2))
          T <- t(combn(letters[1:8],3))
          F <- t(combn(letters[1:8],4))

          a <- vector("list", dim(A)[1])
          d <- vector("list", dim(D)[1])
          t <- vector("list", dim(T)[1])
          f <- vector("list", dim(F)[1])

          for(i in 1:dim(A)[1]) {
            a[[i]] <- which(alleles == as.character(A[i,]), arr.ind = TRUE)[,1]
            a[[i]] <- mean(blups[Reduce(intersect, list(a[[i]]))])
          }
          names(a) <- as.character(A)
          for(i in 1:dim(D)[1]) {
            d[[i]] <- which(apply(alleles == as.character(D[i,1]) | alleles == as.character(D[i,2]), 1, sum) == 2)
            d[[i]] <- mean(blups[Reduce(intersect, list(d[[i]]))])
          }
          names(d) <- apply(D, 1, paste, collapse="")
          for(i in 1:dim(T)[1]) {
            t[[i]] <- which(apply(alleles == as.character(T[i,1]) | alleles == as.character(T[i,2]) | alleles == as.character(T[i,3]), 1, sum) == 3)
            t[[i]] <- mean(blups[Reduce(intersect, list(t[[i]]))])
          }
          names(t) <- apply(T, 1, paste, collapse="")
          for(i in 1:dim(F)[1]) {
            f[[i]] <- which(apply(alleles == as.character(F[i,1]) | alleles == as.character(F[i,2]) | alleles == as.character(F[i,3]) | alleles == as.character(F[i,4]), 1, sum) == 4)
            f[[i]] <- mean(blups[Reduce(intersect, list(f[[i]]))])
          }
          names(f) <- apply(F, 1, paste, collapse="")

          a <- a[!is.nan(unlist(a))]
          d <- d[!is.nan(unlist(d))]
          t <- t[!is.nan(unlist(t))]
          f <- f[!is.nan(unlist(f))]

          for(i in 1:length(d)) {
            d[[i]] <- d[[i]] -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(d), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }

          for(i in 1:length(t)) {
            t[[i]] <- t[[i]] -
              sum(unlist(d[which(lapply(lapply(strsplit(names(d), split = ""), function(x) intersect(strsplit(names(t), split = "")[[i]], x)), function(x) length(x) == 2) == TRUE)])) -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(t), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }

          for(i in 1:length(f)) {
            f[[i]] <- f[[i]] -
              sum(unlist(t[which(lapply(lapply(strsplit(names(t), split = ""), function(x) intersect(strsplit(names(f), split = "")[[i]], x)), function(x) length(x) == 3) == TRUE)])) -
              sum(unlist(d[which(lapply(lapply(strsplit(names(d), split = ""), function(x) intersect(strsplit(names(f), split = "")[[i]], x)), function(x) length(x) == 2) == TRUE)])) -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(f), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }

          effects[[q]] <- list(unlist(a), unlist(d), unlist(t), unlist(f))

        }

      }

    } else {
      effects <- NULL
    }

    results[[p]] <- list(
      pheno.col=fitted$results[[p]]$pheno.col,
      effects=effects)

  }

  structure(list(fitted=deparse(substitute(fitted)),
                 ploidy=ploidy,
                 pheno.col=fitted$pheno.col,
                 results=results),
            class="qtlpoly.effects")

}

#' @rdname qtl_effects
#' @import ggplot2
#' @export

plot.qtlpoly.effects <- function(x, pheno.col = NULL, p1 = "P1", p2 = "P2") {
  if(is.null(pheno.col)) {
    pheno.col <- 1:length(x$results)
  } else {
    pheno.col <- which(x$pheno.col %in% pheno.col)
  }
  for(p in pheno.col) {
    nqtl <- length(x$results[[p]]$effects)
    if(nqtl > 0) {
      for(q in 1:nqtl) {
        if(x$ploidy == 4) {
          data <- unlist(x$results[[p]]$effects[[q]])[1:36]
          data <- data.frame(Estimates=as.numeric(data), Alleles=names(data), Parent=c(rep(p1,4),rep(p2,4),rep(p1,14),rep(p2,14)), Effects=c(rep("Allelic",8),rep("Diallelic",28)))
          data <- data[-c(12:15,18:21,23:30),]
          # plot <- ggplot(data, aes(x = Alleles, y = Estimates, fill = Estimates)) +
          #   geom_bar(stat="identity") +
          #   scale_fill_gradient2(low = "red", high = "blue", guide = FALSE) +
          #   labs(title=names(x$results)[p], subtitle=paste("QTL", q)) +
          #   facet_wrap(Effects ~ Parent, scales="free_x", ncol = 2) +
          #   theme_minimal() +
          #   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), title=element_text(face="bold"))
        }
        if(x$ploidy == 6) {
          data <- unlist(x$results[[p]]$effects[[q]])[-c(18:23,28:33,37:42,45:50,52:63,83:88,92:97,100:105,107:133,137:142,145:150,152:178,181:186,188:214,216:278,299:1763)]
          data <- data.frame(Estimates=as.numeric(data), Alleles=names(data), Parent=c(rep(p1,6),rep(p2,6),rep(p1,15),rep(p2,15),rep(p1,20),rep(p2,20)), Effects=c(rep("Allelic",12),rep("Diallelic",30),rep("Triallelic",40)))
        }
        plot <- ggplot(data, aes(x = Alleles, y = Estimates, fill = Estimates)) +
          geom_bar(stat="identity") +
          scale_fill_gradient2(low = "red", high = "blue", guide = FALSE) +
          labs(title=names(x$results)[p], subtitle=paste("QTL", q, "\n")) +
          facet_wrap(Effects ~ Parent, scales="free_x", ncol = 2) +
          # facet_grid(Effects ~ Parent, scales="free_x", space="free_x") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), title=element_text(face="bold"),
                axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5), strip.text = element_blank(),
                plot.margin = unit(c(0,1,0,0), "lines"))
        print(plot)
        gp <- grid::gpar(fontsize=9, col="grey30")
        grid::grid.text(x=c(.30,.75),y=c(.90,.90), label = c(p1,p2), rot = 0, gp=gp)
        if(x$ploidy == 6) grid::grid.text(x=c(.99,.99,.99),y=c(.74,.49,.21), label = levels(data$Effects), rot = 270, gp=gp) # right
      }
    }
  }
}

