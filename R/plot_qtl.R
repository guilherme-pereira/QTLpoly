#' QTL heritability and significance plot
#'
#' Creates a plot where dot sizes and colors represent the QTLs heritabilities and their \emph{p}-values, respectively.
#'
#' @param data an object of class \code{qtlpoly.data}.
#'
#' @param model an object of class \code{qtlpoly.profile} or \code{qtlpoly.remim}.
#'
#' @param fitted an object of class \code{qtlpoly.fitted}.
#'
#' @param pheno.col the desired phenotype column numbers to be plotted. The order here specifies the order of plotting (from top to bottom.)
#'
#' @param main plot title; if \code{NULL} (the default), no title is shown.
#'
#' @param drop.pheno if \code{FALSE}, shows the names of all traits from \code{pheno.col}, even of those with no QTLs; if \code{TRUE} (the default), shows only the traits with QTL(s).
#'
#' @param drop.lgs if \code{FALSE}, shows all linkage groups, even those with no QTL; if \code{TRUE} (the default), shows only the linkage groups with QTL(s).
#'
#' @return A \pkg{ggplot2} with dots representing the QTLs.
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
#'   # plot qtls
#'   plot_qtl(data = data, model = remim.mod, fitted = fitted.mod)
#'   }
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2019) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{bioRxiv}. \url{doi.org/10.1101/622951}.
#'
#' @export plot_qtl
#' @import ggplot2

plot_qtl <- function(data=data, model=model, fitted=fitted, pheno.col=NULL, main=NULL, drop.pheno=TRUE, drop.lgs=TRUE) {
  trt <- c()
  lgs <- c()
  pos <- c()
  pvl <- c()
  h2q <- c()
  
  if(is.null(pheno.col)) pheno.col <- model$pheno.col
  nphen <- length(pheno.col)
  for(p in 1:nphen) {
    t <- which(model$pheno.col == pheno.col[p])
    nqtl <- dim(model$results[[t]]$qtls)[1]
    pvls <- model$results[[t]]$qtls[,6]
    if(any(pvls == "<2.22e-16")) pvls[which(pvls == "<2.22e-16")] <- 2.22e-16
    if(any(pvls == "<1.00e-16")) pvls[which(pvls == "<1.00e-16")] <- 2.22e-16
    if(!is.null(nqtl)) {
      trt <- c(trt, rep(names(model$results)[[t]], nqtl))
      lgs <- c(lgs, unlist(model$results[[t]]$qtls[1:nqtl,1]))
      pos <- c(pos, unlist(model$results[[t]]$qtls[1:nqtl,2]))
      pvl <- c(pvl, pvls)
      h2q <- c(h2q, unlist(fitted$results[[t]]$qtls[1:nqtl,7]))
    } else if(!drop.pheno) {
      trt <- c(trt, names(model$results)[[t]])
      lgs <- c(lgs, NA)
      pos <- c(pos, 0)
      pvl <- c(pvl, NA)
      h2q <- c(h2q, NA)
    }
  }
  lgs[which(is.na(lgs))] <- lgs[which(!is.na(lgs))][1]
  if(!drop.lgs) {
    trt <- c(trt, rep(trt[1], 2*data$nlgs))
    lgs <- c(lgs, c(1:data$nlgs), c(1:data$nlgs))
    pos <- c(pos, rep(0, data$nlgs), data$lgs.size)
    pvl <- c(pvl, rep(NA, 2*data$nlgs))
    h2q <- c(h2q, rep(NA, 2*data$nlgs))
  } else {
    nlgs <- length(unique(lgs))
    trt <- c(trt, rep(trt[1], 2*nlgs))
    lgs <- c(lgs, unique(lgs), unique(lgs))
    pos <- c(pos, rep(0, nlgs), data$lgs.size[unique(lgs)])
    pvl <- c(pvl, rep(NA, 2*nlgs))
    h2q <- c(h2q, rep(NA, 2*nlgs))
  }
  
  DF <- data.frame(Trait=trt, LG=lgs, Pos=pos, Pval=as.numeric(pvl), H2=h2q)
  DF$Trait = with(DF, factor(Trait, levels = rev(unique(Trait))))
  
  ggplot(data = DF, aes(x = Pos)) +
    facet_grid(~ LG, scales = "free_x", space = "free_x") +
    geom_point(aes(x = Pos, y = Trait, color = Pval, size = H2), na.rm=TRUE) + 
    scale_color_gradient(trans = "log10") +
    scale_radius(range=c(1, 5)) +
    labs(title=main, subtitle = "Linkage group", y = "Trait", x = "Position (cM)", col=expression(paste(italic(P), "-value", sep="")), size=expression(italic(h)[QTL]^{2}))+
    scale_x_continuous(breaks = seq(0, max(data$cum.size), 100), expand = expand_scale(add = 50)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1), panel.spacing = unit(0, "lines"), plot.subtitle = element_text(hjust = 0.5))          
  # suppressWarnings({
  # ggplot(DF, aes(x = Pos, y = Trait, colour = Pval, size = H2)) +
  #   geom_point() +
  #   # {if(!is.null(scale)) scale_radius(range=c(1, scale))} +
  #   scale_radius(range=c(1, 5)) +
  #   scale_color_gradient(trans = "log10") +
  #   # annotate("rect", xmin = data$cum.size[1], xmax = data$cum.size[2], ymin = 3.5, ymax = 3.8, colour=NA, fill="black", alpha=.2) +
  #   labs(title=main, x = "Linkage Group", y = "Trait", col=expression(paste(italic(p), "-value", sep="")), size=expression(italic(h)[QTL]^{2})) +
  #   scale_x_continuous(breaks = c(data$cum.size[c(1:data$nlgs)+1]-(data$lgs.size/2)), labels = c(1:data$nlgs),
  #                      limits = c(0,max(data$cum.size)), minor_breaks = c(data$cum.size), position = "top", expand = c(0, 0)) +
  #   scale_y_discrete(breaks = DF$Trait) +
  #   # theme_minimal() +
  #   # guides(colour = guide_legend(order=1)) + # make sure it's working
  #   # theme(legend.text = element_text(size=10), axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  #   theme_bw() +
  #   theme(axis.title=element_text(face="bold"), plot.title = element_text(face="bold", hjust = 0.5), axis.ticks.x=element_blank(),
  #         axis.text.x = element_text(size=10), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_line(colour = "black"), panel.grid.major.y = element_line(linetype = "dashed"))
  # })
  
}
