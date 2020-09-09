#' Logarithm of \emph{P}-value (LOP) profile plots
#'
#' Plots profiled logarithm of score-based \emph{P}-values (LOP) from individual or combined traits.
#'
#' @param data an object of class \code{qtlpoly.data}.
#'
#' @param model an object of class \code{qtlpoly.profile} or \code{qtlpoly.remim}.
#'
#' @param pheno.col a numeric vector with the phenotype column numbers to be plotted; if \code{NULL}, all phenotypes from \code{'data'} will be included.
#'
#' @param sup.int if \code{TRUE}, support interval are shown as shaded areas; if \code{FALSE} (default), no support interval is show.
#' 
#' @param main a character string with the main title; if \code{NULL}, no title is shown.
#'
#' @param legend legend position (either "bottom", "top", "left" or "right"); if \code{NULL}, no legend is shown.
#'
#' @param ylim a numeric value pair supplying the limits of y-axis, e.g. c(0,10); if \code{NULL} (default), limits will be provided automatically.
#'
#' @param grid if \code{TRUE}, profiles will be organized in rows (one per trait); if \code{FALSE} (default), profiles will appear superimposed. Only effective when plotting profiles from more than one trait.
#'
#' @return A \pkg{ggplot2} with the LOP profiles for each trait.
#'
#' @seealso \code{\link[qtlpoly]{profile_qtl}},  \code{\link[qtlpoly]{remim}}
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
#'   # plot profiles
#'   for (p in remim.mod$pheno.col) { 
#'     plot_profile(data = data, model = remim.mod, pheno.col = p, ylim = c(0, 10))
#'   } # separate plots
#'     
#'   plot_profile(data = data, model = remim.mod, grid = FALSE) # combined plots
#'   }
#'
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2019) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{bioRxiv}. \url{doi.org/10.1101/622951}.
#'
#' @export plot_profile
#' @import ggplot2

plot_profile <- function(data = data, model = model, pheno.col = NULL, sup.int = FALSE, main = NULL, legend="bottom", ylim = NULL, grid = FALSE) {
  
  lines <- points <- thre <- map <- data.frame()
  y.dat <- trait.names <- c()
  count <- 0
  if(is.null(pheno.col)) pheno.col <- model$pheno.col
  nphe <- length(pheno.col)
  LGS <- c(); for(c in 1:length(data$lgs)) LGS <- c(LGS, rep(c, length(data$lgs[[c]])))
  POS <- unlist(data$lgs)
  for(p in 1:nphe) { #lines
    t <- which(model$pheno.col == pheno.col[p])
    TRT <- rep(names(model$results)[t], length(LGS))
    if(any(class(model) == "qtlpoly.feim")) SIG <- model$results[[t]][[3]] else SIG <- -log10(as.numeric(model$results[[t]][[3]]))
    lines <- rbind(lines, data.frame(TRT=as.factor(TRT), LGS=LGS, POS=POS, SIG=SIG))
  }
  for(p in 1:nphe) { #points
    t <- which(model$pheno.col == pheno.col[p])
    trait.names <- c(trait.names, names(model$results)[t])
    if(!is.null(model$results[[t]]$qtls)) {
      nqtls <- dim(model$results[[t]]$qtls)[1]
      TRT <- rep(names(model$results)[t], nqtls)
      LGS <- model$results[[t]]$qtls[,"LG"]
      POS <- model$results[[t]]$qtls[,"Pos"]
      INF <- model$results[[t]]$lower[,"Pos_lower"]
      SUP <- model$results[[t]]$upper[,"Pos_upper"]
      points <- rbind(points, data.frame(TRT=TRT, LGS=LGS, POS=POS, INF=INF, SUP=SUP))
      count <- count+1
      y.dat <- c(y.dat, rep((-0.3*count), nqtls))
    }
  }
  points$TRT <- factor(points$TRT, levels=trait.names)
  if(any(class(model) == "qtlpoly.feim")) {
    for(p in 1:nphe) { #threshold
      t <- which(model$pheno.col == pheno.col[p])
      LGS <- c(1:length(data$lgs))
      TRT <- rep(names(model$results)[t], length(LGS))
      SIG <- rep(model$sig.lod[t], length(LGS))
      thre <- rbind(thre, data.frame(TRT=as.factor(TRT), LGS=LGS, SIG=SIG))
      y.lab <- "LOD"
    }
  } else {
    y.lab <- "LOP"
  }
  if(is.null(y.dat)) y.dat <- ylim[1]
  if(grid) y.dat <- -max(lines$SIG[is.finite(lines$SIG)])/10
  
  for(c in 1:data$nlgs) {
    LGS <- rep(c, length(data$lgs.all[[c]]))
    map <- rbind(map, data.frame(LGS=LGS, POS=data$lgs.all[[c]]))
  }
  if(max(data$lgs.size) > 200) cutx <- 150 else cutx <- 100
  if(max(lines$SIG[is.finite(lines$SIG)])) cuty <- 2 else cuty <- 4
  
  pl <- ggplot(data = lines, aes(x = POS)) +
    {if(grid) facet_grid(TRT ~ LGS, scales = "free", space = "free", shrink = TRUE) else facet_grid(~ LGS, scales = "free_x", space = "free_x")} +
    {if(nrow(points) > 0 & sup.int) geom_rect(data=points, aes(xmin = INF, xmax = SUP, ymin = -Inf, ymax = Inf, fill = TRT), alpha = 0.2)} +
    geom_line(data=lines, aes(y = SIG, color = TRT), size=1, alpha=0.8, lineend = "round") +
    geom_point(data=map, aes(y=0, x=POS), shape="|", alpha=200/dim(map)[1]) +
    {if(nrow(points) > 0) geom_point(data=points, aes(y = y.dat, color = TRT), shape = 2, size = 2, stroke = 1, alpha = 0.8)} +
    scale_x_continuous(breaks=seq(0,max(data$lgs.size),cutx), expand = expansion(add=50)) +
    # {if(!is.null(ylim)) scale_y_continuous(limits = ylim) else scale_y_continuous(limits = c(min(y.dat), 10))} +
    {if(!is.null(ylim)) scale_y_continuous(limits = c(min(y.dat), ylim[2]))} +
    {if(grid) scale_y_continuous(breaks = seq(0, max(lines$SIG[is.finite(lines$SIG)]), cuty), expand = expansion(add=c(0.5, 0.5)))} +
    # {if(grid) scale_y_continuous(breaks = seq(0, max(lines$SIG), cuty), expand = c(0.05,0.15))} +
    # {if((!is.null(ylim) & length(pheno.col) == 1) | (!is.null(ylim) & !grid)) scale_y_continuous(limits = ylim)} +
    {if(nrow(thre) > 0) geom_hline(data=thre, aes(yintercept=SIG, color=TRT), linetype="dashed", size=.5, alpha=0.8)} + #threshold
    guides(color = guide_legend("Trait"), fill = guide_legend("Trait"), shape = guide_legend("Trait")) + 
    labs(title=main, y = y.lab, x = "Position (cM)", subtitle="Linkage group") +
    # {if(!is.null(main)) labs(title=main)} +
    theme_minimal() +
    theme(legend.position=legend, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          panel.spacing.x = unit(0.01, "lines"), panel.spacing.y = unit(0.05, "lines"), strip.text.y = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    {if(is.null(legend)) guides(color = FALSE)}
  print(pl)
}
