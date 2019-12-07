#' QTLs with respective support interval plots
#'
#' Creates a plot where colored bars represent the support intervals for QTL peaks (black dots).
#'
#' @param data an object of class \code{qtlpoly.data}.
#'
#' @param model an object of class \code{qtlpoly.profile} or \code{qtlpoly.remim}.
#'
#' @param pheno.col a numeric vector with the phenotype column numbers to be plotted; if \code{NULL}, all phenotypes from \code{'data'} will be included.
#'
#' @param main a character string with the main title; if \code{NULL}, no title will be shown.
#'
#' @param drop if \code{TRUE}, phenotypes with no QTL will be dropped; if \code{FALSE} (default), all phenotypes will be shown.
#'
#' @return A \pkg{ggplot2} with QTL bars for each linkage group.
#'
#' @seealso \code{\link[qtlpoly]{read_data}}, \code{\link[qtlpoly]{remim}}, \code{\link[qtlpoly]{profile_qtl}}
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
#'   # plot support intervals
#'   plot_sint(data = data, model = remim.mod)
#'   }
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2019) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{bioRxiv}. \url{doi:}.
#'
#' @export plot_sint
#' @import ggplot2

plot_sint <- function(data, model, pheno.col=NULL, main=NULL, drop=FALSE) {
  trait <- c(rep("LGs", data$nlgs))
  lg <- c(1:data$nlgs)
  lg.lab <- c(1:data$nlgs)
  lower <- c(rep(0, data$nlgs))
  pos <- c(rep(-10, data$nlgs))
  upper <- c(data$lgs.size)

  if(is.null(pheno.col)) pheno.col <- model$pheno.col
  nphen <- length(pheno.col)
  for(p in 1:nphen) {
    t <- which(model$pheno.col == pheno.col[p])
    nqtl <- dim(model$results[[t]]$qtls)[1]
    if(!is.null(nqtl)) {
      trait <- c(trait, rep(names(model$results)[[t]], nqtl), rep(names(model$results)[[t]], data$nlgs))
      lg <- c(lg, model$results[[t]]$qtls[,1], c(1:data$nlgs))
      lg.lab <- c(lg.lab, rep(" ", nqtl), rep(" ", data$nlgs))
      lower <- c(lower, model$results[[t]]$lower[,2], rep(-10,data$nlgs))
      pos <- c(pos, model$results[[t]]$qtls[,2], rep(-10,data$nlgs))
      upper <- c(upper, model$results[[t]]$upper[,2], rep(-10,data$nlgs))
    } else if(!drop) {
      trait <- c(trait, rep(names(model$results)[[t]], data$nlgs))
      lg <- c(lg, c(1:data$nlgs))
      lg.lab <- c(lg.lab, rep(" ", data$nlgs))
      lower <- c(lower, c(rep(-10, data$nlgs)))
      pos <- c(pos, c(rep(-10, data$nlgs)))
      upper <- c(upper, c(rep(-10, data$nlgs)))
    }
  }

  DF <- data.frame(trait=trait, lg=as.integer(lg), lg.lab=lg.lab, lower=lower, pos=pos, upper=upper)
  DF$trait <- factor(DF$trait, levels = unique(DF$trait))
  DF$lg <- factor(DF$lg, levels = unique(DF$lg))
  trait.color <- c("black", hcl(h = seq(15, 375, length = nlevels(DF$trait)), l = 65, c = 100)[1:nlevels(DF$trait)])
  bxp.table <- table(interaction(DF$trait, DF$lg))
  bxp.width <- c(bxp.table[unlist(lapply(interaction(DF$trait, DF$lg), function(x) which(names(bxp.table) %in% x)))])
  bxp.width <- ifelse(bxp.width == 1, 1, bxp.width-.25)

  suppressWarnings({
  ggplot(data = DF) +
  facet_grid(~ lg, scales = "free_x", space = "free_x") +
    geom_boxplot(aes(x=trait, ymin = lower, lower = lower, middle = pos, upper = upper, ymax = upper, color=trait, fill = trait, group = interaction(trait, pos)), position = position_dodge(0), stat="identity", size = 0) +
    geom_crossbar(aes(x=trait, y=pos, ymin = -10, ymax = -10), width = bxp.width, position=position_dodge(0)) +
    scale_fill_manual(values = trait.color) +
    scale_color_manual(values = trait.color) +
    coord_cartesian(ylim = c(max(upper), 5)) +
    labs(title=main, x = "Linkage Group", y = "Position (cM)")+
    scale_x_discrete(position="top") +
    scale_y_reverse(expand = c(0,8)) +
    theme(legend.title=element_blank(), axis.text.x=element_blank(), legend.position="bottom",
          title=element_text(face="bold"), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.ticks.x=element_blank(),
          panel.spacing = unit(0.2, "lines"), strip.text.x = element_text(size = 10), plot.title = element_text(face="bold", hjust = 0.5))
  })


  # ggplot(DF, aes(x = lg, ymin = lower, lower = lower, middle = pos, upper = upper, ymax = upper, group = interaction(trait, lg))) +
  #   geom_boxplot(aes(color=trait, fill = trait), position = position_dodge(1), stat="identity", size = 0) +
  #   geom_crossbar(aes(y=pos, ymin = -10, ymax = -10), width = 1, position=position_dodge(1)) +
  #   scale_fill_manual(values = trait.color) +
  #   scale_color_manual(values = trait.color) +
  #   coord_cartesian(ylim = c(max(upper)+9,0)) +
  #   labs(x = "Linkage Groups", y = "Length (cM)") +
  #   # scale_x_discrete(breaks=lg-1, labels = lg, position = "top") +
  #   # scale_x_continuous(limits=c(0,12), breaks = lg, labels = lg.lab, position = "top") +
  #   # scale_x_continuous(breaks = c(unique(lg)-(p/10)), labels = c(unique(lg)), position = "top") +
  #   scale_x_discrete(position="top") +
  #   scale_y_reverse(expand = c(0,9)) +
  #   # theme_minimal() +
  #   geom_text(aes(x = lg, y = lower, label=lg.lab), position = position_dodge(1), vjust = -0.5) +
  #   theme(legend.title=element_blank(), legend.text = element_text(size=10), axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  #   theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())

}
