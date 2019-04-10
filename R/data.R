#' Simulated autohexaploid map
#'
#' A simulated map containing three homology groups of a hypotetical cross between two autohexaploid individuals.
#'
#' @docType data
#'
#' @usage data(maps)
#'
#' @format An object of class \code{"mappoly.map"} from the package \pkg{mappoly}, which is a list of three linkage groups (LGs):
#'
#' \describe{
#'   \item{LG 1}{538 markers distributed along 112.2 cM}
#'   \item{LG 2}{329 markers distributed along 54.6 cM}
#'   \item{LG 3}{443 markers distributed along 98.2 cM}
#' }
#'
#' @keywords datasets
#'
#' @seealso \code{\link[mappoly]{hexafake}}, \code{\link[qtlpoly]{pheno}}
#'
#' @author Marcelo Mollinari, \email{mmollin@@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2019) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{bioRxiv}. \url{doi:}.
#'
#'     Mollinari M, Garcia AAF (2019) Linkage analysis and haplotype phasing in experimental autopolyploid populations with high ploidy level using hidden Markov models, \emph{bioRxiv}. \url{https://doi.org/10.1101/415232}
#'
#' @examples
#' \dontrun{
#' data(maps)
#' library(mappoly)
#' plot(maps)
#' }
"maps"

#' Simulated phenotypes
#'
#' A simulated data set of phenotypes for a hipotetical autohexaploid species map.
#'
#' @docType data
#'
#' @usage data(pheno)
#'
#' @format A data frame of phenotypes with 300 named individuals in rows and three named phenotypes in columns, which are:
#'
#' \describe{
#'   \item{T32}{3 QTLs, with heritabilities of 0.20 (LG 1 at 32.03 cM), 0.15 (LG 1 at 95.02 cM) and 0.30 (LG 2 at 40.01 cM).}
#'   \item{T17}{1 QTL, with heritability of 0.15 (LG 3 at 34.51 cM).}
#'   \item{T45}{no QTLs.}
#' }
#'
#' @keywords datasets
#'
#' @seealso \code{\link[qtlpoly]{simulate_qtl}}, \code{\link[qtlpoly]{pheno}}
#'
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2019) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{bioRxiv}. \url{doi:}.
#'
#' @examples
#' \dontrun{
#' data(pheno)
#' head(pheno)
#' }
"pheno"
