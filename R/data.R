#' Simulated autohexaploid map
#'
#' A simulated map containing three homology groups of a hypotetical cross between two autohexaploid individuals.
#'
#' @docType data
#'
#' @usage data(maps6x)
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
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2019) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{Genetics} 215 (3): 579-595. \url{http://doi.org/10.1534/genetics.120.303080}.
#'
#'     Mollinari M, Garcia AAF (2019) Linkage analysis and haplotype phasing in experimental autopolyploid populations with high ploidy level using hidden Markov models, \emph{G3: Genes|Genomes|Genetics} 9 (10): 3297-3314. \url{https://doi.org/10.1534/g3.119.400378}
#'
#' @examples
#' \dontrun{
#' data(maps)
#' library(mappoly)
#' plot(maps)
#' }
"maps6x"

#' Simulated phenotypes
#'
#' A simulated data set of phenotypes for a hipotetical autohexaploid species map.
#'
#' @docType data
#'
#' @usage data(pheno6x)
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
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2019) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{Genetics} 215 (3): 579-595. \url{http://doi.org/10.1534/genetics.120.303080}.
#'
#' @examples
#' \dontrun{
#' data(pheno6x)
#' head(pheno6x)
#' }
"pheno6x"

#' Tetraploid potato genotype probabilities
#'
#' Genotype probabilities for three chromosomes from a tetraploid potato full-sib population (Atlantic x B1829-5).
#'
#' @docType data
#'
#' @usage data(genoprob4x)
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
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2019) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{Genetics} 215 (3): 579-595. \url{http://doi.org/10.1534/genetics.120.303080}.
#'     
#'     Pereira GS, Mollinari M, Schumann MJ, Clough ME, Yencho C (2020) The recombination landscape and multiple QTL mapping in a \emph{Solanum tuberosum} cv. ‘Atlantic’-derived F_1 population. bioRxiv. \url{https://doi.org/10.1101/2020.08.24.265397}.
#'
#' @examples
#' \dontrun{
#' data(genoprob4x)
#' }
"genoprob4x"

#' Tetraploid potato phenotypes
#'
#' A subset of phenotypes from a tetraploid potato full-sib family (Atlantic x B1829-5).
#'
#' @docType data
#'
#' @usage data(pheno4x)
#'
#' @format A data frame of phenotypes with 156 named individuals in rows and three named phenotypes in columns, which are:
#'
#' \describe{
#'   \item{FM07}{Foliage maturity evaluated in 2007.}
#'   \item{FM08}{Foliage maturity evaluated in 2008.}
#'   \item{FM14}{Foliage maturity evaluated in 2014.}
#' }
#'
#' @keywords datasets
#'
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2020) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{Genetics} 215 (3): 579-595. \url{http://doi.org/10.1534/genetics.120.303080}.
#'     
#'     Pereira GS, Mollinari M, Schumann MJ, Clough ME, Yencho C (2020) The recombination landscape and multiple QTL mapping in a \emph{Solanum tuberosum} cv. ‘Atlantic’-derived F_1 population. bioRxiv. \url{https://doi.org/10.1101/2020.08.24.265397}.
#'
#' @examples
#' \dontrun{
#' data(pheno4x)
#' head(pheno4x)
#' }
"pheno4x"

