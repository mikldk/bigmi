#' Danish Y-STR reference database
#'
#' Note, that the duplications at DYS385a/b have been split into two loci, although
#' there is actually an identifiability issue (e.g. 13,14 is not necessarily equal to 13,14).
#' 
#' Also, note that DYS389II contains DYS389I. Sometimes DYS389II' = DYS389II - DYS389I 
#' can be used instead.
#' 
#' @format A matrix with 171 rows and 12 variables:
#' \describe{
#'   \item{DYS19}{DYS19}
#'   \item{DYS389I}{DYS389I}
#'   \item{DYS389II}{DYS389II}
#'   \item{DYS390}{DYS390}
#'   \item{DYS391}{DYS391}
#'   \item{DYS392}{DYS392}
#'   \item{DYS393}{DYS393}
#'   \item{DYS385a}{DYS385a}
#'   \item{DYS385b}{DYS385b}
#'   \item{DYS439}{DYS439}
#'   \item{DYS438}{DYS438}
#'   \item{DYS437}{DYS437}
#' }
#' @source \url{https://doi.org/10.1016/j.forsciint.2004.12.019}
#' @references C Hallenberg et al. "Y-chromosome STR haplotypes in Danes", Forensic Science International, 155(2-3), 2005.
"danedb"
