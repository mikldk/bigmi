#' Emperical mutual information
#' 
#' Calculate the emperical mutual information between two variables, i.e.  
#' assuming that they follow categorical distributions and use the observed counts.
#' 
#' Note, that this is a symmetric measure, i.e. the ordering of the columns in `x`
#' does not matter.
#' 
#' @param x Integer matrix with 2 columns and no `NA`; one row per observation; integer matrix is required to
#' underline that it is categorical mutual information.
#' 
#' @return The mutual information between columns in `x` assuming 
#' that they follow a categorical distribution.
#' 
MI_categorical_two <- function(x) {
  stopifnot(!anyNA(x))
  stopifnot(is.matrix(x))
  stopifnot(is.integer(x))
  stopifnot(ncol(x) == 2L)
  
  return(MI_categorical_worker_two(x))
}

#' Emperical mutual information
#' 
#' Calculate the emperical mutual information between all variables, i.e.  
#' assuming that they follow categorical distributions and use the observed counts.
#' 
#' Note, that this is a symmetric measure, i.e. the ordering of the columns in `x`
#' does not matter.
#' 
#' @param x Integer matrix with >= 2 columns and no `NA`; one row per observation; integer matrix is required to
#' underline that it is categorical mutual information.
#' 
#' @return The mutual information between columns in `x` assuming 
#' that they follow a categorical distribution. 
#' Note that numbers off-diagonal are mutual information and the numbers in the diagonal are entropies.
#' 
MI_categorical_all_pairwise <- function(x) {
  stopifnot(!anyNA(x))
  stopifnot(is.matrix(x))
  stopifnot(is.integer(x))
  stopifnot(ncol(x) >= 2L)
  
  return(MI_categorical_worker_all(x))
}
