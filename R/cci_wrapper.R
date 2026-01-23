# Wrapper functions to use Rcpp acceleration for CCI inference

#' Use Rcpp for fast permutation if available
#' @importFrom Rcpp sourceCpp
.use_cpp_permutation <- function() {
  # Check if Rcpp functions are available
  return(exists("cpp_fast_permutation") && is.function(cpp_fast_permutation))
}

#' Use Rcpp for fast random walk if available
.use_cpp_random_walk <- function() {
  return(exists("cpp_random_walk") && is.function(cpp_random_walk))
}
