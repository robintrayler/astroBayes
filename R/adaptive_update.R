#' This function is a simple implementation of the adaptive MCMC algorithm of
#' Haario et al. (2001)
#' @name adaptive_update
#' @param chain the current markov chain to update
#' @param i current chain iteration
#' @param start_index when should the proposal start adapting?
#' @param initial_Cd initial guess at the proposal variance
#' @param lower lower truncation for the random normal distribution. Defaults to `-1E-10`
#' @param upper upper truncation for the random normal distribution. Defaults to `1E-10`
#'
#' @import "tidyverse"
#' @import "dplyr"
#' @import "tibble"
#' @importFrom magrittr "%>%"
#'
#' @return a a proposed value from `chain`
#' @md
#' @export

adaptive_update <- function(chain,
                            i,
                            start_index = 1000,
                            initial_Cd = 0.01,
                            lower = -1e10,
                            upper = 1e10) {

  # set some global parameters ------------------------------------------------
  esp <- 1e-5 # keep things from going to zero
  S_d = 2.4^2 # scaling factor
  C_d <- var(chain[1:(i-1)]) * S_d + S_d * esp %*% diag(1)

  if(i >= start_index) {
    # calculate the proposal variance -----------------------------------------
    # per Haario et al., (2001)

    x_proposed <- truncated_random_normal(
      mean = chain[i - 1],
      sd = sqrt(C_d),
      low = lower,
      high = upper)

  } else if (i < start_index) {
    x_proposed <- truncated_random_normal(
      mean = chain[i - 1],
      sd = sqrt(initial_Cd),
      low = lower,
      high = upper)
  }
  return(x_proposed)
}
