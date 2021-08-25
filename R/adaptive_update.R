#' This function is a simple implementation of the adaptive MCMC algorithm of
#' Haario et al. (2001)
#' @name adaptive_update
#' @param chain the current markov chain to update
#' @param i current chain iteration
#' @param start_index when should the proposal start adapting?
#' @param initial_Cd initial guess at the proposal variance
#'
#' @import "tidyverse"
#' @import "dplyr"
#' @import "tibble"
#' @importFrom magrittr "%>%"
#'
#' @return a list with a whole bunch of stuff in it.
#' @md
#' @export

adaptive_update <- function(chain,
                            i,
                            start_index = 1000,
                            initial_Cd = 0.01,
                            distribution = c('gamma', 'gaussian')){

  # set some global parameters ------------------------------------------------
  esp <- 1e-5 # keep things from going to zero
  S_d = 2.4^2 # scaling factor
  if(i > start_index) {
    # calculate the proposal variance -----------------------------------------
    # per Haario et al., (2001)
    C_d <- var(chain[1:(i-1)]) * S_d + S_d * esp %*% diag(1)

    x_proposed <- numeric()
    # propose a new value -----------------------------------------------------
    if(distribution == 'gamma') {
      alpha <- chain[i - 1]^2 / C_d
      beta  <- chain[i - 1] / C_d
      x_proposed <- rgamma(n = 1,
                           shape = alpha,
                           rate = beta)

    } else if (distribution == 'gaussian') {
      mu    <- chain[i - 1]
      sd    <- sqrt(C_d)
      x_proposed <- rnorm(n = 1,
                          mean = mu,
                          sd = sd)

    }
  } else {

    if(distribution == 'gamma') {
      alpha <- chain[i - 1]^2 / initial_Cd
      beta  <- chain[i - 1] / initial_Cd
      x_proposed <- rgamma(n = 1,
                           shape = alpha,
                           rate = beta)
    } else if (distribution == 'gaussian') {
      mu    <- chain[i - 1]
      sd    <- sqrt(initial_Cd)
      x_proposed <- rnorm(n = 1,
                          mean = mu,
                          sd = sd)
    }
  }
  return(x_proposed)
}
