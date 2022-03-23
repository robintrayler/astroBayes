radio_likelihood <- function(anchored_model,
                             position_grid,
                             age,
                             age_sd,
                             position) {
  # INPUTS
  # anchored_model = age model anchored in time
  # position_grid  = grid for interpolation
  # position = position to calculate likelihood of
  # age, age_sd = normal distribution parameters to use for LL calculations
  # create an interpolation function ------------------------------------------
  f <- approxfun(x = position_grid,
                 y = anchored_model)

  # predict age at dated positions and calculate log-likelihood ---------------
  LL <- f(position) %>%
    dnorm(mean = age,
          sd = age_sd,
          log = TRUE) %>%
    sum()

  # return the log likelihood
  return(LL)
}

