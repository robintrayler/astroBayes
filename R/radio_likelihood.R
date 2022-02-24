radio_likelihood <- function(anchored_model,
                                    position_grid,
                                    age,
                                    age_sd,
                                    id,
                                    position) {

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

