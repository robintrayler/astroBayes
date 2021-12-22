hiatus_radio_likelihood <- function(anchored_model,
                                    position_grid,
                                    age,
                                    age_sd,
                                    id,
                                    position) {

  f <- approxfun(x = position_grid,
                 y = anchored_model)

  # predict age at dated positions
  LL <- f(position) %>%
    # calculate probability
    dnorm(mean = age, sd = age_sd) %>%
    log() %>%
    sum()

  return(LL)
}

