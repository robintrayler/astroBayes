radio_likelihood <- function(segment_edges,
                             anchored_model,
                             age,
                             age_sd,
                             id,
                             position) {
  # form approximation function
  f <- approxfun(x = segment_edges,
                 y = anchored_model)
  # predict age at dated positions
  LL <- f(position) %>%
    # calculate probability
    dnorm(age, age_sd) %>%
    log() %>%
    sum()

  return(LL)
}
