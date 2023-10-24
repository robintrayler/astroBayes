#' This function uses the output of an \code{\link{astro_bayes_model}} to summarize
#'a variety of model parameters.
#'s
#' @name summary.astroBayesModel
#' @param age_model the output of \code{\link{astro_bayes_model}}
#' @param type one of:
#' * `sed_rate` calculate the posterior 95% credible interval of sedimentation rate for each layer
#' * `dates`    calculate the posterior 95% credible interval of the dated horizons
#' * `hiatus`   calculate the posterior 95% credible interval of hiatus duration(s)
#' * `anchor`   calculate the posterior 95% credible interval of the anchor point age
#'
#' @import "tidyverse"
#' @import "dplyr"
#' @import "rlang"
#' @md
#' @export
#
summary.astroBayesModel <- function(age_model, type) {
  # choose type of summary to calculate
  switch(type,
         'sed_rate' = summarize_sed_rate(age_model),
         'dates' = summarize_dates(age_model),
         'hiatus' = summarize_hiatus(age_model),
         'anchor' = summarize_anchor(age_model))
}

# sedimentation rate ----------------------------------------------------------
summarize_sed_rate <- function(age_model) {
  # preallocate some storage
  sed_rate_CI <- data.frame(
    min_position = vector(length = ncol(age_model$sed_rate)),
    max_position = vector(length = ncol(age_model$sed_rate)),
    CI_2.5 = vector(length = ncol(age_model$sed_rate)),
    median = vector(length = ncol(age_model$sed_rate)),
    CI_97.5 = vector(length = ncol(age_model$sed_rate)))

  # calculate credible interval for each sed rate
  for(i in 1:ncol(age_model$sed_rate)){
    sed_rate_CI$min_position[i] <- age_model$layer_boundaries$position[i]
    sed_rate_CI$max_position[i] <- age_model$layer_boundaries$position[i+1]
    sed_rate_CI[i, 3:5] <- quantile(age_model$sed_rate[, i],
                                    prob = c(0.25, 0.5, 0.975))
  }
  return(sed_rate_CI)
}

# posterior dates -------------------------------------------------------------

summarize_dates <- function(age_model) {
  # calculate 95% CI of the dates
  dates_CI <- predict(age_model,
                      new_positions = age_model$geochron_data)
  # return it
  return(dates_CI$CI)
}

# hiatus position(s) ----------------------------------------------------------
summarize_hiatus <- function(age_model) {
  # how many hiatuses
  n_hiatus <- sum(age_model$layer_boundaries$hiatus_boundary)

  # preallocate storage
  hiatus_CI <- data.frame(
    position = vector(length = n_hiatus),
    CI_2.5   = vector(length = n_hiatus),
    median   = vector(length = n_hiatus),
    CI_97.5  = vector(length = n_hiatus))

  if(n_hiatus < 1) {
    stop('model does not include any hiatus(es)')
  }

  if(n_hiatus > 0) {

    # loop through the hiatuses
    for(i in 1:n_hiatus) {
      hiatus_CI$position[i] <- age_model$layer_boundaries$position[age_model$layer_boundaries$hiatus_boundary][i]

      hiatus_CI[i, 2:4] <- quantile(age_model$hiatus_durations[, i],
                                    prob = c(0.025, 0.5, 0.975))
    }
    return(hiatus_CI)
  }
}


# posterior anchor point ------------------------------------------------------
summarize_anchor <- function(age_model) {
  # calculate 95% CI of the anchor_point

  anchor_CI <- age_model$anchor_point |>
    quantile(prob = c(0.25, 0.5, 0.975)) |>
    t() |>
    as.data.frame() |>
    set_names(c('CI_2.5', 'median', 'CI_97.5')) |>
    add_column(id = 'anchor_point')

  # return it
  return(anchor_CI)
}
