#' This function uses the output of an \code{\link{astro_bayes_model}} to predict the
#' age of new stratigraphic positions specified in \code{new_positions}
#' @name predict.astroBayesModel
#' @param age_model the output of \code{\link{astro_bayes_model}}
#' @param new_positions a data frame containing the name, position and positional
#' uncertainty (thickness) of stratigraphic points to predict the age of.
#' column headers *must* be named exactly as follows:
#' * id: name
#' * position: stratigraphic position. The position must fall within the original evaluation
#' range of \code{age_model}
#' * thickness: positional uncertainty is treated as a uniform distributions where the total
#' range equals position ± thickness ÷ 2
#'
#' @import "tidyverse"
#' @import "dplyr"
#' @import "astrochron"
#' @import "tibble"
#' @import "rlang"
#' @importFrom magrittr "%>%"
#' @md
#' @export
#

predict.astroBayesModel <- function(age_model, new_positions) {
  # preallocate storage -------------------------------------------------------
  predict_store <- matrix(nrow = age_model$iterations,
                          ncol = nrow(new_positions))

  pb <- progress::progress_bar$new(total = age_model$iterations,
                                   format = '[:bar] :percent eta: :eta')
  for(i in seq_along(1:age_model$iterations)) {
    pb$tick()
    # form interpolation function ---------------------------------------------
    f <- approxfun(x = age_model$CI$position,
                   y = age_model$model_iterations[, i])

    # randomize the positions
    predict_store[i, ] <- runif(n = nrow(new_positions),
                                min = new_positions$position -
                                  new_positions$thickness/2,
                                max = new_positions$position +
                                  new_positions$thickness/2) %>%
      f()
  }
  # organize storage ----------------------------------------------------------
  posterior_sample <- predict_store %>%
    as.data.frame() %>%
    rlang::set_names(nm = new_positions$id)

  # calculate credible interval -----------------------------------------------
  credible_interval <- apply(X = posterior_sample[age_model$burn:age_model$iterations, ],
                             MARGIN = 2,
                             FUN = quantile,
                             prob = c(0.025, 0.5, 0.975)) %>%
    t() %>%
    as.data.frame() %>%
    rlang::set_names(nm = c('CI_2.5', 'median', 'CI_97.5')) %>%
    add_column(position = new_positions$position,
               thickness = new_positions$thickness,
               id = new_positions$id)


  output <- list(CI = credible_interval,
                 posterior = posterior_sample,
                 age_model = age_model)
  class(output) <- "astroBayesPrediction"

  return(output)
}
