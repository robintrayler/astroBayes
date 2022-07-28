#' plots the age model with error bars for new predictions
#' @param age_predictions the output of \code{predict.astroBayesModel}
#'
#'  @import "tidyverse"
#' @import "dplyr"
#' @import "astrochron"
#' @import "cowplot"
#' @import "tibble"
#' @importFrom magrittr "%>%"
#' @md
#' @export
#'

plot.astroBayesPrediction <- function(age_predictions) {
  # plot the original age model -----------------------------------------------
  p <- age_predictions$age_model %>% plot(type = 'age_depth')

  # add the error bars for the new points -------------------------------------
  p <- p + geom_point(data = age_predictions$CI,
                      mapping = aes(x = median,
                                    y = position,
                                    group = id),
                      color = 'tomato',
                      inherit.aes = FALSE) +
    geom_errorbar(data = age_predictions$CI,
                  mapping = aes(
                    xmin = CI_2.5,
                    xmax = CI_97.5,
                    y = position),
                  width = 0,
                  color = 'tomato',
                  inherit.aes = FALSE) +
    geom_errorbar(data = age_predictions$CI,
                  mapping = aes(
                    ymin = position - thickness/2,
                    ymax = position + thickness/2,
                    x = median),
                  width = 0,
                  color = 'tomato',
                  inherit.aes = FALSE)
  p
}



