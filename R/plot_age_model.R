#' plots the age model
#' @param age_model the output of \code{anchored_age_model}
#'
#' @import "ggplot2"
#' @import "cowplot"
#'
#' @return a list with a whole bunch of stuff in it.
#' @md
#' @export
#'

plot_age_model <- function(age_model){
  ggplot(data = age_model$geochron_data,
         mapping = aes(x = position,
                       y = age,
                       color = id)) +
    geom_point(size = 2) +
    geom_errorbar(mapping = aes(ymin = age - age_sd * 2,
                                ymax = age + age_sd * 2),
                  width = 0.75) +
    xlab('Depth (m)') +
    ylab('Age (Ma)') +
    scale_x_reverse() +
    coord_flip() +
    geom_ribbon(data = age_model$CI,
                mapping = aes(ymin = CI_2.5,
                              ymax = CI_97.5,
                              x = grid),
                inherit.aes = FALSE,
                alpha = 0.25) +
    geom_line(data = age_model$CI,
              mapping = aes(x = grid,
                            y = median),
              inherit.aes = FALSE) +
    theme_bw() +
    theme(legend.position = 'top')
}
