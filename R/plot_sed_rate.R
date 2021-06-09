#' plots the sedimentation rate density for each accumulation rate segment
#' @param age_model the output of \code{anchored_age_model}
#'
#' @import "ggplot2"
#' @import "cowplot"
#'
#' @return a plot
#' @md
#' @export
#'

plot_sed_rate <- function(age_model) {
  plots <- list()

  for(k in 1:ncol(age_model$sed_rate)){
    dat <- data.frame(x = 1:age_model$iterations,
                      y = age_model$sed_rate[, k])

    plots[[k]] <- dat %>%
      ggplot(mapping = aes(x = y)) +
      geom_density(fill = 'grey',
                   color = NA,
                   adjust = 2) +
      xlim(age_model$sed_prior_range) +
      xlab('sed rate (m/Ma)') +
      theme_bw()
  }
  plot_grid(plotlist = plots,
            labels = 'AUTO',
            nrow = length(plots))

}
