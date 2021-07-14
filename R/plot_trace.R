#' plots the MCMC trace for each accumulation rate segment
#' @param age_model the output of \code{anchored_age_model}
#'
#' @import "ggplot2"
#' @import "cowplot"
#'
#' @return a plot
#' @md
#' @export
#'

plot_trace <- function(age_model) {
  plots <- list()
  for(k in 1:ncol(age_model$sed_rate)){
    dat <- data.frame(x = 1:age_model$iterations,
                      y = age_model$sed_rate[, k])

    plots[[k]] <- ggplot(data = dat,
                         mapping = aes(x = x,
                                       y = y)) +
      geom_line() +
      xlab('iteration') +
      ylab('sed rate') +
      theme_bw()
  }
  cowplot::plot_grid(plotlist = plots,
            labels = 'AUTO',
            nrow = length(plots))
}
