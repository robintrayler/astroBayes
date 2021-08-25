#' plots the best fit periodogram for each accumulation rate segment
#' @param age_model the output of \code{anchored_age_model}
#'
#' @import "ggplot2"
#' @import "cowplot"
#'
#' @return a plot
#' @md
#' @export
#'

plot_pgram <- function(age_model) {
  plots <- list()
  for(k in 1:ncol(age_model$sed_rate)) {
    pgram <- age_model$cyclostrat_data %>%
      filter(position > age_model$segment_edges$position[k] &
               position < age_model$segment_edges$position[k + 1]) %>%
      periodogram(output = 1,
                  verbose = FALSE,
                  genplot = FALSE,
                  background = 1) %>%
      mutate(probability = Power / AR1_Fit,
             probability = probability / sum(probability),
             time_freq   = Frequency * mean(age_model$sed_rate[age_model$burn:age_model$iterations, k],
                                            na.rm = TRUE))

    plots[[k]] <- pgram %>%
      ggplot(mapping = aes(x = time_freq,
                           y = probability)) +
      geom_line() +
      geom_vline(data = age_model$tuning_frequency,
                 mapping = aes(xintercept = frequency,
                               color = orbital_cycle),
                 linetype = 'dashed') +
      ggtitle(label = paste(age_model$segment_edges$position[k],
                            '-',
                            age_model$segment_edges$position[k + 1],
                            'meters',
                            '; mean sed rate = ',
                            round(mean(age_model$sed_rate[, k],
                                       na.rm = TRUE), 2),
                            'm/Ma')) +
      theme(legend.position = 'none',
            axis.text.y = element_blank()) +
      xlab('frequency (m/Ma)') +
      theme_bw() +
      theme(legend.position = 'none')
  }
  cowplot::plot_grid(plotlist = plots,
            labels = 'AUTO',
            nrow = length(plots))
}
