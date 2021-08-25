#' plots the age model
#' @param age_model the output of \code{astro_bayes_model}
#'
#' @param type type of plot to make. Possible options are
#'
#' * age_depth: plots the geochrologic data where age is on the x axis and depth
#'   is on the y axis.
#' * sed_rate: plots the posterior density of sedimentation rate for each segment
#' * trace: plots the trace plot for each segment
#' * periodogram: plots the  periodogram for each segment scaled to the median
#'   sedimentation rate for that segment
#'
#' @import "ggplot2"
#' @import "cowplot"
#' @import "viridis"
#'
#' @return a list with a whole bunch of stuff in it.
#' @md
#' @export
#'

plot.astroBayesModel <- function(age_model,
                                 type = c('age_depth',
                                          'sed_rate',
                                          'trace',
                                          'periodogram')){

  if(type == 'age_depth'){
    age_depth_plot(age_model)
  } else if(type == 'sed_rate') {
    plot_sed_rate(age_model)
  } else if(type == 'trace') {
    plot_trace(age_model)
  } else if(type == 'periodogram') {
    plot_pgram(age_model)
  }
}

###############################################################################
age_depth_plot <- function(age_model) {
  # make the age_depth plot ---------------------------------------------------
  ggplot(data = age_model$geochron_data,
         mapping = aes(x = position,
                       y = age,
                       color = id)) +
    geom_point(size = 2) +
    geom_errorbar(mapping = aes(ymin = age - age_sd * 2,
                                ymax = age + age_sd * 2),
                  width = 0.75) +
    xlab('Depth') +
    ylab('Age') +
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
    theme(legend.position = 'top') +
    scale_color_viridis(discrete = TRUE, option = 'H')
}
###############################################################################
plot_sed_rate <- function(age_model) {
  plots <- list()
  colors <- viridis(n = ncol(age_model$sed_rate), option = 'D')
  for(k in 1:ncol(age_model$sed_rate)){
    dat <- data.frame(x = age_model$burn:age_model$iterations,
                      y = age_model$sed_rate[age_model$burn:age_model$iterations, k])

    plots[[k]] <- dat %>%
      ggplot(mapping = aes(x = y)) +
      geom_density(fill = colors[k],
                   alpha = 0.75,
                   color = NA) +
      xlim(range(age_model$sed_rate[age_model$burn:age_model$iterations, ])) +
      xlab('sed rate (m/Ma)') +
      theme_bw() +
      ggtitle(label = paste(age_model$segment_edges$position[k],
                            '-',
                            age_model$segment_edges$position[k + 1],
                            'meters',
                            '; median sed rate = ',
                            round(mean(age_model$sed_rate[, k],
                                       na.rm = TRUE), 2),
                            'm/Ma'))
  }
  cowplot::plot_grid(plotlist = plots,
                     nrow = length(plots))

}

###############################################################################
plot_trace <- function(age_model) {
  plots <- list()
  colors <- viridis(n = ncol(age_model$sed_rate), option = 'D')
  for(k in 1:ncol(age_model$sed_rate)){
    dat <- data.frame(x = 1:age_model$iterations,
                      y = age_model$sed_rate[, k])

    plots[[k]] <- ggplot(data = dat,
                         mapping = aes(x = x,
                                       y = y)) +
      geom_line(color = colors[k]) +
      xlab('iteration') +
      ylab('sed rate') +
      theme_bw() +
      ggtitle(label = paste(age_model$segment_edges$position[k],
                            '-',
                            age_model$segment_edges$position[k + 1],
                            'meters',
                            '; median sed rate = ',
                            round(mean(age_model$sed_rate[, k],
                                       na.rm = TRUE), 2),
                            'm/Ma'))
  }
  p <- cowplot::plot_grid(plotlist = plots,
                          nrow = length(plots))
  p
}

###############################################################################
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
             time_freq   = Frequency * quantile(age_model$sed_rate[age_model$burn:age_model$iterations, k],
                                                prob = 0.5,
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
                            '; median sed rate = ',
                            round(mean(age_model$sed_rate[, k],
                                       na.rm = TRUE), 2),
                            'm/Ma')) +
      theme_bw() +
      theme(legend.position = 'none',
            axis.text.y = element_blank()) +
      xlab('frequency (m/Ma)') +
      theme(legend.position = 'none')
  }
  cowplot::plot_grid(plotlist = cowplot::align_plots(plotlist = plots))
}
