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
#' @import "ggridges"
#'
#' @return a list with a whole bunch of stuff in it.
#' @md
#' @export
#'

plot.astroBayesModel <- function(age_model,
                                 type = c('age_depth',
                                          'sed_rate',
                                          'trace',
                                          'periodogram',
                                          'cyclostrat')) {

  switch(type,
         'age_depth' = age_depth_plot(age_model),
         'sed_rate' = plot_sed_rate(age_model),
         'trace' = plot_trace(age_model),
         'periodogram' = plot_pgram(age_model),
         'cyclostrat' = cyclostrat_plot(age_model))
}

###############################################################################
age_depth_plot <- function(age_model) {
  # assemble ridges for plotting
  n = 2000
  ridges <- age_model$geochron_data %>%
    mutate(low = age - age_sd * 4, high = age + age_sd * 4)  %>%
    uncount(n, .id = 'row') %>%
    mutate(x = (1 - row/n) * low + row/n*high,
           density = dnorm(x, age, age_sd))

  # make the plot
  ridges %>%
    ggplot(mapping = aes(x = x,
                         y = position,
                         height = density,
                         group = id,
                         fill = id)) +
    geom_density_ridges(stat = 'identity',
                        scale = 0.25,
                        color = NA,
                        alpha = 0.75) +
    ylab('Depth') +
    xlab('Age') +
    scale_y_reverse() +
    geom_ribbon(data = age_model$CI,
                mapping = aes(xmin = CI_2.5,
                              xmax = CI_97.5,
                              y = position),
                inherit.aes = FALSE,
                alpha = 0.25) +
    geom_line(data = age_model$CI,
              mapping = aes(y = position,
                            x = median),
              inherit.aes = FALSE) +
    theme_bw() +
    theme(legend.position = 'top') +
    scale_fill_viridis(discrete = TRUE, option = 'plasma', end = 0.9)

}
###############################################################################
plot_sed_rate <- function(age_model) {
  plots <- list()
  colors <- viridis(n = ncol(age_model$sed_rate), option = 'D',  end = 0.9)
  for(k in 1:ncol(age_model$sed_rate)){
    dat <- data.frame(x = age_model$burn:age_model$iterations,
                      y = age_model$sed_rate[age_model$burn:age_model$iterations, k])

    plots[[k]] <- dat %>%
      ggplot(mapping = aes(x = y)) +
      geom_density(fill = colors[k],
                   alpha = 0.75,
                   color = NA) +
      xlim(range(age_model$sed_rate[age_model$burn:age_model$iterations, k])) +
      xlab('sed rate (m/Ma)') +
      theme_bw() +
      ggtitle(label = paste(age_model$segment_edges$position[k],
                            '-',
                            age_model$segment_edges$position[k + 1],
                            'meters \n',
                            'median sed rate = \n',
                            round(quantile(age_model$sed_rate[age_model$burn:age_model$iterations, k],
                                           na.rm = TRUE,
                                           prob = 0.5), 2),
                            'm/Ma')) +
      geom_vline(xintercept = quantile(age_model$sed_rate[age_model$burn:age_model$iterations, k],
                                       na.rm = TRUE,
                                       prob = 0.5),
                 linetype = 'dashed',
                 color = 'red')
  }
  cowplot::plot_grid(plotlist = plots)

}

###############################################################################
plot_trace <- function(age_model) {
  plots <- list()
  colors <- viridis(n = ncol(age_model$sed_rate), option = 'D', end = 0.9)
  for(k in 1:ncol(age_model$sed_rate)){
    dat <- data.frame(x = 1:age_model$iterations,
                      y = age_model$sed_rate[, k])

    plots[[k]] <- ggplot(data = dat,
                         mapping = aes(x = x,
                                       y = y)) +
      geom_line(color = colors[k], alpha = 0.5) +
      xlab('iteration') +
      ylab('sed rate') +
      theme_bw() +
      ggtitle(label = paste(age_model$segment_edges$position[k],
                            '-',
                            age_model$segment_edges$position[k + 1],
                            'meters \n',
                            'median sed rate = \n',
                            round(quantile(age_model$sed_rate[age_model$burn:age_model$iterations, k],
                                           na.rm = TRUE,
                                           prob = 0.5),
                                  2),
                            'm/Ma')) +
      geom_hline(yintercept = quantile(age_model$sed_rate[age_model$burn:age_model$iterations, k],
                                       na.rm = TRUE,
                                       prob = 0.5),
                 linetype = 'dashed',
                 color = 'red')
  }
  cowplot::plot_grid(plotlist = plots)
}

###############################################################################
plot_pgram <- function(age_model) {

  # preallocate a list to store plots
  plots <- list()

  # preallocate colors for different segments
  colors <- viridis(n = length(age_model$cyclostrat_data),
                    option = 'plasma',
                    end = 0.75)


  for(i in  seq_along(age_model$cyclostrat_data)) {
    plots[[i]] <- age_model$cyclostrat_CI[[i]] |>
      select(median,
             value) |>
      linterp(genplot = FALSE,
              verbose = FALSE) |>
      periodogram(genplot = FALSE,
                  verbose = FALSE,
                  background = 1,
                  f0 = TRUE,
                  output = 1,
                  padfac = 10) |>
      ggplot(mapping = aes(x = Frequency,
                           y = Power)) +
      geom_line(color = colors[i]) +
      geom_line(mapping = aes(y = AR1_Fit),
                color = 2,
                linetype = 'solid') +
      geom_line(mapping = aes(y = AR1_95_power),
                color = 2,
                linetype = 'dashed') +
      theme_bw() +
      geom_vline(data = age_model$tuning_frequency,
                 mapping = aes(xintercept = frequency),
                 color = 'grey',
                 linetype = 'dashed',
                 size = 0.5) +
      xlab('Frequency (cycles/Ma)') +
      xlim(0, max(age_model$tuning_frequency$frequency))
      ylab('Spectral Power') +
      ggtitle(LETTERS[i]) +
      theme(panel.grid = element_blank())
  }

  cowplot::plot_grid(plotlist = cowplot::align_plots(plotlist = plots),
                     ncol = 1)

}

###############################################################################
cyclostrat_plot <- function(age_model) {

  # assign identifier letters
  for(i in seq_along(age_model$cyclostrat_CI)) {
    age_model$cyclostrat_CI[[i]] <-
      age_model$cyclostrat_CI[[i]] %>%
      add_column(id = LETTERS[i])
  }

  # make the plot
  p <- age_model$cyclostrat_CI %>%
    reduce(rbind) %>%
    ggplot(mapping = aes(x = median,
                         y = value,
                         color = id)) +
    geom_line(size = 0.25) +
    geom_point(size = 0.5) +
    facet_grid(id~.) +
    theme_bw() +
    xlab('Age (Ma)') +
    theme(legend.position = 'none') +
    scale_color_viridis(option = 'plasma',
                        discrete = TRUE,
                        end = 0.75)
  return(p)
}
