#' Helper function to visualize the likelihood functions given a model layering scheme.
#'
#' @description This function *does not* run the full Bayesian MCMC `astro_bayes_model`. Instead it
#' calculate only the astronomical likelihood of the `cyclostrat_data` across a range of
#' sedimentation rates given the `target_frequency` data frame. It is intended to help check for uni modal
#' (or close to uni modal) likelihoods so layer boundary positions can be adjusted as needed.
#'
#' @name visualize_likelihood
#'
#' @import "tidyverse"
#' @import "dplyr"
#' @import "astrochron"
#' @import "tibble"
#' @importFrom magrittr "%>%"
#' @md
#' @export
#'
visualize_likelihood <- function(layer_boundaries,
                                 cyclostrat_data,
                                 target_frequency,
                                 resolution = 0.1,
                                 method = 'malinverno',
                                 plot = TRUE) {

  # make sure cyclostrat_data is a list
  if(class(cyclostrat_data) == 'data.frame') {
    cycle_list <- list()
    cycle_list[[1]] <- cyclostrat_data
    cyclostrat_data <- cycle_list
    rm(cycle_list)
  }

  # calculate the number of layers
  n_segments     <- nrow(layer_boundaries) - 1

  # calculate the number of cyclostrat_records
  n_cyclostrat   <- length(cyclostrat_data)

  # give the segments some names
  segment_names  <- vector(length = n_segments)
  segment_number <- 1:n_segments

  for(z in 1:n_segments) {
    segment_names[z] <- paste(layer_boundaries$position[z],
                              '-',
                              layer_boundaries$position[z + 1],
                              'meters')
  }

  # assign the names
  names(segment_number) = segment_names

  # calculate a grid of sedimentation rates
  sed_grid <- seq(
    min(layer_boundaries$sed_min, na.rm = TRUE),
    max(layer_boundaries$sed_max, na.rm = TRUE),
    by = resolution)

  # set up a matrix for storage
  likelihood <- matrix(ncol = n_segments,
                       nrow = length(sed_grid))

  # loop through the sed grid and calculate each likelihood
  for(i in 1:n_segments) {
    for(j in seq_along(sed_grid)) {
      LL <- vector()
      for(k in seq_along(cyclostrat_data)) {
        LL[k] <- calculate_likelihood(sed_rate = sed_grid[j],
                                      layer_boundaries = layer_boundaries$position[i:(i + 1)],
                                      cyclostrat = cyclostrat_data[[k]],
                                      target_frequency = target_frequency,
                                      method = method) |> exp()
      }

      likelihood[j, i] <- sum(LL)
    }
  }

  # give the matrix names
  colnames(likelihood) <- segment_names

  # reshape into a ggplot friendly data_frame
  likelihood <- as.data.frame(likelihood) %>%
    add_column(sed_rate = sed_grid) %>%
    pivot_longer(cols = -sed_rate,
                 values_to = 'probability')

  likelihood$name <- factor(likelihood$name,
                            levels = segment_names)

  if(plot) {
    # plot the results
    p <- likelihood %>%
      ggplot(mapping = aes(x = sed_rate,
                           y = probability,
                           color = name,
                           fill = name)) +
      geom_area() +
      facet_wrap(~name,
                 scales = 'free_y',
                 ncol = 1,
                 strip.position = 'right') +
      theme_bw() +
      theme(axis.text.y = element_blank(),
            legend.position = 'none')  +
      xlab('sedimentation rate (m/Ma)') +
      ggtitle(paste(method, 'likelihood')) +
      scale_fill_viridis(discrete = TRUE,
                         option = 'plasma',
                         end = 0.9) +
      scale_color_viridis(discrete = TRUE,
                         option = 'plasma',
                         end = 0.9)
    print(p)
  }

  return(likelihood)
}
