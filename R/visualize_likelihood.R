#' helper function to visualize the likeihood functions given a model segementaiton scheme

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
visualize_likelihood <- function(segment_edges,
                                 cyclostrat_data,
                                 tuning_frequency,
                                 resolution = 0.1,
                                 method = 'malinverno',
                                 plot = TRUE) {

  # calculate the number of segments
  n_segments     <- nrow(segment_edges) - 1

  # calculate the number of cyclostrat_records
  n_cyclostrat   <- length(cyclostrat_data)

  # give the segments some names
  segment_names  <- vector(length = n_segments)
  segment_number <- 1:n_segments

  for(z in 1:n_segments) {
    segment_names[z] <- paste(segment_edges$position[z],
                              '-',
                              segment_edges$position[z + 1],
                              'meters')
  }

  # assign the names
  names(segment_number) = segment_names

  # calculate a grid of sedimentation rates
  sed_grid <- seq(
    min(segment_edges$sed_min, na.rm = TRUE),
    max(segment_edges$sed_max, na.rm = TRUE),
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
                                      segment_edges = segment_edges$position[i:(i + 1)],
                                      cyclostrat = cyclostrat_data[[k]],
                                      tuning_frequency = tuning_frequency,
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
      ggtitle(paste(method, 'likelihood'))
    print(p)
  }

  return(likelihood)
}
