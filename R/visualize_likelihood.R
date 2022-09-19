#' helper function to visualize the likeihood functions given a model segementaiton scheme

#'
#' @name vizualize_likelihood
#'
#' @import "tidyverse"
#' @import "dplyr"
#' @import "astrochron"
#' @import "tibble"
#' @importFrom magrittr "%>%"
#' @md
#' @export
#'
vizualize_likelihood <- function(sed_prior_range = c(0.1, 100),
                                 segment_edges,
                                 cyclostrat_data,
                                 tuning_frequency,
                                 resolution = 0.1,
                                 method = 'malinverno') {

  # calculate the number of segments
  n <- nrow(segment_edges) - 1

  # give the segments some names
  segment_names <- vector(length = n)
  segment_number <- 1:n
  for(z in 1:n) {
    segment_names[z] <- paste(segment_edges$position[z],
                              '-',
                              segment_edges$position[z + 1],
                              'meters')
  }

  names(segment_number) = segment_names

  # calculate a grid of sedimentation rates
  sed_grid <- seq(
    min(sed_prior_range),
    max(sed_prior_range),
    by = resolution
  )

  likelihood <- matrix(ncol = n,
                       nrow = length(sed_grid))

  for(i in 1:n) {
    for(j in seq_along(sed_grid)) {
      likelihood[j, i] <- calculate_likelihood(cyclostrat_data = cyclostrat_data,
                                               tuning_frequency = tuning_frequency,
                                               sed_rate = sed_grid[j],
                                               segment_edges = segment_edges$position[i:(i + 1)],
                                               method = method)  %>% exp()
    }
  }

  colnames(likelihood) <- segment_names

  likelihood <- as.data.frame(likelihood) %>%
    add_column(sed_rate = sed_grid) %>%
    pivot_longer(cols = -sed_rate,
                 values_to = 'probability')

  likelihood$name <- factor(likelihood$name, levels = segment_names)


  likelihood %>%
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
    xlab('') +
    ggtitle(paste(method, 'likelihood'))

}
